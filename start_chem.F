!=====================================================================
!
!     SUBROUTINE START_CHEM
!
!        This subroutine is a driver subroutine for gmp_start.F90.
!
!                                            written by Seungun Lee
!
!=====================================================================

      SUBROUTINE START_CHEM

      USE COMCHEM,        ONLY : INIT_COMCHEM
      USE COMODE_MOD,     ONLY : INIT_COMODE
      USE DAO_MOD,        ONLY : INIT_DAO
      USE DIAG_AOD_MOD,   ONLY : INIT_DIAG_AOD
      USE ERROR_MOD,      ONLY : DEBUG_MSG, ERROR_STOP
      USE FLEXAOD,        ONLY : MIE_CALC
      USE INPUT_MOD,      ONLY : READ_INPUT_FILE
      USE LOGICAL_MOD,    ONLY : LCHEM, LEMIS, LPRT, LVARTROP
      USE PRESSURE_MOD,   ONLY : INIT_PRESSURE
      USE RESTART,        ONLY : READ_RESTART
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TROPOPAUSE_MOD, ONLY : READ_TROPOPAUSE

      IMPLICIT NONE

      INTEGER :: RC

      RC=0 ! Error flag. Default 0 is "no error"

      ! Initialize header modules
      CALL INIT_HEADER
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_HEADER' )

      ! Read input file and call init routines from other modules
      CALL READ_INPUT_FILE
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a READ_INPUT_FILE' )

      ! Initialize GRIMs coupling module
      CALL INIT_COMCHEM
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_COMCHEM' )

      ! Initialize AOD diagnostic module
      CALL INIT_DIAG_AOD
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_DIAG_AOD' )

      ! Initialize met field arrays from "dao_mod.f"
      CALL INIT_DAO
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_DAO' )

      ! Initialize diagnostic arrays and counters
      CALL INITIALIZE( 2 )
      CALL INITIALIZE( 3 )
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INITIALIZE' )

      ! Initialize the new hybrid pressure module.  Define Ap and Bp.
      CALL INIT_PRESSURE
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_PRESSURE' )

      ! Read annual mean tropopause if not a variable tropopause
      ! read_tropopause is obsolete with variable tropopause
      IF ( .not. LVARTROP ) THEN
         CALL READ_TROPOPAUSE
         IF ( LPRT ) CALL DEBUG_MSG( '### START: a READ_TROPOPAUSE' )
      ENDIF

      ! Initialize allocatable SMVGEAR/KPP arrays
      IF ( LEMIS .or. LCHEM ) THEN
         IF ( ITS_A_FULLCHEM_SIM() ) CALL INIT_COMODE
         IF ( ITS_AN_AEROSOL_SIM() ) CALL INIT_COMODE
!         IF ( LKPP )  CALL INIT_GCKPP_COMODE( IIPAR, JJPAR, LLTROP,
!     $        ITLOOP, NMTRATE, IGAS, RC )
         IF ( RC == 1 )
     $        CALL ERROR_STOP( "Alloc error", "INIT_GCKPP_COMODE" )
         IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_COMODE' )
      ENDIF

      ! Read land types and fractions from "vegtype.global"
      CALL RDLAND
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a RDLAND' )

      CALL READ_RESTART
      CALL MIE_CALC

      END SUBROUTINE START_CHEM
