!=====================================================================
!
!     SUBROUTINE START_CHEM
!
!        This subroutine is a driver subroutine for gmp_start.F90.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE START_CHEM

      USE COMODE_MOD,       ONLY : INIT_COMODE
      USE DAO_MOD,          ONLY : INIT_DAO
      USE DIAG_AOD_MOD,     ONLY : INIT_DIAG_AOD
      USE ERROR_MOD,        ONLY : DEBUG_MSG, ERROR_STOP
      USE FLEXAOD,          ONLY : READ_MIE
      USE GCKPP_COMODE_MOD, ONLY : INIT_GCKPP_COMODE
      USE INPUT_MOD,        ONLY : READ_INPUT_FILE
      USE LINOZ_MOD,        ONLY : LINOZ_READ
      USE LOGICAL_MOD,      ONLY : LCHEM, LEMIS, LLINOZ, LMEGAN, LPRT
      USE LOGICAL_MOD,      ONLY : LSVMET, LSVMETNEST, LKPP
      USE LOGICAL_MOD,      ONLY : LAERODIR, LAEROIND
      USE MEGAN_MOD,        ONLY : INIT_MEGAN
      USE METEO_MOD,        ONLY : INIT_METEO
      USE PRESSURE_MOD,     ONLY : INIT_PRESSURE
      USE RESTART_MOD,      ONLY : READ_RESTART

      USE CMN_SIZE_MOD    ! Size parameters   
      USE COMODE_LOOP_MOD ! SMVGEAR common blocks

      IMPLICIT NONE

      INTEGER :: RC

      RC=0 ! Error flag. Default 0 is "no error"

      ! Initialize header modules
      CALL INIT_HEADER
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_HEADER' )

      ! Read input file and call init routines from other modules
      CALL READ_INPUT_FILE
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a READ_INPUT_FILE' )

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

      ! Initialize allocatable SMVGEAR/KPP arrays
      IF ( LEMIS .or. LCHEM ) THEN
         CALL INIT_COMODE
         IF ( LKPP )  CALL INIT_GCKPP_COMODE( IIPAR, JJPAR, LLTROP,
     $        ITLOOP, NMTRATE, IGAS, RC )
         IF ( RC == 1 )
     $        CALL ERROR_STOP( "Alloc error", "INIT_GCKPP_COMODE" )
         IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_COMODE' )
      ENDIF

      ! Added to read input file for linoz strat (dbj, jliu, bmy, 10/16/09)
      IF ( LLINOZ ) CALL LINOZ_READ

      ! Initialize MEGAN arrays, get 15-day avg temperatures
      IF ( LMEGAN ) THEN
         CALL INIT_MEGAN
         CALL INITIALIZE( 2 )
         CALL INITIALIZE( 3 )
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INIT_MEGAN' )
      ENDIF

      ! Read restart file from the disk
      CALL READ_RESTART
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a READ_RESTART' )

      ! Read MIE tables for aerosol direct effect
      IF ( LAERODIR ) THEN
         CALL READ_MIE
         IF ( LPRT ) CALL DEBUG_MSG( '### START: a READ_MIE' )
      ENDIF

      ! Initialize METEO module
      IF ( LSVMET .or. LSVMETNEST ) THEN
         CALL INIT_METEO
         IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_METEO' )
      ENDIF

      ! Initialize AOD diagnostic module
      IF ( LAERODIR .or. LAEROIND ) THEN
         CALL INIT_DIAG_AOD
         IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_DIAG_AOD' )
      ENDIF

      END SUBROUTINE START_CHEM
