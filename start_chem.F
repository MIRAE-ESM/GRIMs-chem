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

      USE COMODE_MOD,     ONLY : INIT_COMODE
      USE DAO_MOD,        ONLY : INIT_DAO
      USE DIAG_AOD_MOD,   ONLY : INIT_DIAG_AOD
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE FLEXAOD,        ONLY : MIE_CALC
      USE INPUT_MOD,      ONLY : READ_INPUT_FILE
      USE LOGICAL_MOD,    ONLY : LCHEM, LEMIS, LPRT
      USE PRESSURE_MOD,   ONLY : INIT_PRESSURE
      USE RESTART_MOD,    ONLY : READ_RESTART

      IMPLICIT NONE

      ! Initialize header modules
      CALL INIT_HEADER
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_HEADER' )

      ! Initialize AOD diagnostic module
      CALL INIT_DIAG_AOD
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_DIAG_AOD' )

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
         IF ( LPRT ) CALL DEBUG_MSG( '### START: a INIT_COMODE' )
      ENDIF

      ! Read restart file from the disk
      CALL READ_RESTART
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a READ_RESTART' )

      ! MIE calculation for aerosol direct effect
      CALL MIE_CALC
      IF ( LPRT ) CALL DEBUG_MSG( '### START: a MIE_CALC' )

      END SUBROUTINE START_CHEM
