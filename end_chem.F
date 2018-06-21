!=====================================================================
!
!     SUBROUTINE END_CHEM
!
!        This subroutine is a driver subroutine for gmp_end.F90.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE END_CHEM

      USE DAO_MOD,      ONLY : AD, CONVERT_UNITS
      USE DIAG_AOD_MOD, ONLY : WRITE_DIAG_AOD
      USE DIAG_OH_MOD,  ONLY : PRINT_DIAG_OH
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE FILE_MOD,     ONLY : CLOSE_FILES
      USE FILE_MOD,     ONLY : IU_BPCH, IU_DEBUG, IU_ND48, IU_SMV2LOG
      USE LOGICAL_MOD,  ONLY : LPRT, LSVMET, LSVMETNEST
      USE LOGICAL_MOD,  ONLY : LAERODIR, LAEROIND
      USE METEO_MOD,    ONLY : INIT_A3, INIT_A6
      USE METEO_MOD,    ONLY : SAVE_A3, SAVE_A6, SAVE_I6
      USE METEO_MOD,    ONLY : SAVE_A3_NESTED, SAVE_A6_NESTED
      USE METEO_MOD,    ONLY : SAVE_I6_NESTED
      USE RESTART_MOD,  ONLY : WRITE_RESTART
      USE TIME_MOD,     ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,     ONLY : GET_MONTH, GET_TAU
      USE TIME_MOD,     ONLY : GET_ELAPSED_SEC, GET_DAY_OF_YEAR
      USE TIME_MOD,     ONLY : GET_SEASON, GET_YEAR
      USE TIME_MOD,     ONLY : ITS_TIME_FOR_BPCH
      USE TIME_MOD,     ONLY : ITS_TIME_FOR_UNIT
      USE TIME_MOD,     ONLY : SET_DIAGb, SET_DIAGe
      USE TIME_MOD,     ONLY : SET_CURRENT_TIME, PRINT_CURRENT_TIME
      USE TIME_MOD,     ONLY : GET_HOUR, GET_MINUTE
      USE TRACER_MOD,   ONLY : N_TRACERS, STT, TCVV

      IMPLICIT NONE

      INTEGER            :: YEAR,          MONTH, DAY_OF_YEAR
      INTEGER            :: SEASON,        NYMD,  NHMS
      INTEGER            :: ELAPSED_SEC
      INTEGER            :: ELAPSED_TODAY, HOUR,  MINUTE
      REAL*8             :: TAU

      ! Compute & print time quantities at start of dyn step
      CALL SET_CURRENT_TIME
      CALL PRINT_CURRENT_TIME

      ! Set time variables for dynamic loop
      DAY_OF_YEAR   = GET_DAY_OF_YEAR()
      ELAPSED_SEC   = GET_ELAPSED_SEC()
      MONTH         = GET_MONTH()
      NHMS          = GET_NHMS()
      NYMD          = GET_NYMD()
      HOUR          = GET_HOUR()
      MINUTE        = GET_MINUTE()
      TAU           = GET_TAU()
      YEAR          = GET_YEAR()
      SEASON        = GET_SEASON()
      ELAPSED_TODAY = ( HOUR * 3600 ) + ( MINUTE * 60 )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### END: a SET_CURRENT_TIME' )

      !==============================================================
      !   ***** U N I T   C O N V E R S I O N  ( v/v -> kg ) *****
      !==============================================================
      IF ( ITS_TIME_FOR_UNIT() ) THEN
         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### END: a CONVERT_UNITS:1' )
      ENDIF

      !==============================================================
      !   ***** W R I T E   D I A G N O S T I C   F I L E S *****
      !==============================================================

      ! Set time at end of diagnostic timestep
      CALL SET_DIAGe( TAU )

      ! Write bpch file
      CALL DIAG3

      ! Write AOD, SSA, ASM
      IF ( LAERODIR .or. LAEROIND ) THEN
         CALL WRITE_DIAG_AOD
      ENDIF

      ! Flush file units
      CALL FLUSH( IU_ND48    )
      CALL FLUSH( IU_BPCH    )
      CALL FLUSH( IU_SMV2LOG )
      CALL FLUSH( IU_DEBUG   )

      ! Set time at beginning of next diagnostic timestep
      CALL SET_DIAGb( TAU )

      !===========================================================
      !        ***** Z E R O   D I A G N O S T I C S *****
      !===========================================================
      CALL INITIALIZE( 2 ) ! Zero arrays
      CALL INITIALIZE( 3 ) ! Zero counters

      !==============================================================
      !   ***** U N I T   C O N V E R S I O N  ( kg -> v/v ) *****
      !==============================================================
      IF ( ITS_TIME_FOR_UNIT() ) THEN
         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### END: a CONVERT_UNITS:2' )
      ENDIF

      ! Save GRIMs meteorology for offline GRIMs/GEOS-Chem
      IF ( LSVMET .or. LSVMETNEST ) THEN

         ! Save A3 fields
         IF ( MOD(NHMS,010000) .EQ. 0 ) THEN
            IF ( LSVMET     ) CALL SAVE_A3
            IF ( LSVMETNEST ) CALL SAVE_A3_NESTED
            CALL INIT_A3
         ENDIF

         ! Save A6 fields
         IF ( MOD(NHMS,010000) .EQ. 0 ) THEN
            IF ( LSVMET     ) CALL SAVE_A6
            IF ( LSVMETNEST ) CALL SAVE_A6_NESTED
            CALL INIT_A6
         ENDIF

      ENDIF

      ! Write restart file
      CALL WRITE_RESTART

      ! Print the mass-weighted mean OH concentration (if applicable)
      CALL PRINT_DIAG_OH

      ! Close all files
      CALL CLOSE_FILES
      IF ( LPRT ) CALL DEBUG_MSG( '### END: a CLOSE_FILES' )

      ! Deallocate dynamic module arrays
      CALL CLEANUP

      IF ( LPRT ) CALL DEBUG_MSG( '### END: a CLEANUP' )

      END SUBROUTINE END_CHEM
