!=====================================================================
!
!     SUBROUTINE END_CHEM
!
!        This subroutine is a driver subroutine for gmp_end.F90.
!
!                                            written by Seungun Lee
!
!=====================================================================

      SUBROUTINE END_CHEM

      USE DIAG_AOD_MOD, ONLY : WRITE_DIAG_AOD
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE FILE_MOD,     ONLY : CLOSE_FILES
      USE FILE_MOD,     ONLY : IU_BPCH, IU_DEBUG, IU_ND48, IU_SMV2LOG
      USE LOGICAL_MOD,  ONLY : LPRT
      USE METEO_MOD,    ONLY : SAVE_A3, SAVE_I6, SAVE_A6
      USE TIME_MOD,     ONLY : SET_CURRENT_TIME, PRINT_CURRENT_TIME
      USE TIME_MOD,     ONLY : GET_TAU, SET_DIAGe, GET_NHMS
      USE RESTART_MOD,  ONLY : WRITE_RESTART

      IMPLICIT NONE

      INTEGER :: NHMS

      ! Compute & print time quantities at start of dyn step
      CALL SET_CURRENT_TIME
      CALL PRINT_CURRENT_TIME

      NHMS = GET_NHMS()

      !==============================================================
      !   ***** W R I T E   D I A G N O S T I C   F I L E S *****
      !==============================================================

      ! Set time at end of diagnostic timestep
      CALL SET_DIAGe( GET_TAU() )

      ! Write bpch file
      CALL DIAG3

      ! Write AOD, SSA, ASM
      CALL WRITE_DIAG_AOD

      ! Flush file units
      CALL FLUSH( IU_ND48    )
      CALL FLUSH( IU_BPCH    )
      CALL FLUSH( IU_SMV2LOG )
      CALL FLUSH( IU_DEBUG   )

      !===========================================================
      !        ***** Z E R O   D I A G N O S T I C S *****
      !===========================================================
      CALL INITIALIZE( 2 ) ! Zero arrays
      CALL INITIALIZE( 3 ) ! Zero counters

      ! Meteo
      IF ( MOD(NHMS,030000) .EQ. 0 ) THEN
         CALL SAVE_A3
      ENDIF
      IF ( MOD(NHMS,060000) .EQ. 0 ) THEN
         CALL SAVE_I6
         CALL SAVE_A6
      ENDIF

      ! Write restart file
      CALL WRITE_RESTART

      ! Close all files
      CALL CLOSE_FILES
      IF ( LPRT ) CALL DEBUG_MSG( '### END: a CLOSE_FILES' )

      ! Deallocate dynamic module arrays
      CALL CLEANUP
      IF ( LPRT ) CALL DEBUG_MSG( '### END: a CLEANUP' )

      END SUBROUTINE END_CHEM
