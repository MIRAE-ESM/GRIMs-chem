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

      USE FILE_MOD, ONLY : IU_BPCH, IU_DEBUG, IU_ND48, IU_SMV2LOG
      USE TIME_MOD, ONLY : SET_CURRENT_TIME, PRINT_CURRENT_TIME
      USE TIME_MOD, ONLY : GET_TAU, SET_DIAGe, SET_DIAGb
      USE RESTART,  ONLY : WRITE_RESTART

      IMPLICIT NONE

      REAL*8 :: TAU

      ! Compute & print time quantities at start of dyn step
      CALL SET_CURRENT_TIME
      CALL PRINT_CURRENT_TIME

      TAU = GET_TAU()

      !==============================================================
      !   ***** W R I T E   D I A G N O S T I C   F I L E S *****
      !==============================================================

      ! Set time at end of diagnostic timestep
      CALL SET_DIAGe( TAU )

      ! Write bpch file
      CALL DIAG3

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

      CALL WRITE_RESTART

      END SUBROUTINE END_CHEM
