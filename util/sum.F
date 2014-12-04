#include "define.h"
!=====================================================================
!
!     FUNCTION MAKE_SUM
!
!        This function returns the sum of input array.
!
!     FUNCTION MPISUM
!
!        This function gathers the value saved in each pes,
!        and then returns the sum of them.
!
!                                            Written by Seungun Lee
!
!=====================================================================
      FUNCTION MAKE_SUM( ARRAY, VERT ) RESULT( TOTAL )

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      INTEGER :: VERT
      REAL*8 :: ARRAY(IIPAR,JJPAR,VERT)
      REAL*8 :: F_ARRAY(IGLOB*2,JGLOB/2,VERT)
      REAL*8 :: TOTAL

#ifdef MP
      CALL MPGP2F(ARRAY,IIPAR,JJPAR,F_ARRAY,IGLOB*2,JGLOB/2,VERT)
#else
      F_ARRAY = ARRAY
#endif

      TOTAL = SUM( F_ARRAY )

      END FUNCTION MAKE_SUM

#ifdef MP
      SUBROUTINE MPISUM( TOTAL )

      USE COMMPI

      IMPLICIT NONE

      REAL*8 :: TOTAL
      REAL*8 :: BUFFER(NPES)
      INTEGER :: IERR

      CALL MPI_GATHER( TOTAL, 1, MPI_REAL8, BUFFER, 1, MPI_REAL8,
     &                 0, MPI_COMM_WORLD, IERR )

      TOTAL = SUM( BUFFER )

      END SUBROUTINE MPISUM
#endif
