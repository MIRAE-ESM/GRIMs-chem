#include "define.h"
!=====================================================================
!
!     FUNCTION SUM_ARRAY
!
!        This function returns the sum of input array.
!
!                                                       Seungun Lee
!
!=====================================================================

      FUNCTION SUM_ARRAY( ARRAY, NL ) RESULT( TOTAL )

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      INTEGER :: NL
      REAL*8  :: ARRAY(IIPAR,JJPAR,NL)
      REAL*8  :: TOTAL
      REAL*8  :: SUM_VALUE

      ! Make sum at each node and then gather the values
      TOTAL = SUM_VALUE( SUM( ARRAY ) )

      END FUNCTION SUM_ARRAY

!=====================================================================
!
!     FUNCTION SUM_VALUE
!
!        This function returns the sum of input value.
!
!                                                       Seungun Lee
!
!=====================================================================

      FUNCTION SUM_VALUE( VALUE ) RESULT( TOTAL )

#ifdef MP
      USE commpi
#endif

      IMPLICIT NONE

      REAL*8  :: VALUE, TOTAL, L_TOT
#ifdef MP
      INTEGER :: IERR
      REAL*8  :: BUFFER(npes)
#endif

#ifdef MP
      ! Gather values from all PE
      CALL MPI_GATHER( VALUE, 1, MPI_REAL8, BUFFER, 1, MPI_REAL8,
     &                 0, comm_world, IERR )

      ! Sum
      TOTAL = SUM( BUFFER )
#else
      ! Sum
      TOTAL = VALUE
#endif

      ! Remove floating-point error
      IF ( TOTAL .GT. 0. ) THEN
         L_TOT = LOG(TOTAL)
         L_TOT = REAL(L_TOT-INT(L_TOT),4) + INT(L_TOT)
         TOTAL = EXP(L_TOT)
      ENDIF

      END FUNCTION SUM_VALUE
