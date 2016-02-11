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
      USE paramodel

      IMPLICIT NONE

      INTEGER :: NL
      INTEGER :: I, J, L
      REAL*8  :: ARRAY(IIPAR,JJPAR,NL)
      REAL*8  :: P_ARRAY(LONF2S,LATG2S,NL)
      REAL*8  :: F_ARRAY(lonf2_,latg2_,NL)
      REAL*8  :: TOTAL

      ! Copy
      DO L = 1, NL
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         P_ARRAY(I,J,L) = ARRAY(I,J,L)
      ENDDO
      ENDDO
      ENDDO

#ifdef MP
      ! MPI partial grid to full grid
      CALL mpgp2f(P_ARRAY,LONF2S,LATG2S,F_ARRAY,lonf2_,latg2_,NL)
#else
      F_ARRAY = P_ARRAY
#endif

      ! Sum
      TOTAL = SUM( F_ARRAY )

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

      REAL*8  :: VALUE, TOTAL
#ifdef MP
      INTEGER :: IERR
      REAL*8  :: BUFFER(npes)
#endif

#ifdef MP
      ! Gather values from all PE
      CALL MPI_GATHER( VALUE, 1, MPI_REAL8, BUFFER, 1, MPI_REAL8,
     &                 0, MPI_COMM_WORLD, IERR )

      ! Sum
      TOTAL = SUM( BUFFER )
#else
      ! Sum
      TOTAL = VALUE
#endif

      ! Remove numerical noise
      TOTAL = DBLE( REAL( TOTAL, 4 ) )

      END FUNCTION SUM_VALUE
