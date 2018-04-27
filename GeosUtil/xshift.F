!=====================================================================
!
!     SUBROUTINE XSHIFT
!
!        This subroutine shifts an array in x-direction by a half of
!        the size of x-dimension.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE XSHIFT( ARRAY, NI, NJ, NL )

      IMPLICIT NONE

      INTEGER :: NI, NJ, NL
      INTEGER :: I, J, L
      REAL*8  :: ARRAY(NI,NJ,NL)
      REAL*8  :: TMP(NI/2)

      ! Shift
      DO L = 1, NL
      DO J = 1, NJ
         DO I = 1, NI/2
            TMP(I) = ARRAY(I,J,L)
         ENDDO
         DO I = 1, NI/2
            ARRAY(I,J,L) = ARRAY(NI/2+I,J,L)
         ENDDO
         DO I = 1, NI/2
            ARRAY(NI/2+I,J,L) = TMP(I)
         ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE XSHIFT
