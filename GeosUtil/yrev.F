!=====================================================================
!
!     SUBROUTINE YREV
!
!        This subroutine makes an array reversed in y-direction.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE YREV( ARRAY, NI, NJ, NL )

      IMPLICIT NONE

      INTEGER :: NI, NJ, NL
      INTEGER :: I, J, L
      REAL*8  :: ARRAY(NI,NJ,NL)
      REAL*8  :: TMP(NI,NJ,NL)

      ! Reverse
      DO L = 1, NL
      DO J = 1, NJ
      DO I = 1, NI
         TMP(I,J,L) = ARRAY(I,NJ+1-J,L)
      ENDDO
      ENDDO
      ENDDO

      ! Return
      ARRAY = TMP

      END SUBROUTINE YREV
