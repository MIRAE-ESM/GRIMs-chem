!=====================================================================
!
!     SUBROUTINE MAKE_CLDTOPS
!
!        This subroutine makes CLDTOPS met-field using CLDMAS.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE MAKE_CLDTOPS

      USE DAO_MOD, ONLY : CLDMAS, CLDTOPS

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      INTEGER :: I, J, L

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         CLDTOPS(I,J) = 1
         DO L = LLPAR, 1, -1
            IF ( CLDMAS(I,J,L) > 0d0 ) THEN
               CLDTOPS(I,J) = L + 1
               EXIT
            ENDIF
         ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE MAKE_CLDTOPS
