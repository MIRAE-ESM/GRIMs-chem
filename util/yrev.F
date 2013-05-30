!=====================================================================
!
!     SUBROUTINE YREV
!
!        This subroutine makes an array reversed in y-direction.
!
!                                            written by Seungun Lee
!
!=====================================================================

      SUBROUTINE YREV( IODATA, XSIZE, YSIZE, ZSIZE )

      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: XSIZE, YSIZE, ZSIZE
      REAL*8, INTENT(INOUT) :: IODATA(XSIZE,YSIZE,ZSIZE)
      REAL*8  :: ARRAY(XSIZE,YSIZE,ZSIZE)
      INTEGER :: J

      ARRAY = 0d0

      DO J=1, YSIZE
         ARRAY(:,J,:) = IODATA(:,YSIZE-J+1,:)
      ENDDO

      IODATA = ARRAY

      END SUBROUTINE YREV
