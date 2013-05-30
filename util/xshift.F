!=====================================================================
!
!     SUBROUTINE XSHIFT
!
!        This subroutine shifts an array in x-direction by a half of
!        the size of x-dimension.
!
!                                            written by Seungun Lee
!
!=====================================================================

      SUBROUTINE XSHIFT( IODATA, XSIZE, YSIZE, ZSIZE )

      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: XSIZE, YSIZE, ZSIZE
      REAL*8, INTENT(INOUT) :: IODATA(XSIZE,YSIZE,ZSIZE)
      REAL*8  :: ARRAY(XSIZE,YSIZE,ZSIZE)
      INTEGER :: I

      ARRAY = 0d0

      DO I=1, XSIZE/2
         ARRAY(I,:,:) = IODATA(I+XSIZE/2,:,:)
      ENDDO

      DO I=XSIZE/2+1, XSIZE
         ARRAY(I,:,:) = IODATA(I-XSIZE/2,:,:)
      ENDDO

      IODATA = ARRAY

      END SUBROUTINE XSHIFT
