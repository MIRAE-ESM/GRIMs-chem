#include "define.h"
!=====================================================================
!
!     SUBROUTINE TRANS_INDEX
!
!        This subroutine performs regridding and transforming of
!        index data by matching nearest point.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE TRANS_INDEX( INDATA, OUTDATA, IN_NI, IN_NJ )

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      INTEGER :: IN_NI, IN_NJ
      INTEGER :: INDATA(IN_NI,IN_NJ)
      INTEGER :: OUTDATA(IIPAR,JJPAR)
      INTEGER :: I, J, XID, YID
      INTEGER :: NEAR_X, NEAR_Y

#ifdef MP
      ! Broadcast
#ifndef RMP
      CALL mpbcasti(INDATA,IN_NI*IN_NJ)
#else
      CALL rmpbcasti(INDATA,IN_NI*IN_NJ)
#endif
#endif

      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Get nearest point
         XID = NEAR_X(I,J)
         YID = NEAR_Y(I,J)

         ! Get data
         OUTDATA(I,J) = INDATA(XID,YID)

      ENDDO
      ENDDO

      END SUBROUTINE TRANS_INDEX

!=====================================================================
!
!     FUNCTION NEAR_X
!
!        This function returns the nearest GRIMs x-index from
!        GMAO 2x2.5 grids.
!
!                                                       Seungun Lee
!
!=====================================================================

      FUNCTION NEAR_X( X, Y ) RESULT( XID )

      USE GRID_MOD, ONLY : GET_XMID
      USE CMN_SIZE_MOD

      IMPLICIT NONE

      INTEGER :: X, Y, XID
      REAL*8  :: LON

      ! Get longitude
      LON = GET_XMID(X,Y)

      ! Get nearest XID
      XID = NINT( ( LON + 180d0 ) / 2.5d0 ) + 1

      ! Cyclic longitude
      IF ( XID .GT. I2x25 ) XID = XID - I2x25

      END FUNCTION NEAR_X

!=====================================================================
!
!     FUNCTION NEAR_Y
!
!        This function returns the nearest GRIMs y-index from
!        GMAO 2x2.5 grids.
!
!                                                       Seungun Lee
!
!=====================================================================

      FUNCTION NEAR_Y( X, Y ) RESULT( YID )

      USE GRID_MOD, ONLY : GET_YMID

      IMPLICIT NONE

      INTEGER :: X, Y, YID
      REAL*8  :: LAT

      ! Get latitude
      LAT = GET_YMID(X,Y)

      ! Get nearest YID
      YID = NINT( ( LAT + 90d0 ) / 2.0d0 ) + 1

      END FUNCTION NEAR_Y
