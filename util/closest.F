!=====================================================================
!
!     FUNCTION CLX
!
!        This function returns the closest GRIMs x-index from
!        GMAO 2x2.5 grids.
!
!     FUNCTION CLY
!
!        This function returns the closest GRIMs y-index from
!        GMAO 2x2.5 grids.
!
!                                            written by Seungun Lee
!
!=====================================================================

      FUNCTION CLX( X, Y ) RESULT( CLX )

      USE GRID_MOD, ONLY : GET_XMID

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: X
      INTEGER, INTENT(IN)  :: Y
      INTEGER, INTENT(OUT) :: CLX
      INTEGER :: I2x25=144
      INTEGER :: I
      REAL*8 :: LON1, LON2, DIFF, DIFF_SAVE

      LON1 = GET_XMID(X,Y)

      DO I = 1, I2x25

         LON2 = -180. + DBLE(I-1) * 2.5

         DIFF = ABS(LON1-LON2)

         IF (I .EQ. 1) DIFF_SAVE = DIFF

         IF (DIFF .LE. DIFF_SAVE) THEN
            DIFF_SAVE = DIFF
            CLX = I
         ENDIF

      ENDDO

      END FUNCTION CLX

      FUNCTION CLY( X, Y ) RESULT( CLY )

      USE GRID_MOD, ONLY : GET_YMID

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: X
      INTEGER, INTENT(IN)  :: Y
      INTEGER, INTENT(OUT) :: CLY
      INTEGER :: J2x25=91
      INTEGER :: J
      REAL*8 :: LAT1, LAT2, DIFF, DIFF_SAVE

      LAT1 = GET_YMID(X,Y)

      DO J = 1, J2x25

         LAT2 = -90. + DBLE(J-1) * 2.

         DIFF = ABS(LAT1-LAT2)

         IF (J .EQ. 1) DIFF_SAVE = DIFF

         IF (DIFF .LE. DIFF_SAVE) THEN
            DIFF_SAVE = DIFF
            CLY = J
         ENDIF

      ENDDO

      END FUNCTION CLY
