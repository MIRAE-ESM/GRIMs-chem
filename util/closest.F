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
!     FUNCTION CLZ_GEOS4
!
!        This function returns the closest GRIMs z-index from
!        GEOS-4 vertical grids.
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

      FUNCTION CLZ_GEOS4( X, Y, Z ) RESULT( CLZ )

      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_PEDGE
      USE CMN_SIZE_MOD, ONLY : LLPAR

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: X
      INTEGER, INTENT(IN)  :: Y
      INTEGER, INTENT(IN)  :: Z
      INTEGER, INTENT(OUT) :: CLZ
      INTEGER :: L
      REAL*8 :: AP(56), BP(56)
      REAL*8 :: PEDGEDN, PEDGEUP, PCENTER1, PCENTER2, DIFF, DIFF_SAVE

      !-----------------------------------------------------------------
      ! GEOS-4 55 level grid
      !-----------------------------------------------------------------
      
      ! AP [hPa] for 55 levels (56 edges)  
      AP = (/ 0.000000d0,   0.000000d0,  12.704939d0,  35.465965d0, 
     &       66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0,
     &      198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0,
     &      216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0,
     &      127.837006d0, 108.663429d0,  92.365662d0,  78.512299d0, 
     &       66.603378d0,  56.387939d0,  47.643932d0,  40.175419d0,
     &       33.809956d0,  28.367815d0,  23.730362d0,  19.791553d0,
     &       16.457071d0,  13.643393d0,  11.276889d0,   9.292943d0,
     &        7.619839d0,   6.216800d0,   5.046805d0,   4.076567d0,
     &        3.276433d0,   2.620212d0,   2.084972d0,   1.650792d0,
     &        1.300508d0,   1.019442d0,   0.795134d0,   0.616779d0,
     &        0.475806d0,   0.365041d0,   0.278526d0,   0.211349d0,
     &        0.159495d0,   0.119703d0,   0.089345d0,   0.066000d0,
     &        0.047585d0,   0.032700d0,   0.020000d0,   0.010000d0 /)
     
      ! BP [unitless] for 55 levels (56 edges)
      BP = (/  1.000000d0,  0.985110d0,   0.943290d0,   0.867830d0,
     &         0.764920d0,  0.642710d0,   0.510460d0,   0.378440d0,
     &         0.270330d0,  0.183300d0,   0.115030d0,   0.063720d0,
     &         0.028010d0,  0.006960d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0 /)

      PEDGEDN = AP(Z)   + ( BP(Z)   * GET_PEDGE(X,Y,1) )
      PEDGEUP = AP(Z+1) + ( BP(Z+1) * GET_PEDGE(X,Y,1) )

      PCENTER1 = 0.5d0 * ( PEDGEDN + PEDGEUP )

      DO L = 1, LLPAR

         PCENTER2 = GET_PCENTER( X, Y, L )

         DIFF = ABS(PCENTER1-PCENTER2)

         IF (L .EQ. 1) DIFF_SAVE = DIFF

         IF (DIFF .LE. DIFF_SAVE) THEN
            DIFF_SAVE = DIFF
            CLZ = L
         ENDIF

      ENDDO

      END FUNCTION CLZ_GEOS4
