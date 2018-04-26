#include "define.h"
!=====================================================================
!
!     SUBROUTINE TRANS_INPUT
!
!        This subroutine performs horizontal and vertical
!        regridding and transforms an array from various
!        input grid to model grid.
!
!     VERT=0 : no vertical regrid
!     VERT=1 : z weighted vertical regrid
!     VERT=2 : p weighted vertical regrid
!     MASS=0 : conc regrid
!     MASS=1 : mass regrid
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE TRANS_INPUT_1x1(INDATA,OUTDATA,IN_NL,OUT_NL,VERT,MASS)

#ifndef RMP
      USE CMN_SIZE_MOD
      USE module_trans
      USE paramodel
#ifdef MP
      USE commpi
#endif

      IMPLICIT NONE

      INTEGER :: IN_NL, OUT_NL, VERT, MASS
      INTEGER :: I, J, L
      REAL*8  :: INDATA(I1x1,J1x1,IN_NL)
      REAL*8  :: OUTDATA(IIPAR,JJPAR,OUT_NL)
      REAL*8  :: F_ARRAY(LONF2F,LATG2F,IN_NL)
      REAL*8  :: P_ARRAY(LONF2S,LATG2S,IN_NL)
      REAL*8  :: TMP(IIPAR,JJPAR,IN_NL)

      IF (LMASTER) THEN

      ! Perform horizontal regridding
      DO L = 1, IN_NL

         ! 1x1 -> model grid
         CALL REGRID_1x1( INDATA(:,:,L), F_ARRAY(:,:,L), MASS )

      ENDDO

      ! Transform to model grid
      CALL XSHIFT(F_ARRAY,LONFD,LATGD,IN_NL)
      CALL YREV(F_ARRAY,LONFD,LATGD,IN_NL)
      CALL dyn_trans2model_grid(F_ARRAY,IN_NL)
      ENDIF

#ifdef MP
      ! MPI full grid to partial grid
      CALL mpgf2p(F_ARRAY,LONF2F,LATG2F,P_ARRAY,LONF2S,LATG2S,IN_NL)
#else
      P_ARRAY = F_ARRAY
#endif

      ! No vertical regridding
      IF ( IN_NL .EQ. OUT_NL ) THEN

         DO L = 1, IN_NL
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            OUTDATA(I,J,L) = P_ARRAY(I,J,L)
         ENDDO
         ENDDO
         ENDDO

      ! Perform vertical regridding
      ELSE

         ! Initialize
         TMP = 0d0

         ! Copy
         DO L = 1, IN_NL
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TMP(I,J,L) = P_ARRAY(I,J,L)
         ENDDO
         ENDDO
         ENDDO

         IF ( MASS == 1 ) THEN
            WRITE(6,*) 'Error: VERT /= 0 and MASS == 1 in trans_regrid'
            STOP 9999
         ENDIF

         ! 38 layers -> model grid
         IF ( IN_NL .EQ. 38 ) THEN
            CALL INTERP_VER_T( TMP, OUTDATA, VERT, MASS )
         ! 47 layers -> model grid
         ELSE IF ( IN_NL .EQ. 47 ) THEN
            CALL INTERP_VER_R( TMP, OUTDATA, VERT, MASS )
         ! 72 layers -> model grid
         ELSE IF ( IN_NL .EQ. 72 ) THEN
            CALL INTERP_VER_F( TMP, OUTDATA, VERT, MASS )
         ELSE
            WRITE(6,*) 'Error: Unknown IN_NL in trans_regrid'
            STOP 9999
         ENDIF

      ENDIF
#endif

      END SUBROUTINE TRANS_INPUT_1x1

      SUBROUTINE TRANS_INPUT_2x25(INDATA,OUTDATA,IN_NL,OUT_NL,VERT,MASS)

#ifndef RMP
      USE CMN_SIZE_MOD
      USE module_trans
      USE paramodel
#ifdef MP
      USE commpi
#endif

      IMPLICIT NONE

      INTEGER :: IN_NL, OUT_NL, VERT, MASS
      INTEGER :: I, J, L
      REAL*8  :: INDATA(I2x25,J2x25,IN_NL)
      REAL*8  :: OUTDATA(IIPAR,JJPAR,OUT_NL)
      REAL*8  :: F_ARRAY(LONF2F,LATG2F,IN_NL)
      REAL*8  :: P_ARRAY(LONF2S,LATG2S,IN_NL)
      REAL*8  :: TMP(IIPAR,JJPAR,IN_NL)

      IF (LMASTER) THEN

      ! Perform horizontal regridding
      DO L = 1, IN_NL

         ! 2x25 -> model grid
         CALL REGRID_2x25( INDATA(:,:,L), F_ARRAY(:,:,L), MASS )

      ENDDO

      ! Transform to model grid
      CALL XSHIFT(F_ARRAY,LONFD,LATGD,IN_NL)
      CALL YREV(F_ARRAY,LONFD,LATGD,IN_NL)
      CALL dyn_trans2model_grid(F_ARRAY,IN_NL)
      ENDIF

#ifdef MP
      ! MPI full grid to partial grid
      CALL mpgf2p(F_ARRAY,LONF2F,LATG2F,P_ARRAY,LONF2S,LATG2S,IN_NL)
#else
      P_ARRAY = F_ARRAY
#endif

      ! No vertical regridding
      IF ( IN_NL .EQ. OUT_NL ) THEN

         DO L = 1, IN_NL
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            OUTDATA(I,J,L) = P_ARRAY(I,J,L)
         ENDDO
         ENDDO
         ENDDO

      ! Perform vertical regridding
      ELSE

         ! Initialize
         TMP = 0d0

         ! Copy
         DO L = 1, IN_NL
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TMP(I,J,L) = P_ARRAY(I,J,L)
         ENDDO
         ENDDO
         ENDDO

         IF ( MASS == 1 ) THEN
            WRITE(6,*) 'Error: VERT /= 0 and MASS == 1 in trans_regrid'
            STOP 9999
         ENDIF

         ! 38 layers -> model grid
         IF ( IN_NL .EQ. 38 ) THEN
            CALL INTERP_VER_T( TMP, OUTDATA, VERT, MASS )
         ! 47 layers -> model grid
         ELSE IF ( IN_NL .EQ. 47 ) THEN
            CALL INTERP_VER_R( TMP, OUTDATA, VERT, MASS )
         ! 72 layers -> model grid
         ELSE IF ( IN_NL .EQ. 72 ) THEN
            CALL INTERP_VER_F( TMP, OUTDATA, VERT, MASS )
         ELSE
            WRITE(6,*) 'Error: Unknown IN_NL in trans_regrid'
            STOP 9999
         ENDIF

      ENDIF
#endif

      END SUBROUTINE TRANS_INPUT_2x25

      SUBROUTINE TRANS_INPUT_RST(INDATA,OUTDATA,IN_NL,OUT_NL,VERT,MASS)

#ifndef RMP
      USE CMN_SIZE_MOD
      USE module_trans
      USE paramodel
#ifdef MP
      USE commpi
#endif

      IMPLICIT NONE

      INTEGER :: IN_NL, OUT_NL, VERT, MASS
      INTEGER :: I, J, L
      REAL*8  :: INDATA(LONF2F,LATG2F,IN_NL)
      REAL*8  :: OUTDATA(IIPAR,JJPAR,OUT_NL)
      REAL*8  :: F_ARRAY(LONF2F,LATG2F,IN_NL)
      REAL*8  :: P_ARRAY(LONF2S,LATG2S,IN_NL)
      REAL*8  :: TMP(IIPAR,JJPAR,IN_NL)

      IF (LMASTER) THEN

      ! No horizontal regridding
      F_ARRAY = INDATA

      ! Transform to model grid
      CALL XSHIFT(F_ARRAY,LONFD,LATGD,IN_NL)
      CALL YREV(F_ARRAY,LONFD,LATGD,IN_NL)
      CALL dyn_trans2model_grid(F_ARRAY,IN_NL)
      ENDIF

#ifdef MP
      ! MPI full grid to partial grid
      CALL mpgf2p(F_ARRAY,LONF2F,LATG2F,P_ARRAY,LONF2S,LATG2S,IN_NL)
#else
      P_ARRAY = F_ARRAY
#endif

      ! No vertical regridding
      IF ( IN_NL .EQ. OUT_NL ) THEN

         DO L = 1, IN_NL
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            OUTDATA(I,J,L) = P_ARRAY(I,J,L)
         ENDDO
         ENDDO
         ENDDO

      ! Perform vertical regridding
      ELSE

         ! Initialize
         TMP = 0d0

         ! Copy
         DO L = 1, IN_NL
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TMP(I,J,L) = P_ARRAY(I,J,L)
         ENDDO
         ENDDO
         ENDDO

         IF ( MASS == 1 ) THEN
            WRITE(6,*) 'Error: VERT /= 0 and MASS == 1 in trans_regrid'
            STOP 9999
         ENDIF

         ! 38 layers -> model grid
         IF ( IN_NL .EQ. 38 ) THEN
            CALL INTERP_VER_T( TMP, OUTDATA, VERT, MASS )
         ! 47 layers -> model grid
         ELSE IF ( IN_NL .EQ. 47 ) THEN
            CALL INTERP_VER_R( TMP, OUTDATA, VERT, MASS )
         ! 72 layers -> model grid
         ELSE IF ( IN_NL .EQ. 72 ) THEN
            CALL INTERP_VER_F( TMP, OUTDATA, VERT, MASS )
         ELSE
            WRITE(6,*) 'Error: Unknown IN_NL in trans_regrid'
            STOP 9999
         ENDIF

      ENDIF
#endif

      END SUBROUTINE TRANS_INPUT_RST

#ifndef RMP
!=====================================================================
!
!     SUBROUTINE REGRID_1x1
!
!        This subroutine performs horizontal regridding.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE REGRID_1x1( INDATA, OUTDATA, MASS )

      USE comfgrid
      USE paramodel

      INTEGER, PARAMETER  :: IDM1=360,     JDM1=181,
     &                       IDM1P=IDM1+1, JDM1P=JDM1+1
      INTEGER             :: IDM2, JDM2, IDM2P, JDM2P
      INTEGER, INTENT(IN) :: MASS
      REAL, INTENT(IN)    :: INDATA(IDM1,JDM1)
      REAL, INTENT(OUT)   :: OUTDATA(LONFD,LATGD)
      REAL                :: RLON(IDM1P), RLAT(JDM1P)
      REAL, ALLOCATABLE   :: ALON(:), ALAT(:), BLAT(:)
      REAL, ALLOCATABLE   :: IR1(:), IR2(:), JR1(:), JR2(:)
      REAL, ALLOCATABLE   :: FSUM(:,:)
      REAL*8, PARAMETER   :: PI = 3.14159265358979323d0

      IDM2 = LONFD
      JDM2 = LATGD
      IDM2P = IDM2+1
      JDM2P = JDM2+1

      ALLOCATE(ALON(IDM2P), ALAT(JDM2P), BLAT(JDM2))
      ALLOCATE(IR1(IDM2), IR2(IDM2), JR1(JDM2), JR2(JDM2))
      ALLOCATE(FSUM(IDM2,JDM2))

      !==============================================================
      ! Set up input data lon, lat of each grid box edges in radians
      !==============================================================
      DRLON = 2d0 * PI / DBLE(IDM1)
      DRLAT =       PI / DBLE(JDM1-1)
      DO I=1, IDM1P
         RLON(I) = DRLON * ( I - 1.5 )
      ENDDO

      RLAT(1) = PI / 2d0
      DO J=2, JDM1
         RLAT(J) = PI / 2d0 - DRLAT * ( J - 1.5 )
      ENDDO
      RLAT(JDM1P) = - ( PI / 2d0 )

      !==============================================================
      ! Set up model data lon, lat of each grid box edges in radians
      !==============================================================
      DALON = 2d0 * PI / DBLE(IDM2)
      DO I=1, IDM2P
         ALON(I) = DALON * ( I - 1.5 )
      ENDDO

      DO J=1, JDM2/2
         BLAT(J) = PI / 2d0 - colrad(J)
      ENDDO
      DO J=JDM2/2+1, JDM2
         BLAT(J) = - BLAT(JDM2+1-J)
      ENDDO

      ALAT(1) = PI / 2d0
      DO J=2, JDM2
         ALAT(J) = 0.5d0 * ( BLAT(J-1) + BLAT(J) )
      ENDDO
      ALAT(JDM2P) = - ( PI / 2d0 )

      !========================================================
      ! Find starting and ending indices of the input data for
      ! each of the model grid boxes
      !========================================================
      IR = 1
      DO 100 IA=1, IDM2
         IR1(IA) = 0
         IR2(IA) = 0
         IF (ALON(IA) .LT. 0.0) THEN
            XLON = 2d0 * PI + ALON(IA)
            IR0 = IDM1
80          IF (XLON .GE. RLON(IR0) .AND.
     &          XLON .LE. RLON(IR0+1)) THEN
               IR1(IA) = IR0
               GO TO 90
            ELSE
               IR0 = IR0 - 1
               IF ( IR0 .EQ. 1 ) THEN
                  IR1(IA) = IR0
                  GO TO 90
               ENDIF
               GO TO 80
            ENDIF
         ENDIF

90       IF (IR1(IA) .EQ. 0 .AND. ALON(IA) .GE. RLON(IR)
     &       .AND. ALON(IA) .LE. RLON(IR+1)) THEN
            IR1(IA) = IR
         ENDIF
         IF (IR2(IA) .EQ. 0 .AND. ALON(IA+1) .GE. RLON(IR)
     &       .AND. ALON(IA+1) .LE. RLON(IR+1)) THEN
            IR2(IA) = IR
            GO TO 100
         ELSE
            IF (IR .GE. IDM1) THEN
               IF (IR1(IA) .EQ. 0) IR1(IA) = 1
               IF (IR2(IA) .EQ. 0) IR2(IA) = 1
               GO TO 100
            ELSE
               IR = IR + 1
               GO TO 90
            ENDIF
         ENDIF
100   CONTINUE

      JR = 1
      DO 160 JA=1, JDM2
         JR1(JA) = 0
         JR2(JA) = 0
150      IF (JR1(JA) .EQ. 0 .AND. ALAT(JA) .LE. RLAT(JR)
     &       .AND. ALAT(JA) .GE. RLAT(JR+1)) THEN
            JR1(JA) = JR
         ENDIF
         IF (JR2(JA) .EQ. 0 .AND. ALAT(JA+1) .LE. RLAT(JR)
     &       .AND. ALAT(JA+1) .GE. RLAT(JR+1)) THEN
            JR2(JA) = JR
            GO TO 160
         ELSE
            IF (JR .GE. JDM1) THEN
               JR2(JA) = JR
               GO TO 180
            ELSE
               JR = JR + 1
               GO TO 150
            ENDIF
         ENDIF
160   CONTINUE

180   CONTINUE

      !================
      ! Initialization
      !================
      DO J=1, JDM2
         DO I=1, IDM2
            OUTDATA(I,J) = 0d0
         ENDDO
      ENDDO

      !========================
      ! Map data to model grid
      !========================
      DO JA=1, JDM2

      !===================================================
      ! Fraction of coverage in j direction (cos weightd)
      !===================================================
         JJR1 = JR1(JA)
         JJR2 = JR2(JA)

      IF (JJR1 .EQ. JJR2) THEN
         FJR1 = (COS(ALAT(JA)) - COS(ALAT(JA+1)))
     &        / (COS(RLAT(JJR1)) - COS(RLAT(JJR1+1)))
         IF (COS(RLAT(JJR1)) - COS(RLAT(JJR1+1)) .EQ. 0) THEN
            FJR1 = ABS(COS(ALAT(JA)) - COS(ALAT(JA+1)))
     &           / ( 2d0 * ( 1d0 - COS(RLAT(JJR1)) ) )
         ENDIF
         FJR2 = 0d0
      ELSE
         FJR1 = (COS(ALAT(JA)) - COS(RLAT(JJR1+1)))
     &        / (COS(RLAT(JJR1)) - COS(RLAT(JJR1+1)))
         IF (COS(RLAT(JJR1)) - COS(RLAT(JJR1+1)) .EQ. 0) THEN
            FJR1 = ABS(COS(ALAT(JA)) - COS(RLAT(JJR1+1)))
     &           / ( 2d0 * ( 1d0 - COS(RLAT(JJR1)) ) )
         ENDIF
         FJR2 = (COS(RLAT(JJR2)) - COS(ALAT(JA+1)))
     &        / (COS(RLAT(JJR2)) - COS(RLAT(JJR2+1)))
         IF (COS(RLAT(JJR2)) - COS(RLAT(JJR2+1)) .EQ. 0) THEN
            FJR2 = ABS(COS(RLAT(JJR2)) - COS(ALAT(JA+1)))
     &           / ( 2d0 * ( 1d0 - COS(RLAT(JJR2)) ) )
         ENDIF
      ENDIF

      DO IA=1, IDM2

      !=====================================
      ! Fraction of coverage in i direction
      !=====================================
         IIR1 = IR1(IA)
         IIR2 = IR2(IA)

         IF (IIR1 .EQ. IIR2) THEN
            FIR1 = (ALON(IA+1) - ALON(IA))
     &           / (RLON(IIR1+1) - RLON(IIR1))
            FIR2 = 0d0
         ELSE
            FIR1 = (RLON(IIR1+1) - ALON(IA))
     &           / (RLON(IIR1+1) - RLON(IIR1))
            IF (RLON(IIR1+1) - ALON(IA) .GT. 2d0 * PI ) THEN
               FIR1 = (RLON(IIR1+1) - (2d0*PI+ALON(IA)))
     &              / (RLON(IIR1+1) - RLON(IIR1))
            ENDIF
            FIR2 = (ALON(IA+1) - RLON(IIR2))
     &           / (RLON(IIR2+1) - RLON(IIR2))
            IF (ALON(IA+1) - RLON(IIR2) .GT. 2d0 * PI ) THEN
               FIR2 = (ALON(IA+1) - (2d0*PI+RLON(IIR2)))
     &              / (RLON(IIR2+1) - RLON(IIR2))
            ENDIF
         ENDIF

         !=====================
         ! First latitude band
         !=====================
         OUTDATA(IA,JA) = INDATA(IIR1,JJR1) * FIR1 * FJR1
         FSUM(IA,JA) = FIR1 * FJR1
         IF ( (IIR2 - IIR1) .LT. 0 ) THEN
            DO IR=IIR1+1, IDM1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR1) * FJR1
               FSUM(IA,JA) = FSUM(IA,JA) + FJR1
            ENDDO
            DO IR=1, IIR2-1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR1) * FJR1
               FSUM(IA,JA) = FSUM(IA,JA) + FJR1
            ENDDO
            GO TO 200
         ENDIF
         DO IR=IIR1+1, IIR2-1
            OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR1) * FJR1
            FSUM(IA,JA) = FSUM(IA,JA) + FJR1
         ENDDO

200      CONTINUE

         OUTDATA(IA,JA) = OUTDATA(IA,JA)
     &                    + INDATA(IIR2,JJR1) * FIR2 * FJR1
         FSUM(IA,JA) = FSUM(IA,JA) + FIR2 * FJR1

         !=======================
         ! Consecutive latitudes
         !=======================
         DO JR=JJR1+1, JJR2-1

            OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IIR1,JR) * FIR1
            FSUM(IA,JA) = FSUM(IA,JA) + FIR1
            IF ( (IIR2 - IIR1) .LT. 0 ) THEN
               DO IR=IIR1+1, IDM1
                  OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JR)
                  FSUM(IA,JA) = FSUM(IA,JA) + 1.0
               ENDDO
               DO IR=1, IIR2-1
                  OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JR)
                  FSUM(IA,JA) = FSUM(IA,JA) + 1.0
               ENDDO
               GO TO 210
            ENDIF
            DO IR=IIR1+1, IIR2-1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JR)
               FSUM(IA,JA) = FSUM(IA,JA) + 1.0
            ENDDO

210         CONTINUE

            OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IIR2,JR) * FIR2
            FSUM(IA,JA) = FSUM(IA,JA) + FIR2

         ENDDO

         !===============
         ! Last latitude
         !===============
         OUTDATA(IA,JA) = OUTDATA(IA,JA) +
     &                    INDATA(IIR1,JJR2) * FIR1 * FJR2
         FSUM(IA,JA) = FSUM(IA,JA) + FIR1 * FJR2
         IF ( (IIR2 - IIR1) .LT. 0 ) THEN
            DO IR=IIR1+1, IDM1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR2) * FJR2
               FSUM(IA,JA) = FSUM(IA,JA) + FJR2
            ENDDO
            DO IR=1, IIR2-1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR2) * FJR2
               FSUM(IA,JA) = FSUM(IA,JA) + FJR2
            ENDDO
            GO TO 220
         ENDIF
         DO IR=IIR1+1, IIR2-1
            OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR2) * FJR2
            FSUM(IA,JA) = FSUM(IA,JA) + FJR2
         ENDDO

220      CONTINUE

         OUTDATA(IA,JA) = OUTDATA(IA,JA) +
     &                    INDATA(IIR2,JJR2) * FIR2 * FJR2
         FSUM(IA,JA) = FSUM(IA,JA) + FIR2 * FJR2

      ENDDO
      ENDDO

      !===========================================
      ! Average if it is concentration (MASS = 0)
      !===========================================
      IF (MASS .EQ. 0) OUTDATA(:,:) = OUTDATA(:,:) / FSUM(:,:)

      DEALLOCATE(ALON, ALAT, BLAT)
      DEALLOCATE(IR1, IR2, JR1, JR2)
      DEALLOCATE(FSUM)

      END SUBROUTINE REGRID_1x1

!=====================================================================
!
!     SUBROUTINE REGRID_2x25
!
!        This subroutine performs horizontal regridding.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE REGRID_2x25( INDATA, OUTDATA, MASS )

      USE comfgrid
      USE paramodel

      INTEGER, PARAMETER  :: IDM1=144,     JDM1=91,
     &                       IDM1P=IDM1+1, JDM1P=JDM1+1
      INTEGER             :: IDM2, JDM2, IDM2P, JDM2P
      INTEGER, INTENT(IN) :: MASS
      REAL, INTENT(IN)    :: INDATA(IDM1,JDM1)
      REAL, INTENT(OUT)   :: OUTDATA(LONFD,LATGD)
      REAL                :: RLON(IDM1P), RLAT(JDM1P)
      REAL, ALLOCATABLE   :: ALON(:), ALAT(:), BLAT(:)
      REAL, ALLOCATABLE   :: IR1(:), IR2(:), JR1(:), JR2(:)
      REAL, ALLOCATABLE   :: FSUM(:,:)
      REAL*8, PARAMETER   :: PI = 3.14159265358979323d0

      IDM2 = LONFD
      JDM2 = LATGD
      IDM2P = IDM2+1
      JDM2P = JDM2+1

      ALLOCATE(ALON(IDM2P), ALAT(JDM2P), BLAT(JDM2))
      ALLOCATE(IR1(IDM2), IR2(IDM2), JR1(JDM2), JR2(JDM2))
      ALLOCATE(FSUM(IDM2,JDM2))

      !==============================================================
      ! Set up input data lon, lat of each grid box edges in radians
      !==============================================================
      DRLON = 2d0 * PI / DBLE(IDM1)
      DRLAT =       PI / DBLE(JDM1-1)
      DO I=1, IDM1P
         RLON(I) = DRLON * ( I - 1.5 )
      ENDDO

      RLAT(1) = PI / 2d0
      DO J=2, JDM1
         RLAT(J) = PI / 2d0 - DRLAT * ( J - 1.5 )
      ENDDO
      RLAT(JDM1P) = - ( PI / 2d0 )

      !==============================================================
      ! Set up model data lon, lat of each grid box edges in radians
      !==============================================================
      DALON = 2d0 * PI / DBLE(IDM2)
      DO I=1, IDM2P
         ALON(I) = DALON * ( I - 1.5 )
      ENDDO

      DO J=1, JDM2/2
         BLAT(J) = PI / 2d0 - colrad(J)
      ENDDO
      DO J=JDM2/2+1, JDM2
         BLAT(J) = - BLAT(JDM2+1-J)
      ENDDO

      ALAT(1) = PI / 2d0
      DO J=2, JDM2
         ALAT(J) = 0.5d0 * ( BLAT(J-1) + BLAT(J) )
      ENDDO
      ALAT(JDM2P) = - ( PI / 2d0 )

      !========================================================
      ! Find starting and ending indices of the input data for
      ! each of the model grid boxes
      !========================================================
      IR = 1
      DO 100 IA=1, IDM2
         IR1(IA) = 0
         IR2(IA) = 0
         IF (ALON(IA) .LT. 0.0) THEN
            XLON = 2d0 * PI + ALON(IA)
            IR0 = IDM1
80          IF (XLON .GE. RLON(IR0) .AND.
     &          XLON .LE. RLON(IR0+1)) THEN
               IR1(IA) = IR0
               GO TO 90
            ELSE
               IR0 = IR0 - 1
               IF ( IR0 .EQ. 1 ) THEN
                  IR1(IA) = IR0
                  GO TO 90
               ENDIF
               GO TO 80
            ENDIF
         ENDIF

90       IF (IR1(IA) .EQ. 0 .AND. ALON(IA) .GE. RLON(IR)
     &       .AND. ALON(IA) .LE. RLON(IR+1)) THEN
            IR1(IA) = IR
         ENDIF
         IF (IR2(IA) .EQ. 0 .AND. ALON(IA+1) .GE. RLON(IR)
     &       .AND. ALON(IA+1) .LE. RLON(IR+1)) THEN
            IR2(IA) = IR
            GO TO 100
         ELSE
            IF (IR .GE. IDM1) THEN
               IF (IR1(IA) .EQ. 0) IR1(IA) = 1
               IF (IR2(IA) .EQ. 0) IR2(IA) = 1
               GO TO 100
            ELSE
               IR = IR + 1
               GO TO 90
            ENDIF
         ENDIF
100   CONTINUE

      JR = 1
      DO 160 JA=1, JDM2
         JR1(JA) = 0
         JR2(JA) = 0
150      IF (JR1(JA) .EQ. 0 .AND. ALAT(JA) .LE. RLAT(JR)
     &       .AND. ALAT(JA) .GE. RLAT(JR+1)) THEN
            JR1(JA) = JR
         ENDIF
         IF (JR2(JA) .EQ. 0 .AND. ALAT(JA+1) .LE. RLAT(JR)
     &       .AND. ALAT(JA+1) .GE. RLAT(JR+1)) THEN
            JR2(JA) = JR
            GO TO 160
         ELSE
            IF (JR .GE. JDM1) THEN
               JR2(JA) = JR
               GO TO 180
            ELSE
               JR = JR + 1
               GO TO 150
            ENDIF
         ENDIF
160   CONTINUE

180   CONTINUE

      !================
      ! Initialization
      !================
      DO J=1, JDM2
         DO I=1, IDM2
            OUTDATA(I,J) = 0d0
         ENDDO
      ENDDO

      !========================
      ! Map data to model grid
      !========================
      DO JA=1, JDM2

      !===================================================
      ! Fraction of coverage in j direction (cos weightd)
      !===================================================
         JJR1 = JR1(JA)
         JJR2 = JR2(JA)

      IF (JJR1 .EQ. JJR2) THEN
         FJR1 = (COS(ALAT(JA)) - COS(ALAT(JA+1)))
     &        / (COS(RLAT(JJR1)) - COS(RLAT(JJR1+1)))
         IF (COS(RLAT(JJR1)) - COS(RLAT(JJR1+1)) .EQ. 0) THEN
            FJR1 = ABS(COS(ALAT(JA)) - COS(ALAT(JA+1)))
     &           / ( 2d0 * ( 1d0 - COS(RLAT(JJR1)) ) )
         ENDIF
         FJR2 = 0d0
      ELSE
         FJR1 = (COS(ALAT(JA)) - COS(RLAT(JJR1+1)))
     &        / (COS(RLAT(JJR1)) - COS(RLAT(JJR1+1)))
         IF (COS(RLAT(JJR1)) - COS(RLAT(JJR1+1)) .EQ. 0) THEN
            FJR1 = ABS(COS(ALAT(JA)) - COS(RLAT(JJR1+1)))
     &           / ( 2d0 * ( 1d0 - COS(RLAT(JJR1)) ) )
         ENDIF
         FJR2 = (COS(RLAT(JJR2)) - COS(ALAT(JA+1)))
     &        / (COS(RLAT(JJR2)) - COS(RLAT(JJR2+1)))
         IF (COS(RLAT(JJR2)) - COS(RLAT(JJR2+1)) .EQ. 0) THEN
            FJR2 = ABS(COS(RLAT(JJR2)) - COS(ALAT(JA+1)))
     &           / ( 2d0 * ( 1d0 - COS(RLAT(JJR2)) ) )
         ENDIF
      ENDIF

      DO IA=1, IDM2

      !=====================================
      ! Fraction of coverage in i direction
      !=====================================
         IIR1 = IR1(IA)
         IIR2 = IR2(IA)

         IF (IIR1 .EQ. IIR2) THEN
            FIR1 = (ALON(IA+1) - ALON(IA))
     &           / (RLON(IIR1+1) - RLON(IIR1))
            FIR2 = 0d0
         ELSE
            FIR1 = (RLON(IIR1+1) - ALON(IA))
     &           / (RLON(IIR1+1) - RLON(IIR1))
            IF (RLON(IIR1+1) - ALON(IA) .GT. 2d0 * PI ) THEN
               FIR1 = (RLON(IIR1+1) - (2d0*PI+ALON(IA)))
     &              / (RLON(IIR1+1) - RLON(IIR1))
            ENDIF
            FIR2 = (ALON(IA+1) - RLON(IIR2))
     &           / (RLON(IIR2+1) - RLON(IIR2))
            IF (ALON(IA+1) - RLON(IIR2) .GT. 2d0 * PI ) THEN
               FIR2 = (ALON(IA+1) - (2d0*PI+RLON(IIR2)))
     &              / (RLON(IIR2+1) - RLON(IIR2))
            ENDIF
         ENDIF

         !=====================
         ! First latitude band
         !=====================
         OUTDATA(IA,JA) = INDATA(IIR1,JJR1) * FIR1 * FJR1
         FSUM(IA,JA) = FIR1 * FJR1
         IF ( (IIR2 - IIR1) .LT. 0 ) THEN
            DO IR=IIR1+1, IDM1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR1) * FJR1
               FSUM(IA,JA) = FSUM(IA,JA) + FJR1
            ENDDO
            DO IR=1, IIR2-1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR1) * FJR1
               FSUM(IA,JA) = FSUM(IA,JA) + FJR1
            ENDDO
            GO TO 200
         ENDIF
         DO IR=IIR1+1, IIR2-1
            OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR1) * FJR1
            FSUM(IA,JA) = FSUM(IA,JA) + FJR1
         ENDDO

200      CONTINUE

         OUTDATA(IA,JA) = OUTDATA(IA,JA)
     &                    + INDATA(IIR2,JJR1) * FIR2 * FJR1
         FSUM(IA,JA) = FSUM(IA,JA) + FIR2 * FJR1

         !=======================
         ! Consecutive latitudes
         !=======================
         DO JR=JJR1+1, JJR2-1

            OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IIR1,JR) * FIR1
            FSUM(IA,JA) = FSUM(IA,JA) + FIR1
            IF ( (IIR2 - IIR1) .LT. 0 ) THEN
               DO IR=IIR1+1, IDM1
                  OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JR)
                  FSUM(IA,JA) = FSUM(IA,JA) + 1.0
               ENDDO
               DO IR=1, IIR2-1
                  OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JR)
                  FSUM(IA,JA) = FSUM(IA,JA) + 1.0
               ENDDO
               GO TO 210
            ENDIF
            DO IR=IIR1+1, IIR2-1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JR)
               FSUM(IA,JA) = FSUM(IA,JA) + 1.0
            ENDDO

210         CONTINUE

            OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IIR2,JR) * FIR2
            FSUM(IA,JA) = FSUM(IA,JA) + FIR2

         ENDDO

         !===============
         ! Last latitude
         !===============
         OUTDATA(IA,JA) = OUTDATA(IA,JA) +
     &                    INDATA(IIR1,JJR2) * FIR1 * FJR2
         FSUM(IA,JA) = FSUM(IA,JA) + FIR1 * FJR2
         IF ( (IIR2 - IIR1) .LT. 0 ) THEN
            DO IR=IIR1+1, IDM1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR2) * FJR2
               FSUM(IA,JA) = FSUM(IA,JA) + FJR2
            ENDDO
            DO IR=1, IIR2-1
               OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR2) * FJR2
               FSUM(IA,JA) = FSUM(IA,JA) + FJR2
            ENDDO
            GO TO 220
         ENDIF
         DO IR=IIR1+1, IIR2-1
            OUTDATA(IA,JA) = OUTDATA(IA,JA) + INDATA(IR,JJR2) * FJR2
            FSUM(IA,JA) = FSUM(IA,JA) + FJR2
         ENDDO

220      CONTINUE

         OUTDATA(IA,JA) = OUTDATA(IA,JA) +
     &                    INDATA(IIR2,JJR2) * FIR2 * FJR2
         FSUM(IA,JA) = FSUM(IA,JA) + FIR2 * FJR2

      ENDDO
      ENDDO

      !===========================================
      ! Average if it is concentration (MASS = 0)
      !===========================================
      IF (MASS .EQ. 0) OUTDATA(:,:) = OUTDATA(:,:) / FSUM(:,:)

      DEALLOCATE(ALON, ALAT, BLAT)
      DEALLOCATE(IR1, IR2, JR1, JR2)
      DEALLOCATE(FSUM)

      END SUBROUTINE REGRID_2x25

!=====================================================================
!
!     SUBROUTINE INTERP_VER_T
!
!        This subroutine is for interpolating GEOS-Chem reduced
!        tropospheric vertical pressure grid to GRIMs vertical
!        pressure grid.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE INTERP_VER_T( INDATA, OUTDATA, VERT, MASS )

      USE CMN_SIZE_MOD
      USE comfver

      IMPLICIT NONE

      INTEGER :: VERT, MASS
      INTEGER :: I, J, L, ID, ID1, ID2
      REAL*8  :: INDATA(IIPAR,JJPAR,38)
      REAL*8  :: OUTDATA(IIPAR,JJPAR,LLTROP)
      REAL*8  :: AP(39), BP(39)
      REAL*8  :: PEDGE1(39), PEDGE2(LLTROP+1)
      REAL*8  :: FSUM(LLTROP)
      REAL*8  :: FR, FR1, FR2

      AP = (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,
     &        1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,
     &        4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,
     &        7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,
     &        1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,
     &        1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,
     &        2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,
     &        2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02,
     &        1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01,
     &        7.851231d+01, 5.638791d+01, 4.017541d+01 /)

      BP = (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01,
     &        9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01,
     &        8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01,
     &        7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01,
     &        6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01,
     &        4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01,
     &        2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01,
     &        6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00 /)

      !================
      ! Initialization
      !================
      OUTDATA = 0d0
      FSUM = 0d0

      !===============================================
      ! Set up GEOS-5 reduced vertical pressure level
      !===============================================
      DO L=1, 39
         PEDGE1(L) = AP(L) + BP(L) * 1013.25
      ENDDO

      !======================================
      ! Set up GRIMs vertical pressure level
      !======================================
      DO L=1, LLTROP+1
         PEDGE2(L) = ak5(LLPAR+2-L) * 1d1 + bk5(LLPAR+2-L) * 1013.25
         IF (PEDGE2(L) .LE. 0) PEDGE2(L) = 0.01d0
      ENDDO

      DO L=1, LLTROP

         ID1=0
         ID2=0

         !=====================================================
         ! Finding GEOS-5 reduced pressure level indices
         ! which correspond to the tip of GRIMs pressure level
         !=====================================================
         DO ID=1, 38
         IF (PEDGE2(L) .LE. PEDGE1(ID) .AND.
     &       PEDGE2(L) .GT. PEDGE1(ID+1)) THEN
            ID1=ID
         ENDIF
         IF (PEDGE2(L+1) .LT. PEDGE1(ID) .AND.
     &       PEDGE2(L+1) .GE. PEDGE1(ID+1)) THEN
            ID2=ID
         ENDIF
         ENDDO

         !=======================================================
         ! Calculate coverage at the tip of GRIMs pressure level
         !=======================================================
         IF (VERT .EQ. 1) THEN

            ! z weighted
            IF (ID1 .EQ. ID2) THEN
               FR1 = LOG(PEDGE2(L)) - LOG(PEDGE2(L+1))
               FR2 = 0d0
            ELSE
               FR1 = LOG(PEDGE2(L)) - LOG(PEDGE1(ID1+1))
               FR2 = LOG(PEDGE1(ID2)) - LOG(PEDGE2(L+1))
            ENDIF

         ELSEIF (VERT .EQ. 2) THEN

            ! p weighted
            IF (ID1 .EQ. ID2) THEN
               FR1 = PEDGE2(L) - PEDGE2(L+1)
               FR2 = 0d0
            ELSE
               FR1 = PEDGE2(L) - PEDGE1(ID1+1)
               FR2 = PEDGE1(ID2) - PEDGE2(L+1)
            ENDIF

         ENDIF

         !============================
         ! Interpolate pressure level
         !============================
         OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID1) * FR1
         FSUM(L) = FSUM(L) + FR1

         DO ID=ID1+1, ID2-1
            FR = PEDGE1(ID) - PEDGE1(ID+1)
            OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID) * FR
            FSUM(L) = FSUM(L) + FR
         ENDDO

         OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID2) * FR2
         FSUM(L) = FSUM(L) + FR2

      ENDDO

      !================================
      ! Average if it is concentration
      !================================
      IF (MASS .EQ. 0) THEN
         DO L=1, LLTROP
         DO J=1, JJPAR
         DO I=1, IIPAR
            OUTDATA(I,J,L) = OUTDATA(I,J,L) / FSUM(L)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE INTERP_VER_T

!=====================================================================
!
!     SUBROUTINE INTERP_VER_R
!
!        This subroutine is for interpolating GEOS-Chem reduced
!        vertical pressure grid to GRIMs vertical pressure grid.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE INTERP_VER_R( INDATA, OUTDATA, VERT, MASS )

      USE CMN_SIZE_MOD
      USE comfver

      IMPLICIT NONE

      INTEGER :: VERT, MASS
      INTEGER :: I, J, L, ID, ID1, ID2
      REAL*8  :: INDATA(IIPAR,JJPAR,47)
      REAL*8  :: OUTDATA(IIPAR,JJPAR,LLPAR)
      REAL*8  :: AP(48), BP(48)
      REAL*8  :: PEDGE1(48), PEDGE2(LLPAR+1)
      REAL*8  :: FSUM(LLPAR)
      REAL*8  :: FR, FR1, FR2

      AP = (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,
     &        1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,
     &        4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,
     &        7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,
     &        1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,
     &        1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,
     &        2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,
     &        2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02,
     &        1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01,
     &        7.851231d+01, 5.638791d+01, 4.017541d+01, 2.836781d+01,
     &        1.979160d+01, 9.292942d+00, 4.076571d+00, 1.650790d+00,
     &        6.167791d-01, 2.113490d-01, 6.600001d-02, 1.000000d-02 /)

      BP = (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01,
     &        9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01,
     &        8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01,
     &        7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01,
     &        6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01,
     &        4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01,
     &        2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01,
     &        6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00 /)

      !================
      ! Initialization
      !================
      OUTDATA = 0d0
      FSUM = 0d0

      !===============================================
      ! Set up GEOS-5 reduced vertical pressure level
      !===============================================
      DO L=1, 48
         PEDGE1(L) = AP(L) + BP(L) * 1013.25
      ENDDO

      !======================================
      ! Set up GRIMs vertical pressure level
      !======================================
      DO L=1, LLPAR+1
         PEDGE2(L) = ak5(LLPAR+2-L) * 1d1 + bk5(LLPAR+2-L) * 1013.25
         IF (PEDGE2(L) .LE. 0) PEDGE2(L) = 0.01d0
      ENDDO

      DO L=1, LLPAR

         ID1=0
         ID2=0

         !=====================================================
         ! Finding GEOS-5 reduced pressure level indices
         ! which correspond to the tip of GRIMs pressure level
         !=====================================================
         DO ID=1, 47
         IF (PEDGE2(L) .LE. PEDGE1(ID) .AND.
     &       PEDGE2(L) .GT. PEDGE1(ID+1)) THEN
            ID1=ID
         ENDIF
         IF (PEDGE2(L+1) .LT. PEDGE1(ID) .AND.
     &       PEDGE2(L+1) .GE. PEDGE1(ID+1)) THEN
            ID2=ID
         ENDIF
         ENDDO

         !=======================================================
         ! Calculate coverage at the tip of GRIMs pressure level
         !=======================================================
         IF (VERT .EQ. 1) THEN

            ! z weighted
            IF (ID1 .EQ. ID2) THEN
               FR1 = LOG(PEDGE2(L)) - LOG(PEDGE2(L+1))
               FR2 = 0d0
            ELSE
               FR1 = LOG(PEDGE2(L)) - LOG(PEDGE1(ID1+1))
               FR2 = LOG(PEDGE1(ID2)) - LOG(PEDGE2(L+1))
            ENDIF

         ELSEIF (VERT .EQ. 2) THEN

            ! p weighted
            IF (ID1 .EQ. ID2) THEN
               FR1 = PEDGE2(L) - PEDGE2(L+1)
               FR2 = 0d0
            ELSE
               FR1 = PEDGE2(L) - PEDGE1(ID1+1)
               FR2 = PEDGE1(ID2) - PEDGE2(L+1)
            ENDIF

         ENDIF

         !============================
         ! Interpolate pressure level
         !============================
         OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID1) * FR1
         FSUM(L) = FSUM(L) + FR1

         DO ID=ID1+1, ID2-1
            FR = PEDGE1(ID) - PEDGE1(ID+1)
            OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID) * FR
            FSUM(L) = FSUM(L) + FR
         ENDDO

         OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID2) * FR2
         FSUM(L) = FSUM(L) + FR2

      ENDDO

      !================================
      ! Average if it is concentration
      !================================
      IF (MASS .EQ. 0) THEN
         DO L=1, LLPAR
         DO J=1, JJPAR
         DO I=1, IIPAR
            OUTDATA(I,J,L) = OUTDATA(I,J,L) / FSUM(L)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE INTERP_VER_R

!=====================================================================
!
!     SUBROUTINE INTERP_VER_F
!
!        This subroutine is for interpolating GEOS-Chem full
!        vertical pressure grid to GRIMs vertical pressure grid.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE INTERP_VER_F( INDATA, OUTDATA, VERT, MASS )

      USE CMN_SIZE_MOD
      USE comfver

      IMPLICIT NONE

      INTEGER :: VERT, MASS
      INTEGER :: I, J, L, ID, ID1, ID2
      REAL*8  :: INDATA(IIPAR,JJPAR,72)
      REAL*8  :: OUTDATA(IIPAR,JJPAR,LLPAR)
      REAL*8  :: AP(73), BP(73)
      REAL*8  :: PEDGE1(73), PEDGE2(LLPAR+1)
      REAL*8  :: FSUM(LLPAR)
      REAL*8  :: FR, FR1, FR2

      AP = (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,
     &        1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,
     &        4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,
     &        7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,
     &        1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,
     &        1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,
     &        2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,
     &        2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02,
     &        1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01,
     &        7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01,
     &        4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01,
     &        1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01,
     &        9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00,
     &        4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00,
     &        1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01,
     &        6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01,
     &        2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02,
     &        6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02,
     &        1.000000d-02 /)

      BP = (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01,
     &        9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01,
     &        8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01,
     &        7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01,
     &        6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01,
     &        4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01,
     &        2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01,
     &        6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00 /)

      !================
      ! Initialization
      !================
      OUTDATA = 0d0
      FSUM = 0d0

      !============================================
      ! Set up GEOS-5 full vertical pressure level
      !============================================
      DO L=1, 73
         PEDGE1(L) = AP(L) + BP(L) * 1013.25
      ENDDO

      !======================================
      ! Set up GRIMs vertical pressure level
      !======================================
      DO L=1, LLPAR+1
         PEDGE2(L) = ak5(LLPAR+2-L) * 1d1 + bk5(LLPAR+2-L) * 1013.25
         IF (PEDGE2(L) .LE. 0) PEDGE2(L) = 0.01d0
      ENDDO

      DO L=1, LLPAR

         ID1=0
         ID2=0

         !=====================================================
         ! Finding GEOS-5 full pressure level indices
         ! which correspond to the tip of GRIMs pressure level
         !=====================================================
         DO ID=1, 72
         IF (PEDGE2(L) .LE. PEDGE1(ID) .AND.
     &       PEDGE2(L) .GT. PEDGE1(ID+1)) THEN
            ID1=ID
         ENDIF
         IF (PEDGE2(L+1) .LT. PEDGE1(ID) .AND.
     &       PEDGE2(L+1) .GE. PEDGE1(ID+1)) THEN
            ID2=ID
         ENDIF
         ENDDO

         !=======================================================
         ! Calculate coverage at the tip of GRIMs pressure level
         !=======================================================
         IF (VERT .EQ. 1) THEN

            ! z weighted
            IF (ID1 .EQ. ID2) THEN
               FR1 = LOG(PEDGE2(L)) - LOG(PEDGE2(L+1))
               FR2 = 0d0
            ELSE
               FR1 = LOG(PEDGE2(L)) - LOG(PEDGE1(ID1+1))
               FR2 = LOG(PEDGE1(ID2)) - LOG(PEDGE2(L+1))
            ENDIF

         ELSEIF (VERT .EQ. 2) THEN

            ! p weighted
            IF (ID1 .EQ. ID2) THEN
               FR1 = PEDGE2(L) - PEDGE2(L+1)
               FR2 = 0d0
            ELSE
               FR1 = PEDGE2(L) - PEDGE1(ID1+1)
               FR2 = PEDGE1(ID2) - PEDGE2(L+1)
            ENDIF

         ENDIF

         !============================
         ! Interpolate pressure level
         !============================
         OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID1) * FR1
         FSUM(L) = FSUM(L) + FR1

         DO ID=ID1+1, ID2-1
            FR = PEDGE1(ID) - PEDGE1(ID+1)
            OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID) * FR
            FSUM(L) = FSUM(L) + FR
         ENDDO

         OUTDATA(:,:,L) = OUTDATA(:,:,L) + INDATA(:,:,ID2) * FR2
         FSUM(L) = FSUM(L) + FR2

      ENDDO

      !================================
      ! Average if it is concentration
      !================================
      IF (MASS .EQ. 0) THEN
         DO L=1, LLPAR
         DO J=1, JJPAR
         DO I=1, IIPAR
            OUTDATA(I,J,L) = OUTDATA(I,J,L) / FSUM(L)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE INTERP_VER_F
#endif
