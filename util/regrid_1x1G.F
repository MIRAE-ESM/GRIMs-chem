      !=========================================================
      ! This code is originally from the GRIMs 'cvnaer.f90'.
      ! It is modified for use of regridding emission data.
      !                                 - Seungun Lee -
      !=========================================================
      SUBROUTINE REGRID_1x1G( INDATA, OUTDATA, FLAG )

      USE CMN_SIZE_MOD

      INTEGER, PARAMETER  :: IDM1=360,     JDM1=180,
     &                       IDM1P=IDM1+1, JDM1P=JDM1+1
      INTEGER             :: IDM2, JDM2, IDM2P, JDM2P
      REAL, INTENT(IN)    :: INDATA(IDM1,JDM1)
      REAL, INTENT(OUT)   :: OUTDATA(IGLOB,JGLOB)
      REAL, ALLOCATABLE   :: FSUM(:,:)
      REAL                :: RLON(IDM1P), RLAT(JDM1P)
      REAL, ALLOCATABLE   :: ALON(:), ALAT(:), BLAT(:)
      REAL, ALLOCATABLE   :: IR1(:), IR2(:), JR1(:), JR2(:)
      INTEGER, INTENT(IN) :: FLAG
      REAL*8, PARAMETER   :: PI = 3.14159265358979323d0

      IDM2 = IGLOB
      JDM2 = JGLOB
      IDM2P = IDM2+1
      JDM2P = JDM2+1

      ALLOCATE(FSUM(IDM2,JDM2))
      ALLOCATE(ALON(IDM2P), ALAT(JDM2P), BLAT(JDM2))
      ALLOCATE(IR1(IDM2), IR2(IDM2), JR1(JDM2), JR2(JDM2))

      !==============================================================
      ! Set up input data lon, lat of each grid box edges in radians
      !==============================================================
      DRLON = 2d0 * PI / DBLE(IDM1)
      DRLAT =       PI / DBLE(JDM1)
      DO I=1, IDM1P
         RLON(I) = DRLON * ( I - 1 )
      ENDDO

      RLAT(1) = PI / 2d0
      DO J=2, JDM1
         RLAT(J) = PI / 2d0 - DRLAT * ( J - 1 )
      ENDDO
      RLAT(JDM1P) = - ( PI / 2d0 )

      !==============================================================
      ! Set up model data lon, lat of each grid box edges in radians
      !==============================================================
      DALON = 2d0 * PI / DBLE(IDM2)
      DO I=1, IDM2P
         ALON(I) = DALON * ( I - 1.5 )
      ENDDO

      CALL CNV_GAUSSIAN_LAT(JDM2, BLAT)
      ALAT(1) = PI / 2d0
      DO J=1, JDM2-1
         ALAT(J+1) = 0.5 * ( BLAT(J) + BLAT(J+1) )
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
               IR2(IA) = IR
               GO TO 120
            ELSE
               IR = IR + 1
               GO TO 90
            ENDIF
         ENDIF
100   CONTINUE

120   CONTINUE

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
      ! Average if it is concentration (FLAG = 2)
      !===========================================
      IF (FLAG .EQ. 2) OUTDATA(:,:) = OUTDATA(:,:) / FSUM(:,:)

      DEALLOCATE(ALON, ALAT, BLAT)
      DEALLOCATE(FSUM)
      DEALLOCATE(IR1, IR2, JR1, JR2)

      END SUBROUTINE REGRID_1x1G
