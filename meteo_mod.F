      MODULE METEO_MOD

      USE DAO_MOD
      USE MODULE_TRANS, ONLY : DYN_TRANS2OUTPUT_GRID
      USE TIME_MOD,     ONLY : GET_TIME_AHEAD
      USE CMN_SIZE_MOD

#     include "define.h"

      ! Parameter
      INTEGER, PARAMETER :: N_A3_FIELDS = 18
      INTEGER, PARAMETER :: N_I6_FIELDS = 3
      INTEGER, PARAMETER :: N_A6_FIELDS = 10

      ! A3, I6, A6
      REAL*8, ALLOCATABLE :: A3_FIELDS(:,:,:)
      REAL*8, ALLOCATABLE :: I6_FIELDS(:,:,:)
      REAL*8, ALLOCATABLE :: A6_FIELDS(:,:,:,:)

      ! Count
      INTEGER :: A3_COUNT, A6_COUNT

      CONTAINS

!------------------------------------------------------------------------

      SUBROUTINE DO_METEO

      IMPLICIT NONE

      INTEGER :: L

      A3_FIELDS(:,:,1)  = A3_FIELDS(:,:,1)  + ALBD
      A3_FIELDS(:,:,2)  = A3_FIELDS(:,:,2)  + CLDFRC
      A3_FIELDS(:,:,3)  = A3_FIELDS(:,:,3)  + GWETTOP
      A3_FIELDS(:,:,4)  = A3_FIELDS(:,:,4)  + HFLUX
      A3_FIELDS(:,:,5)  = A3_FIELDS(:,:,5)  + PARDF
      A3_FIELDS(:,:,6)  = A3_FIELDS(:,:,6)  + PARDR
      A3_FIELDS(:,:,7)  = A3_FIELDS(:,:,7)  + PBL
      A3_FIELDS(:,:,8)  = A3_FIELDS(:,:,8)  + PREACC
      A3_FIELDS(:,:,9)  = A3_FIELDS(:,:,9)  + PRECON
      A3_FIELDS(:,:,10) = A3_FIELDS(:,:,10) + RADLWG
      A3_FIELDS(:,:,11) = A3_FIELDS(:,:,11) + RADSWG
      A3_FIELDS(:,:,12) = A3_FIELDS(:,:,12) + SNOW
      A3_FIELDS(:,:,13) = A3_FIELDS(:,:,13) + TS
      A3_FIELDS(:,:,14) = A3_FIELDS(:,:,14) + TSKIN
      A3_FIELDS(:,:,15) = A3_FIELDS(:,:,15) + U10M
      A3_FIELDS(:,:,16) = A3_FIELDS(:,:,16) + USTAR
      A3_FIELDS(:,:,17) = A3_FIELDS(:,:,17) + V10M
      A3_FIELDS(:,:,18) = A3_FIELDS(:,:,18) + Z0

      DO L = 1, LLPAR
         A6_FIELDS(:,:,L,1)  = A6_FIELDS(:,:,L,1)  + CLDF(L,:,:)
         A6_FIELDS(:,:,L,2)  = A6_FIELDS(:,:,L,2)  + CLDMAS(:,:,L)
         A6_FIELDS(:,:,L,3)  = A6_FIELDS(:,:,L,3)  + DTRAIN(:,:,L)
         A6_FIELDS(:,:,L,4)  = A6_FIELDS(:,:,L,4)  + MOISTQ(L,:,:)
         A6_FIELDS(:,:,L,5)  = A6_FIELDS(:,:,L,5)  + OPTDEP(L,:,:)
         A6_FIELDS(:,:,L,6)  = A6_FIELDS(:,:,L,6)  + SPHU(:,:,L)
         A6_FIELDS(:,:,L,7)  = A6_FIELDS(:,:,L,7)  + T(:,:,L)
         A6_FIELDS(:,:,L,8)  = A6_FIELDS(:,:,L,8)  + UWND(:,:,L)
         A6_FIELDS(:,:,L,9)  = A6_FIELDS(:,:,L,9)  + VWND(:,:,L)
         A6_FIELDS(:,:,L,10) = A6_FIELDS(:,:,L,10) + OMEGA(:,:,L)
      ENDDO

      A3_COUNT = A3_COUNT + 1
      A6_COUNT = A6_COUNT + 1

      END SUBROUTINE DO_METEO

!------------------------------------------------------------------------

      SUBROUTINE SAVE_A3

      IMPLICIT NONE

      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''
      CHARACTER(LEN=8) :: IDENT, NAME(N_A3_FIELDS)
      REAL*8 :: F_ARRAY(IGLOB*JGLOB,N_A3_FIELDS)
      REAL*8 :: ARRAY8(I2x25,J2x25,N_A3_FIELDS)
      REAL*4 :: ARRAY4(I2x25,J2x25,N_A3_FIELDS)
      INTEGER :: DATE(2), I, J, N

      NAME = (/ 'ALBEDO', 'CLDFRC', 'GWETTOP', 'HFLUX',
     &          'PARDF', 'PARDR', 'PBL', 'PREACC', 'PRECON',
     &          'RADLWG', 'RADSWG', 'SNOW', 'TS', 'TSKIN',
     &          'U10M', 'USTAR', 'V10M', 'Z0' /)

      DATE = GET_TIME_AHEAD(-90)

      WRITE(FILENAME, '(i8,".a3.2x25")') DATE(1)

      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (MASTER) THEN
         CLOSE(3001)
         OPEN(3001, FILE=TRIM(FILENAME), STATUS='NEW',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 22",1x,i2.2)') N_A3_FIELDS
         WRITE(3001) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      A3_FIELDS = A3_FIELDS / A3_COUNT

#ifdef MP
      CALL MPGP2F(A3_FIELDS,IIPAR,JJPAR,
     &            F_ARRAY,IGLOB*2,JGLOB/2,N_A3_FIELDS)
#else
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         F_ARRAY(IIPAR*(J-1)+I,:) = A3_FIELDS(I,J,:)
      ENDDO
      ENDDO
#endif

      IF (MASTER) THEN
      CALL DYN_TRANS2OUTPUT_GRID(F_ARRAY,N_A3_FIELDS)

      CALL XSHIFT(F_ARRAY,IGLOB,JGLOB,N_A3_FIELDS)
      CALL YREV(F_ARRAY,IGLOB,JGLOB,N_A3_FIELDS)

      DO N = 1, N_A3_FIELDS
         CALL REGRID_METEO( F_ARRAY(:,N), ARRAY8(:,:,N) )
      ENDDO

      ARRAY4 = ARRAY8

      DO N = 1, N_A3_FIELDS
         WRITE(3001) NAME(N)
         WRITE(3001) DATE(1), DATE(2), ARRAY4(:,:,N)
      ENDDO
      ENDIF

      A3_FIELDS = 0d0
      A3_COUNT = 0

      END SUBROUTINE SAVE_A3

!------------------------------------------------------------------------

      SUBROUTINE SAVE_I6

      IMPLICIT NONE

      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''
      CHARACTER(LEN=8) :: IDENT, NAME(N_I6_FIELDS)
      REAL*8 :: F_ARRAY(IGLOB*JGLOB,N_I6_FIELDS)
      REAL*8 :: ARRAY8(I2x25,J2x25,N_I6_FIELDS)
      REAL*4 :: ARRAY4(I2x25,J2x25,N_I6_FIELDS)
      INTEGER :: DATE(2), I, J, N

      NAME = (/ 'LWI', 'PS', 'TROPP' /)

      DATE = GET_TIME_AHEAD(0)

      WRITE(FILENAME, '(i8,".i6.2x25")') DATE(1)

      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (MASTER) THEN
         CLOSE(3002)
         OPEN(3002, FILE=TRIM(FILENAME), STATUS='NEW',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 22",1x,i2.2)') N_I6_FIELDS
         WRITE(3002) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      I6_FIELDS(:,:,1) = LWI
      I6_FIELDS(:,:,2) = PS1
      I6_FIELDS(:,:,3) = TROPP

#ifdef MP
      CALL MPGP2F(I6_FIELDS,IIPAR,JJPAR,
     &            F_ARRAY,IGLOB*2,JGLOB/2,N_I6_FIELDS)
#else
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         F_ARRAY(IIPAR*(J-1)+I,:) = I6_FIELDS(I,J,:)
      ENDDO
      ENDDO
#endif

      IF (MASTER) THEN
      CALL DYN_TRANS2OUTPUT_GRID(F_ARRAY,N_I6_FIELDS)

      CALL XSHIFT(F_ARRAY,IGLOB,JGLOB,N_I6_FIELDS)
      CALL YREV(F_ARRAY,IGLOB,JGLOB,N_I6_FIELDS)

      DO N = 1, N_I6_FIELDS
         CALL REGRID_METEO( F_ARRAY(:,N), ARRAY8(:,:,N) )
      ENDDO

      ARRAY4 = ARRAY8

      DO N = 1, N_I6_FIELDS
         WRITE(3002) NAME(N)
         WRITE(3002) DATE(1), DATE(2), ARRAY4(:,:,N)
      ENDDO
      ENDIF

      END SUBROUTINE SAVE_I6

!------------------------------------------------------------------------

      SUBROUTINE SAVE_A6

      IMPLICIT NONE

      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''
      CHARACTER(LEN=8) :: IDENT, NAME(N_A6_FIELDS)
      REAL*8 :: F_ARRAY(IGLOB*JGLOB,LLPAR,N_A6_FIELDS)
      REAL*8 :: ARRAY8(I2x25,J2x25,LLPAR,N_A6_FIELDS)
      REAL*4 :: ARRAY4(I2x25,J2x25,LLPAR,N_A6_FIELDS)
      INTEGER :: DATE(2), I, J, L, N

      NAME = (/ 'CLDTOT', 'CLDMAS', 'DTRAIN', 'MOISTQ',
     &          'OPTDEPTH', 'Q', 'T', 'U', 'V', 'OMEGA' /)

      DATE = GET_TIME_AHEAD(-180)

      WRITE(FILENAME, '(i8,".a6.2x25")') DATE(1)

      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (MASTER) THEN
         CLOSE(3003)
         OPEN(3003, FILE=TRIM(FILENAME), STATUS='NEW',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 22",1x,i2.2)') N_A6_FIELDS
         WRITE(3003) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      A6_FIELDS = A6_FIELDS / A6_COUNT

#ifdef MP
      CALL MPGP2F(A6_FIELDS,IIPAR,JJPAR,
     &            F_ARRAY,IGLOB*2,JGLOB/2,LLPAR*N_A6_FIELDS)
#else
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         F_ARRAY(IIPAR*(J-1)+I,:,:) = A6_FIELDS(I,J,:,:)
      ENDDO
      ENDDO
#endif

      IF (MASTER) THEN
      CALL DYN_TRANS2OUTPUT_GRID(F_ARRAY,LLPAR*N_A6_FIELDS)

      CALL XSHIFT(F_ARRAY,IGLOB,JGLOB,LLPAR*N_A6_FIELDS)
      CALL YREV(F_ARRAY,IGLOB,JGLOB,LLPAR*N_A6_FIELDS)

      DO N = 1, N_A6_FIELDS
      DO L = 1, LLPAR
         CALL REGRID_METEO( F_ARRAY(:,L,N), ARRAY8(:,:,L,N) )
      ENDDO
      ENDDO

      ARRAY4 = ARRAY8

      DO N = 1, N_A6_FIELDS
         WRITE(3003) NAME(N)
         WRITE(3003) DATE(1), DATE(2), ARRAY4(:,:,:,N)
      ENDDO
      ENDIF

      A6_FIELDS = 0d0
      A6_COUNT = 0

      END SUBROUTINE SAVE_A6

!------------------------------------------------------------------------

      SUBROUTINE REGRID_METEO( INDATA, OUTDATA )

      USE REGRID_A2A_MOD, ONLY : MAP_A2A
      USE COMFGRID,       ONLY : COLRAD

      USE CMN_GCTM_MOD

      REAL*8 :: INDATA(IGLOB,JGLOB)
      REAL*8 :: OUTDATA(I2x25,J2x25)
      REAL*8, ALLOCATABLE, DIMENSION(:)   :: LON1, LON2, SIN1, SIN2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Q1, Q2
      INTEGER :: IM, JM, IN, JN
      INTEGER :: I, J, L, N

      IM = IGLOB
      JM = JGLOB

      ALLOCATE( LON1(IM+1) )
      ALLOCATE( SIN1(IM+1) )
      ALLOCATE( Q1(IM,JM) )

      LON1(1) = -180d0 - ( DISIZE / 2d0 )
      DO I = 2, IM+1
         LON1(I) = LON1(I-1) + DISIZE
      ENDDO

      SIN1(1) = -1d0
      DO J = 2, JM/2
         SIN1(J) = -1d0 * COS( ( COLRAD(J-1) + COLRAD(J) ) / 2d0 )
      ENDDO

      SIN1(JM/2+1) = 0d0
      DO J = JM/2+2, JM+1
         SIN1(J) = -1d0 * SIN1(JM+2-J)
      ENDDO

      IN = I2x25
      JN = J2x25

      ALLOCATE( LON2(IN+1) )
      ALLOCATE( SIN2(JN+1) )
      ALLOCATE( Q2(IN,JN) )

      LON2(1) = -180d0 - ( 2.5d0 / 2d0 )
      DO I = 2, IN+1
         LON2(I) = LON2(I-1) + 2.5d0
      ENDDO

      SIN2(1) = -1d0
      DO J = 2, JN
         SIN2(J) = SIND( -90d0 + 2d0 * (J-1) )
      ENDDO
      SIN2(JN+1) = 1d0

      Q1 = INDATA
      CALL MAP_A2A( IM, JM, LON1, SIN1, Q1,
     &              IN, JN, LON2, SIN2, Q2, 0, 1 )
      OUTDATA = Q2

      DEALLOCATE( LON1, LON2, SIN1, SIN2, Q1, Q2 )

      END SUBROUTINE REGRID_METEO

!------------------------------------------------------------------------

      SUBROUTINE INIT_METEO

      IMPLICIT NONE

      ALLOCATE( A3_FIELDS(IIPAR,JJPAR,N_A3_FIELDS) )
      A3_FIELDS = 0d0

      ALLOCATE( I6_FIELDS(IIPAR,JJPAR,N_I6_FIELDS) )
      I6_FIELDS = 0d0

      ALLOCATE( A6_FIELDS(IIPAR,JJPAR,LLPAR,N_A6_FIELDS) )
      A6_FIELDS = 0d0

      A3_COUNT = 0
      A6_COUNT = 0

      END SUBROUTINE INIT_METEO

!------------------------------------------------------------------------

      END MODULE METEO_MOD
