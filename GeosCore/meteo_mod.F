#include <define.h>
!=====================================================================
!
!     MODULE METEO_MOD
!
!        This module accumulates and saves GRIMs meteorology for
!        offline GRIMs/GEOS-Chem.
!
!                                                       Seungun Lee
!
!=====================================================================

      MODULE METEO_MOD

      USE DAO_MOD
      USE TIME_MOD, ONLY : GET_TIME_AHEAD

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      ! Parameter
      INTEGER, PARAMETER :: N_A3_FIELDS = 18
      INTEGER, PARAMETER :: N_A6_FIELDS = 10
      INTEGER, PARAMETER :: N_I6_FIELDS = 3
      INTEGER, PARAMETER :: NI_NEST=101
      INTEGER, PARAMETER :: NJ_NEST=101
      REAL*8,  PARAMETER :: DI_NEST=0.1d0
      REAL*8,  PARAMETER :: DJ_NEST=0.1d0
      REAL*8,  PARAMETER :: I1_NEST=122d0
      REAL*8,  PARAMETER :: J1_NEST=32d0

      ! Accumulated A3, A6
      REAL*8, ALLOCATABLE :: A3_FIELDS(:,:,:)
      REAL*8, ALLOCATABLE :: A6_FIELDS(:,:,:,:)

      ! Count
      INTEGER :: A3_COUNT, A6_COUNT

      CONTAINS

!=====================================================================

      SUBROUTINE DO_METEO

      INTEGER :: L

      ! Accumulate A3 fields
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

      ! Accumulate A6 fields
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

      ! Count
      A3_COUNT = A3_COUNT + 1
      A6_COUNT = A6_COUNT + 1

      END SUBROUTINE DO_METEO

!=====================================================================

      SUBROUTINE SAVE_A3

      INTEGER :: DATE(2), N
      REAL*8  :: P_ARRAY(IIPAR,JJPAR,N_A3_FIELDS)
      REAL*8  :: F_ARRAY(IGLOB,JGLOB,N_A3_FIELDS)
      REAL*8  :: ARRAY(I2x25,J2x25,N_A3_FIELDS)
      CHARACTER(LEN=8) :: IDENT, NAME(N_A3_FIELDS)
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''

      IF ( A3_COUNT .EQ. 0 ) RETURN

      NAME = (/ 'ALBEDO', 'CLDFRC', 'GWETTOP', 'HFLUX',
     &          'PARDF', 'PARDR', 'PBL', 'PREACC', 'PRECON',
     &          'RADLWG', 'RADSWG', 'SNOW', 'TS', 'TSKIN',
     &          'U10M', 'USTAR', 'V10M', 'Z0' /)

      ! For timestamp
      DATE = GET_TIME_AHEAD(-90)

      ! Filename
      WRITE(FILENAME, '(i8,".a3.2x25")') DATE(1)

      ! Open file if it is a new one
      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (LMASTER) THEN
         CLOSE(3001)
         OPEN(3001, FILE=TRIM(FILENAME), STATUS='REPLACE',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 22",1x,i2.2)') N_A3_FIELDS
         WRITE(3001) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      ! Divide by accumulated count
      P_ARRAY = A3_FIELDS / DBLE( A3_COUNT )

      ! Transform
      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, N_A3_FIELDS )

      IF (LMASTER) THEN
      ! Regrid
      DO N = 1, N_A3_FIELDS
         CALL REGRID_METEO( F_ARRAY(:,:,N), ARRAY(:,:,N) )
      ENDDO

      ! Write
      DO N = 1, N_A3_FIELDS
         WRITE(3001) NAME(N)
         WRITE(3001) DATE(1), DATE(2), REAL(ARRAY(:,:,N),4)
      ENDDO
      ENDIF

      END SUBROUTINE SAVE_A3

!=====================================================================

      SUBROUTINE SAVE_A6

      INTEGER :: DATE(2), L, N
      REAL*8  :: P_ARRAY(IIPAR,JJPAR,LLPAR,N_A6_FIELDS)
      REAL*8  :: F_ARRAY(IGLOB,JGLOB,LLPAR,N_A6_FIELDS)
      REAL*8  :: ARRAY(I2x25,J2x25,LLPAR,N_A6_FIELDS)
      CHARACTER(LEN=8) :: IDENT, NAME(N_A6_FIELDS)
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''

      IF ( A6_COUNT .EQ. 0 ) RETURN

      NAME = (/ 'CLDTOT', 'CLDMAS', 'DTRAIN', 'MOISTQ',
     &          'OPTDEPTH', 'Q', 'T', 'U', 'V', 'OMEGA' /)

      ! For timestamp
      DATE = GET_TIME_AHEAD(-180)

      ! Filename
      WRITE(FILENAME, '(i8,".a6.2x25")') DATE(1)

      ! Open file if it is a new one
      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (LMASTER) THEN
         CLOSE(3002)
         OPEN(3002, FILE=TRIM(FILENAME), STATUS='REPLACE',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 22",1x,i2.2)') N_A6_FIELDS
         WRITE(3002) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      ! Divide by accumulated count
      P_ARRAY = A6_FIELDS / DBLE( A6_COUNT )

      ! Transform
      DO N = 1, N_A6_FIELDS
         CALL TRANS_OUTPUT( P_ARRAY(:,:,:,N), F_ARRAY(:,:,:,N), LLPAR )
      ENDDO

      IF (LMASTER) THEN
      ! Regrid
      DO N = 1, N_A6_FIELDS
      DO L = 1, LLPAR
         CALL REGRID_METEO( F_ARRAY(:,:,L,N), ARRAY(:,:,L,N) )
      ENDDO
      ENDDO

      ! Write
      DO N = 1, N_A6_FIELDS
         WRITE(3002) NAME(N)
         WRITE(3002) DATE(1), DATE(2), REAL(ARRAY(:,:,:,N),4)
      ENDDO
      ENDIF

      END SUBROUTINE SAVE_A6

!=====================================================================

      SUBROUTINE SAVE_I6

      INTEGER :: DATE(2), N
      REAL*8  :: P_ARRAY(IIPAR,JJPAR,N_I6_FIELDS)
      REAL*8  :: F_ARRAY(IGLOB,JGLOB,N_I6_FIELDS)
      REAL*8  :: ARRAY(I2x25,J2x25,N_I6_FIELDS)
      CHARACTER(LEN=8) :: IDENT, NAME(N_I6_FIELDS)
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''

      NAME = (/ 'LWI', 'PS', 'TROPP' /)

      ! For timestamp
      DATE = GET_TIME_AHEAD(0)

      ! Filename
      WRITE(FILENAME, '(i8,".i6.2x25")') DATE(1)

      ! Open file if it is a new one
      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (LMASTER) THEN
         CLOSE(3003)
         OPEN(3003, FILE=TRIM(FILENAME), STATUS='REPLACE',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 22",1x,i2.2)') N_I6_FIELDS
         WRITE(3003) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      ! Get meteorology
      P_ARRAY(:,:,1) = LWI
      P_ARRAY(:,:,2) = PSC2
      P_ARRAY(:,:,3) = TROPP

      ! Transform
      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, N_I6_FIELDS )

      IF (LMASTER) THEN
      ! Regrid
      DO N = 1, N_I6_FIELDS
         CALL REGRID_METEO( F_ARRAY(:,:,N), ARRAY(:,:,N) )
      ENDDO

      ! Write
      DO N = 1, N_I6_FIELDS
         WRITE(3003) NAME(N)
         WRITE(3003) DATE(1), DATE(2), REAL(ARRAY(:,:,N),4)
      ENDDO
      ENDIF

      END SUBROUTINE SAVE_I6

!=====================================================================

      SUBROUTINE SAVE_A3_NESTED

      INTEGER :: DATE(2), N
      REAL*8  :: P_ARRAY(IIPAR,JJPAR,N_A3_FIELDS)
      REAL*8  :: F_ARRAY(IGLOB,JGLOB,N_A3_FIELDS)
      REAL*8  :: ARRAY(NI_NEST,NJ_NEST,N_A3_FIELDS)
      CHARACTER(LEN=8) :: IDENT, NAME(N_A3_FIELDS)
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''

      IF ( A3_COUNT .EQ. 0 ) RETURN

      NAME = (/ 'ALBEDO', 'CLDFRC', 'GWETTOP', 'HFLUX',
     &          'PARDF', 'PARDR', 'PBL', 'PREACC', 'PRECON',
     &          'RADLWG', 'RADSWG', 'SNOW', 'TS', 'TSKIN',
     &          'U10M', 'USTAR', 'V10M', 'Z0' /)

      ! For timestamp
      DATE = GET_TIME_AHEAD(-90)

      ! Filename
      WRITE(FILENAME, '(i8,".a3.01x01")') DATE(1)

      ! Open file if it is a new one
      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (LMASTER) THEN
         CLOSE(4001)
         OPEN(4001, FILE=TRIM(FILENAME), STATUS='REPLACE',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 11",1x,i2.2)') N_A3_FIELDS
         WRITE(4001) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      ! Divide by accumulated count
      P_ARRAY = A3_FIELDS / DBLE( A3_COUNT )

      ! Transform
      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, N_A3_FIELDS )

      IF (LMASTER) THEN
      ! Regrid
      DO N = 1, N_A3_FIELDS
         CALL REGRID_METEO_NESTED( F_ARRAY(:,:,N), ARRAY(:,:,N) )
      ENDDO

      ! Write
      DO N = 1, N_A3_FIELDS
         WRITE(4001) NAME(N)
         WRITE(4001) DATE(1), DATE(2), REAL(ARRAY(:,:,N),4)
      ENDDO
      ENDIF

      END SUBROUTINE SAVE_A3_NESTED

!=====================================================================

      SUBROUTINE SAVE_A6_NESTED

      INTEGER :: DATE(2), L, N
      REAL*8  :: P_ARRAY(IIPAR,JJPAR,LLPAR,N_A6_FIELDS)
      REAL*8  :: F_ARRAY(IGLOB,JGLOB,LLPAR,N_A6_FIELDS)
      REAL*8  :: ARRAY(NI_NEST,NJ_NEST,LLPAR,N_A6_FIELDS)
      CHARACTER(LEN=8) :: IDENT, NAME(N_A6_FIELDS)
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''

      IF ( A6_COUNT .EQ. 0 ) RETURN

      NAME = (/ 'CLDTOT', 'CLDMAS', 'DTRAIN', 'MOISTQ',
     &          'OPTDEPTH', 'Q', 'T', 'U', 'V', 'OMEGA' /)

      ! For timestamp
      DATE = GET_TIME_AHEAD(-180)

      ! Filename
      WRITE(FILENAME, '(i8,".a6.01x01")') DATE(1)

      ! Open file if it is a new one
      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (LMASTER) THEN
         CLOSE(4002)
         OPEN(4002, FILE=TRIM(FILENAME), STATUS='REPLACE',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 11",1x,i2.2)') N_A6_FIELDS
         WRITE(4002) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      ! Divide by accumulated count
      P_ARRAY = A6_FIELDS / DBLE( A6_COUNT )

      ! Transform
      DO N = 1, N_A6_FIELDS
         CALL TRANS_OUTPUT( P_ARRAY(:,:,:,N), F_ARRAY(:,:,:,N), LLPAR )
      ENDDO

      IF (LMASTER) THEN
      ! Regrid
      DO N = 1, N_A6_FIELDS
      DO L = 1, LLPAR
         CALL REGRID_METEO_NESTED( F_ARRAY(:,:,L,N), ARRAY(:,:,L,N) )
      ENDDO
      ENDDO

      ! Write
      DO N = 1, N_A6_FIELDS
         WRITE(4002) NAME(N)
         WRITE(4002) DATE(1), DATE(2), REAL(ARRAY(:,:,:,N),4)
      ENDDO
      ENDIF

      END SUBROUTINE SAVE_A6_NESTED

!=====================================================================

      SUBROUTINE SAVE_I6_NESTED

      INTEGER :: DATE(2), N
      REAL*8  :: P_ARRAY(IIPAR,JJPAR,N_I6_FIELDS)
      REAL*8  :: F_ARRAY(IGLOB,JGLOB,N_I6_FIELDS)
      REAL*8  :: ARRAY(NI_NEST,NJ_NEST,N_I6_FIELDS)
      CHARACTER(LEN=8) :: IDENT, NAME(N_I6_FIELDS)
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255), SAVE :: FILENAME_SAVE=''

      NAME = (/ 'LWI', 'PS', 'TROPP' /)

      ! For timestamp
      DATE = GET_TIME_AHEAD(0)

      ! Filename
      WRITE(FILENAME, '(i8,".i6.01x01")') DATE(1)

      ! Open file if it is a new one
      IF (FILENAME .NE. FILENAME_SAVE) THEN

         IF (LMASTER) THEN
         CLOSE(4003)
         OPEN(4003, FILE=TRIM(FILENAME), STATUS='REPLACE',
     &              ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         WRITE(IDENT, '("GR 11",1x,i2.2)') N_I6_FIELDS
         WRITE(4003) IDENT
         ENDIF

         FILENAME_SAVE = FILENAME

      ENDIF

      ! Get meteorology
      P_ARRAY(:,:,1) = LWI
      P_ARRAY(:,:,2) = PSC2
      P_ARRAY(:,:,3) = TROPP

      ! Transform
      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, N_I6_FIELDS )

      IF (LMASTER) THEN
      ! Regrid
      DO N = 1, N_I6_FIELDS
         CALL REGRID_METEO_NESTED( F_ARRAY(:,:,N), ARRAY(:,:,N) )
      ENDDO

      ! Write
      DO N = 1, N_I6_FIELDS
         WRITE(4003) NAME(N)
         WRITE(4003) DATE(1), DATE(2), REAL(ARRAY(:,:,N),4)
      ENDDO
      ENDIF

      END SUBROUTINE SAVE_I6_NESTED

!=====================================================================

      SUBROUTINE REGRID_METEO( INDATA, OUTDATA )

      USE REGRID_A2A_MOD, ONLY : MAP_A2A
#ifndef RMP
      USE comfgrid,       ONLY : colrad
#endif

      USE CMN_GCTM_MOD

      INTEGER :: IM, JM, IN, JN
      INTEGER :: I, J, L, N
      REAL*8  :: INDATA(IGLOB,JGLOB)
      REAL*8  :: OUTDATA(I2x25,J2x25)
      REAL*8, ALLOCATABLE, DIMENSION(:)   :: LON1, LON2, SIN1, SIN2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Q1, Q2

#ifndef RMP
      ! Set dimension
      IM = IGLOB
      JM = JGLOB

      ! Allocate variables
      ALLOCATE( LON1(IM+1) )
      ALLOCATE( SIN1(JM+1) )
      ALLOCATE( Q1(IM,JM) )

      ! Longitude
      DO I = 1, IM+1
         LON1(I) = -180d0 + DISIZE * ( DBLE(I) - 1.5d0 )
      ENDDO

      ! Latitude
      SIN1(1) = -1d0
      DO J = 2, JM/2
         SIN1(J) = -1d0 * COS( ( colrad(J-1) + colrad(J) ) / 2d0 )
      ENDDO
      SIN1(JM/2+1) = 0d0
      DO J = JM/2+2, JM+1
         SIN1(J) = -1d0 * SIN1(JM+2-J)
      ENDDO

      ! Set dimension
      IN = I2x25
      JN = J2x25

      ! Allocate variables
      ALLOCATE( LON2(IN+1) )
      ALLOCATE( SIN2(JN+1) )
      ALLOCATE( Q2(IN,JN) )

      ! Longitude
      DO I = 1, IN+1
         LON2(I) = -180d0 + 2.5d0 * ( DBLE(I) - 1.5d0 )
      ENDDO

      ! Latitude
      SIN2(1) = -1d0
      DO J = 2, JN
         SIN2(J) = SIND( -90d0 + 2d0 * ( DBLE(J) - 1.5d0 ) )
      ENDDO
      SIN2(JN+1) = 1d0

      ! Regrid
      Q1 = INDATA
      CALL MAP_A2A( IM, JM, LON1, SIN1, Q1,
     &              IN, JN, LON2, SIN2, Q2, 0, 1 )
      OUTDATA = Q2

      ! Deallocate variables
      DEALLOCATE( LON1, LON2, SIN1, SIN2, Q1, Q2 )
#endif

      END SUBROUTINE REGRID_METEO

!=====================================================================

      SUBROUTINE REGRID_METEO_NESTED( INDATA, OUTDATA )

      USE REGRID_A2A_MOD, ONLY : MAP_A2A
#ifndef RMP
      USE comfgrid,       ONLY : colrad
#endif

      USE CMN_GCTM_MOD

      INTEGER :: IM, JM, IN, JN
      INTEGER :: I, J, L, N
      REAL*8  :: INDATA(IGLOB,JGLOB)
      REAL*8  :: OUTDATA(NI_NEST,NJ_NEST)
      REAL*8, ALLOCATABLE, DIMENSION(:)   :: LON1, LON2, SIN1, SIN2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Q1, Q2
#ifdef RMP
      REAL*8  :: flat(IGLOB,JGLOB), flon(IGLOB,JGLOB)
      REAL*8  :: delx, dely, dlamda0
#endif

      ! Set dimension
      IM = IGLOB
      JM = JGLOB

      ! Allocate variables
      ALLOCATE( LON1(IM+1) )
      ALLOCATE( SIN1(JM+1) )
      ALLOCATE( Q1(IM,JM) )

#ifndef RMP
      ! Longitude
      DO I = 1, IM+1
         LON1(I) = -180d0 + DISIZE * ( DBLE(I) - 1.5d0 )
      ENDDO

      ! Latitude
      SIN1(1) = -1d0
      DO J = 2, JM/2
         SIN1(J) = -1d0 * COS( ( colrad(J-1) + colrad(J) ) / 2d0 )
      ENDDO
      SIN1(JM/2+1) = 0d0
      DO J = JM/2+2, JM+1
         SIN1(J) = -1d0 * SIN1(JM+2-J)
      ENDDO
#else
      CALL rmp_setup_rsm_grid(flat,flon,delx,dely,dlamda0)

      ! Longitude
      LON1(1) = flon(1,1) - ( flon(2,1) - flon(1,1) ) * 0.5d0
      DO I = 1, IM-1
         LON1(I+1) = ( flon(I,1) + flon(I+1,1) ) * 0.5d0
      ENDDO
      LON1(IM+1) = flon(IM,1) +
     &                    ( flon(IM,1) - flon(IM-1,1) ) * 0.5d0

      DO I = 1, IM+1
         LON1(I) = LON1(I) * 180d0 / PI
      ENDDO

      ! Latitude
      SIN1(1) = flat(1,1) - ( flat(1,3) - flat(1,1) ) * 0.5d0
      SIN1(JM+1) = flat(1,2) + ( flat(1,2) - flat(1,4) ) * 0.5d0
      SIN1(JM/2+1) = ( flat(1,JM/2) + flat(1,JM/2+1) ) * 0.5d0
      DO J = 1, JM/2-1
         SIN1(J+1) = ( flat(1,2*J-1) + flat(1,2*J+1) ) * 0.5d0
         SIN1(JM+1-J) = ( flat(1,2*J) + flat(1,2*J+2) ) * 0.5d0
      ENDDO

      DO J = 1, JM+1
         SIN1(J) = SIN( SIN1(J) )
      ENDDO
#endif

      ! Set dimension
      IN = NI_NEST
      JN = NJ_NEST

      ! Allocate variables
      ALLOCATE( LON2(IN+1) )
      ALLOCATE( SIN2(JN+1) )
      ALLOCATE( Q2(IN,JN) )

      ! Longitude
      DO I = 1, IN+1
         LON2(I) = I1_NEST + DI_NEST * ( DBLE(I) - 1.5d0 )
      ENDDO

      ! Latitude
      DO J = 1, JN+1
         SIN2(J) = SIND( J1_NEST + DJ_NEST * ( DBLE(J) - 1.5d0 ) )
      ENDDO

      ! Regrid
      Q1 = INDATA
      CALL MAP_A2A( IM, JM, LON1, SIN1, Q1,
     &              IN, JN, LON2, SIN2, Q2, 0, 1 )
      OUTDATA = Q2

      ! Deallocate variables
      DEALLOCATE( LON1, LON2, SIN1, SIN2, Q1, Q2 )

      END SUBROUTINE REGRID_METEO_NESTED

!=====================================================================

      SUBROUTINE INIT_A3

      ! Initialize A3
      A3_FIELDS = 0d0
      A3_COUNT  = 0

      END SUBROUTINE INIT_A3

!=====================================================================

      SUBROUTINE INIT_A6

      ! Initialize A6
      A6_FIELDS = 0d0
      A6_COUNT  = 0

      END SUBROUTINE INIT_A6

!=====================================================================

      SUBROUTINE INIT_METEO

      ! Allocate module variables
      ALLOCATE( A3_FIELDS(IIPAR,JJPAR,N_A3_FIELDS) )
      ALLOCATE( A6_FIELDS(IIPAR,JJPAR,LLPAR,N_A6_FIELDS) )

      ! Initialize module variables
      A3_FIELDS = 0d0
      A6_FIELDS = 0d0
      A3_COUNT  = 0
      A6_COUNT  = 0

      END SUBROUTINE INIT_METEO

!=====================================================================

      END MODULE METEO_MOD
