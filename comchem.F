!=====================================================================
!
!     MODULE COMCHEM
!
!        This module converts and saves the variables which
!        are needed in the chemistry processes.
!
!                                            written by Seungun Lee
!
!=====================================================================

      MODULE COMCHEM

      USE COMFCST,     ONLY : FCD, FCU
      USE COMFPHYS,    ONLY : U10M_ => U10M
      USE COMFPHYS,    ONLY : V10M_ => V10M
      USE COMFPHYS,    ONLY : HPBL, PSURF
      USE COMFSPEC_VR, ONLY : TEM
      USE COMFVER,     ONLY : DELTIM
      USE COMGRAD,     ONLY : ALNDR, ALVDR
      USE COMSFC,      ONLY : SMC, SLMSK, SNOWEQ, TSEA, UUSTAR, Z0CM
      USE DAO_MOD,     ONLY : ALBD, CLDFRC, CMFMC, DTRAIN, GWETTOP
      USE DAO_MOD,     ONLY : HFLUX, LWI, MOISTQ, PBL, PREACC, PRECON
      USE DAO_MOD,     ONLY : PS1, RADSWG, SNOMAS, SPHU, T, TROPP, TS
      USE DAO_MOD,     ONLY : U10M, USTAR, V10M, Z0
      USE NISLQ,       ONLY : SLQ_Q1
      USE PARAMODEL,   ONLY : LEVSS
      USE TRACER_MOD,  ONLY : N_TRACERS

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      REAL*8, DIMENSION(:,:),   ALLOCATABLE :: _CLDSA
      REAL*8, DIMENSION(:,:),   ALLOCATABLE :: _DTSFC
      REAL*8, DIMENSION(:,:),   ALLOCATABLE :: _RAINCPS
      REAL*8, DIMENSION(:,:),   ALLOCATABLE :: _RAINTOT
      REAL*8, DIMENSION(:,:),   ALLOCATABLE :: _SSWDN
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: _MOISTQ1
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: _MOISTQ2
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: _DETRN
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: DP1
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: DP2
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: DP3

      INTEGER :: LEVH_, NCLD, NLEVS, NLEVSP

      CONTAINS

!==========================================================

      SUBROUTINE UPDATE_COMCHEM

      REAL*8  :: WORK(IIPAR,LLPAR,JJPAR)
      REAL*8  :: WORK1(IIPAR,JJPAR,LLPAR)
      REAL*8  :: WORK2(IIPAR,JJPAR,LLPAR)
      INTEGER :: L

      ! Surface Albedo [unitless]
      ALBD = 0.466d0 * ALVDR + 0.529d0 * ALNDR

      ! 2-D column cloud fraction [unitless]
      CLDFRC = _CLDSA

      ! Cloud mass flux [kg/m2/s]
      CALL TRANSFORM_3D( FCU, WORK1 )
      CALL TRANSFORM_3D( FCD, WORK2 )
      CMFMC = WORK1 + WORK2

      ! Detrainment rate
      DTRAIN = _DETRN

      ! Soil moisture [unitless]
      GWETTOP = SMC(:,:,1)

      ! Sensible heat flux [W/m2]
      HFLUX = _DTSFC

      ! Land/water flags [unitless]
      LWI = SLMSK

      ! Tendency of specific humidity [kg/kg/s]
      CALL TRANSFORM_3D( _MOISTQ1, WORK1 )
      CALL TRANSFORM_3D( _MOISTQ2, WORK2 )
      DO L=1, LLPAR
         MOISTQ(L,:,:) = WORK2(:,:,L) - WORK1(:,:,L)
      ENDDO
      MOISTQ = - MOISTQ / ( 2 * DELTIM )

      ! Mixed layer depth [m]
      PBL = HPBL

      ! Total precip at ground [mm/day]
      PREACC = _RAINTOT / DELTIM * 1d3 * 86400d0

      ! Convective precip at ground [mm/day]
      PRECON = _RAINCPS / DELTIM * 1d3 * 86400d0

      ! Surface pressure [hPa]
      PS1 = PSURF * 1d1

      ! Shortwave rad at ground [W/m2]
      RADSWG = _SSWDN

      ! Water equivalent snow depth [mm H2O]
      SNOMAS = SNOWEQ

      ! Specific humidity [g/kg]
      CALL TRANSFORM_3D( SLQ_Q1(:,1:LLPAR,:), SPHU )
      SPHU = SPHU * 1d3

      ! Temperature [K]
      CALL SPH_FFT_DRIVER(-1,WORK,TEM,1)
      CALL TRANSFORM_3D( WORK, T )

      ! Tropopause [hPa]
      CALL TPP_REICHLER( T, TROPP )
      TROPP = TROPP * 1d-2

      ! Surface temperature [K]
      TS = TSEA

      ! 10m U-winds [m/s]
      U10M = U10M_

      ! Friction velocity [m/s]
      USTAR = UUSTAR

      ! 10m V-winds [m/s]
      V10M = V10M_

      ! Roughness height [m]
      Z0 = Z0CM * 1d-2

      END SUBROUTINE UPDATE_COMCHEM

!==========================================================

      SUBROUTINE TRANSFORM_3D( INDATA, OUTDATA )

      REAL*8, INTENT(IN)  :: INDATA(IIPAR,LLPAR,JJPAR)
      REAL*8, INTENT(OUT) :: OUTDATA(IIPAR,JJPAR,LLPAR)
      INTEGER :: L

      DO L=1, LLPAR
         OUTDATA(:,:,L) = INDATA(:,L,:)
      ENDDO

      END SUBROUTINE TRANSFORM_3D

!=====================================================================

      SUBROUTINE INIT_COMCHEM

      ALLOCATE( _CLDSA   (IIPAR,JJPAR) )
      ALLOCATE( _DTSFC   (IIPAR,JJPAR) )
      ALLOCATE( _RAINCPS (IIPAR,JJPAR) )
      ALLOCATE( _RAINTOT (IIPAR,JJPAR) )
      ALLOCATE( _SSWDN   (IIPAR,JJPAR) )
      ALLOCATE( _MOISTQ1 (IIPAR,LLPAR,JJPAR) )
      ALLOCATE( _MOISTQ2 (IIPAR,LLPAR,JJPAR) )
      ALLOCATE( _DETRN   (IIPAR,JJPAR,LLPAR) )

      ALLOCATE( DP1      (IIPAR,LLPAR,JJPAR) )
      ALLOCATE( DP2      (IIPAR,LLPAR,JJPAR) )
      ALLOCATE( DP3      (IIPAR,LLPAR,JJPAR) )

      _CLDSA     = 0
      _DTSFC     = 0
      _RAINCPS   = 0
      _RAINTOT   = 0
      _SSWDN     = 0
      _MOISTQ1   = 0
      _MOISTQ2   = 0
      _DETRN     = 0

      DP1 = 0d0
      DP2 = 0d0
      DP3 = 0d0

      LEVH_  = N_TRACERS * LLPAR
      NCLD   = N_TRACERS
      NLEVS  = N_TRACERS * LLPAR
      NLEVSP = N_TRACERS * LEVSS

      END SUBROUTINE INIT_COMCHEM

!==========================================================

      END MODULE COMCHEM
