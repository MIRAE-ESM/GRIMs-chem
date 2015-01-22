! $Id: dust_shao_mod.f,v 1.3 2010/02/02 16:57:53 bmy Exp $
      MODULE DUST_SHAO_MOD
!
!  Module Variables:
!  ============================================================================
!  (10) TPT_FRZ_PNT      (REAL*8)  : Freezing point of water        [K       ]
!  (25) MSS_FRC_CLY      (REAL*8 ) : Mass fraction of clay          [fraction]
!  (26) MSS_FRC_SND      (REAL*8 ) : Mass fraction of sand          [fraction]
!(26-1) MSS_FRC_LOA      (REAL*8 ) : Mass fraction of loam          [fraction]
!(26-2) MSS_FRC_SLO      (REAL*8 ) : Mass fraction of sandy loam    [fraction]
!  (27) SFC_TYP          (INTEGER) : Surface type index (0..28)     [unitless]
!  (30) TPT_GND          (REAL*8 ) : Ground temperature             [K       ]
!  (31) TPT_SOI          (REAL*8 ) : Soil temperature               [K       ]
!  (33) VAI_DST          (REAL*8 ) : Vegetation area index          [m2/m2   ]
!  (44) OROGRAPHY        (INTEGER) : 0=ocean; 1=land; 2=ice         [unitless]
!
!  Module Routines:
!  ============================================================================
!  (1 ) QF_SHAO                       : Driver routine for dust mobilization
!  (2 ) SOI_TXT_GET                   : Gets latitude slice of soil texture
!  (3 ) SFC_TYP_GET                   : Gets latitude slice of surface type
!  (4 ) TPT_GND_SOI_GET               : Gets latitude slice of soil & gnd tmp
!  (11) GET_ORO                       : Gets 2-D orography array
!  (23) SNW_FRC_GET                   : Converts LW snow depth to snow cover
!  (39) DST_TVBDS_GET                 : Gets a latitude slice of VAI data
!  (44) GET_TIME_INVARIANT_DATA       : Reads time-invariant fields from disk
!  (45) GET_MONTHLY_DATA              : Reads monthly fields from disk
!  (46) INIT_DUST_SHAO                : Allocates & zeroes module arrays
!  (47) CLEANUP_DUST_SHAO             : Deallocates
!
!  GEOS-CHEM modules referenced by dust_shao_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f       : Module containing arrays for GMAO met fields
!  (3 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (4 ) error_mod.f     : Module containing I/O error and NaN check routines
!  (5 ) grid_mod.f      : Module containing horizontal grid information
!  (6 ) time_mod.f      : Module containing routines for computing time & date
!  (7 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Added parallel DO loop in GET_ORO (bmy, 4/14/04)
!  (2 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Fixed typo in ORO_IS_LND for PGI compiler (bmy, 3/1/05)
!  (4 ) Modified for GEOS-5 and GCAP met fields (swu, bmy, 8/16/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now uses GOCART source function (tdf, bmy, 1/25/07)
!  (7 ) Modifications for 0.5 x 0.667 grid (yxw, dan, bmy, 11/6/08)
!  (8 ) Updates for nested grids (amv, bmy, 12/18/09)
!******************************************************************************
!
!****************************byk(090519)for shao scheme - start  ****************
! Y. Shao, 29 Jan 2004 
!-----------------------------------------------------------------------------------
! Calculate sediment flux for multi-particle size soils as a weighted average of Q(d) 
!           dust emission F(d) for covered and moisture soil
!
! Options for dust calculation:
! 1 Shao et al. (1996)
! 2 Marticorena and Bergametti (1995)
! 3 Lu and Shao (1999)
! 4 Shao (2000) + Lu and Shao (1999)
! 5 Shao (2001, 2003)
!
!--------------------------------------------------------------------------------------
!
!      input: 
!      n:          number of particle size ranges.                 ---> N
!      dm:         median diameter of each particle size.  [m]     
!      m_fract:    Weight fraction of each particle range. Sum m_fract = 1
!      u_s:        Friciton velocity.        [m/s]                 ---> FRC_VEL
!      cf:         fraction area covered by roughness elements     ---> VAI_DST_SLICE
!      w:          surface soil moisture contains [m^3/m^3]        ---> SOI_MST
!      c:          Owen's coefficient         
!
!      output:  
!      ustart:     Mean threshold velocity of each particle range.
!      q:          Sand flux from each size range.
!      ffq:        Weighted sand flux from each size range.
!      qtotal:     Total sand flux (weighted average).
!      f:          Weighted dust flux ejected by single sand size range for given d_d. 
!      fff:        Weighted dust flux ejected by all sand size range from each dust size range.
!      ftotal:     Total dust flux (weighted average).
!--------------------------------------------------------------------------------------
!****************************byk(090519)for shao scheme - end  *******************
!
      USE TRANS_REGRID_MOD
      USE CMN_SIZE_MOD, ONLY : MASTER

      IMPLICIT NONE
#     include "define.h"

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "dust_shao_mod.f"
      !=================================================================

      ! Make everything PRIVATE....
      PRIVATE

      ! Except these routines
      PUBLIC :: QF_SHAO
      PUBLIC :: CLEANUP_DUST_SHAO
      PUBLIC :: GET_ORO
      PUBLIC :: GET_TIME_INVARIANT_DATA
      PUBLIC :: GET_MONTHLY_DATA

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Fundamental physical constants
      REAL*8,  PARAMETER   :: GAS_CST_UNV      = 8.31441d0
      REAL*8,  PARAMETER   :: MMW_H2O          = 1.8015259d-02
      REAL*8,  PARAMETER   :: MMW_DRY_AIR      = 28.9644d-3
      REAL*8,  PARAMETER   :: CST_VON_KRM      = 0.4d0
      REAL*8,  PARAMETER   :: GRV_SFC          = 9.80616d0
      REAL*8,  PARAMETER   :: GAS_CST_DRY_AIR  = 287.05d0
      REAL*8,  PARAMETER   :: RDS_EARTH        = 6.37122d+6
      REAL*8,  PARAMETER   :: GAS_CST_H2O      = 461.65D0
      REAL*8,  PARAMETER   :: SPC_HEAT_DRY_AIR = 1005.0d0
      REAL*8,  PARAMETER   :: TPT_FRZ_PNT      = 273.15d0

      ! Derived quantities
      REAL*8,  PARAMETER   :: GRV_SFC_RCP      = 1.0d0   / GRV_SFC
      REAL*8,  PARAMETER   :: CST_VON_KRM_RCP  = 1.0d0   / CST_VON_KRM
      REAL*8,  PARAMETER   :: EPS_H2O          = MMW_H2O / MMW_DRY_AIR
      REAL*8,  PARAMETER   :: EPS_H2O_RCP_M1   = -1.0d0  + MMW_DRY_AIR
     &                                                   / MMW_H2O
      REAL*8,  PARAMETER   :: KAPPA_DRY_AIR    = GAS_CST_DRY_AIR
     &                                         / SPC_HEAT_DRY_AIR

      ! Land surface that is not lake or wetland (by area)
      REAL*8,  ALLOCATABLE :: LND_FRC_DRY(:,:)
      REAL*8,  ALLOCATABLE :: MSS_FRC_CLY(:,:)
      REAL*8,  ALLOCATABLE :: MSS_FRC_SND(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S01(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S02(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S03(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S04(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S05(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S06(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S07(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S08(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S09(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S10(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S11(:,:)
      REAL*8,  ALLOCATABLE :: FRC_S12(:,:)
      INTEGER, ALLOCATABLE :: SFC_TYP(:,:)

      ! Time-varying surface info from CTM
      REAL*8,  ALLOCATABLE :: TPT_GND(:,:)
      REAL*8,  ALLOCATABLE :: TPT_SOI(:,:)

      ! Variables initialized in dst_tvbds_ntp() and dst_tvbds_ini()
      REAL*8,  ALLOCATABLE :: VAI_DST(:,:)

      ! Other fields
      INTEGER, ALLOCATABLE :: OROGRAPHY(:,:)
      REAL*8,  ALLOCATABLE :: DMT_MIN(:)
      REAL*8,  ALLOCATABLE :: DMT_MAX(:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE QF_SHAO( DOY,         HGT_MDP,     LAT_IDX,
     &                    LAT_RDN,     ORO,         PRS_DLT,
     &                    PRS_MDP,     Q_H2O_VPR,   DSRC,
     &                    SNW_HGT_LQD, TM_ADJ,      TPT_MDP,
     &                    TPT_PTN_MDP, WND_MRD_MDP, WND_ZNL_MDP,
     &                    FIRST,       NSTEP,       FRC_VEL,
     &                    SOI_MST,     AREA_M2 )
!
!******************************************************************************
!  Subroutine QF_SHAO is the driver for aerosol mobilization (SHAO model).
!  It is designed to require only single layer surface fields, allowing for
!  easier implementation.  QF_SHAO is called once per latitude.  Modified
!  for GEOS-CHEM by Duncan Fairlie and Bob Yantosca.
!  (tdf, bmy, 1/25/07, 12/18/09)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DOY         (REAL*8 ) : Day of year [1.0..366.0)            [unitless]
!  (2 ) HGT_MDP     (REAL*8 ) : Midpoint height above surface       [m       ]
!  (3 ) LAT_IDX     (INTEGER) : Model latitude index                [unitless]
!  (4 ) LAT_RDN     (REAL*8 ) : Model latitude                      [radians ]
!  (5 ) ORO         (REAL*8 ) : Orography                           [fraction]
!  (6 ) PRS_DLT     (REAL*8 ) : Pressure thickness of grid box      [Pa      ]
!  (7 ) PRS_MDP     (REAL*8 ) : Pressure @ midpoint of grid box     [Pa      ]
!  (8 ) Q_H2O_VPR,  (REAL*8 ) : Water vapor mixing ratio            [kg/kg   ]
!  (9 ) SNW_HGT_LQD (REAL*8 ) : Equivalent liquid water snow depth  [m       ]
!  (10) TM_ADJ,     (REAL*8 ) : Adjustment timestep                 [s       ]
!  (11) TPT_MDP,    (REAL*8 ) : Temperature                         [K       ]
!  (12) TPT_PTN_MDP (REAL*8 ) : Midlayer local potential temp.      [K       ]
!  (13) WND_MRD_MDP (REAL*8 ) : Meridional wind component (V-wind)  [m/s     ]
!  (14) WND_ZNL_MDP (REAL*8 ) : Zonal wind component (U-wind)       [m/s     ]
!  (15) FIRST,      (LOGICAL) : Logical used ot open output dataset [unitless]
!  (16) NSTEP       (INTEGER) : Iteration counter                   [unitless]
!  (17) FRC_VEL     (REAL*8)  : Friction velocity                   [m/s     ]
!  (18) SOI_MST     (REAL*8)  : Surface soil moisture contains      [m3/m3   ]
!  (19) AREA_M2     (REAL*8)  : Grid box surface area               [m2      ]
!
!  Arguments as Output:
!  ============================================================================
!  (10) DSRC                ! O [kg kg-1] Dust mixing ratio increment
!
!  NOTES:
!  (1 ) Cleaned up and added comments.  Also force double precision with
!        "D" exponents. (bmy, 3/30/04)
!  (2 ) Now get GOCART source function. (tdf, bmy, 1/25/07)      
!  (3 ) Tune nested-domain emissions dust to the same as 2x2.5 simulation
!        Also tune GEOS-3 1x1 N. America nested-grid dust emissions to
!        the 4x5 totals from the GEOS-5 4x5 v8-01-01-Run0 benchmark. 
!        (yxw, bmy, dan, 11/6/08)
!  (4 ) New scale parameter for 2x2.5 GEOS-5 (tdf, jaf, phs, 10/30/09)
!  (5 ) Defined FLX_MSS_FDG_FCT for GEOS_4 2x2.5, GEOS_5 2x2.5, NESTED_NA and 
!        NESTED_EU.  Redefined FLX_MSS_FDG_FCT for NESTED_CH, based upon above
!        changes. (amv, bmy, 12/18/09)
!  (6 ) For now treat MERRA like GEOS-5 (bmy, 8/13/10)
!  29 Oct 2010 - T. D. Fairlie, R. Yantosca - Retune dust for MERRA 4x5
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,   ONLY : USTAR, Z0
      USE GRID_MOD,  ONLY : GET_AREA_M2
      USE ERROR_MOD, ONLY : ERROR_STOP

      USE CMN_SIZE_MOD      ! Size parameters
      USE CMN_GCTM_MOD      ! PI, g0

      ! Arguments
      INTEGER, INTENT(IN)    :: LAT_IDX
      REAL*8,  INTENT(IN)    :: DOY
      REAL*8,  INTENT(IN)    :: HGT_MDP(IIPAR)
      REAL*8,  INTENT(IN)    :: LAT_RDN
      REAL*8,  INTENT(IN)    :: ORO(IIPAR)
      REAL*8,  INTENT(IN)    :: PRS_DLT(IIPAR)
      REAL*8,  INTENT(IN)    :: PRS_MDP(IIPAR)
      REAL*8,  INTENT(IN)    :: Q_H2O_VPR(IIPAR)
      REAL*8,  INTENT(IN)    :: SNW_HGT_LQD(IIPAR)
      REAL*8,  INTENT(IN)    :: TM_ADJ
      REAL*8,  INTENT(IN)    :: TPT_MDP(IIPAR)
      REAL*8,  INTENT(IN)    :: TPT_PTN_MDP(IIPAR)
      REAL*8,  INTENT(IN)    :: WND_MRD_MDP(IIPAR)
      REAL*8,  INTENT(IN)    :: WND_ZNL_MDP(IIPAR)
      INTEGER, INTENT(IN)    :: NSTEP
      LOGICAL, INTENT(IN)    :: FIRST
      REAL*8,  INTENT(IN)    :: FRC_VEL(IIPAR)
      REAL*8,  INTENT(IN)    :: SOI_MST(IIPAR)
      REAL*8,  INTENT(IN)    :: AREA_M2
      REAL*8,  INTENT(INOUT) :: DSRC(IIPAR,NDSTBIN)

      !--------------
      ! Parameters
      !--------------
      ! [kg m-3] Density of optimal saltation particles,
      ! MBA97 p. 4388 friction velocity for saltation
      REAL*8,  PARAMETER   :: DNS_SLT     = 2650.0d0

      !--------------
      ! Local Output
      !--------------
      REAL*8 LND_FRC_MBL_SLICE(IIPAR)   ! [frc] Bare ground fraction
      REAL*8 SNW_FRC(IIPAR)             ! [frc] Fraction of surface covered
                                        !       by snow
      REAL*8 WND_MDP(IIPAR)             ! [m/s] Surface layer mean wind speed

      LOGICAL FLG_MBL_SLICE(IIPAR) ! [flg] Mobilization candidates
      INTEGER I                    ! [idx] Counting index
      INTEGER IJLOOP               ! [idx] counting index
      INTEGER M                    ! [idx] Counting index
      INTEGER MBL_NBR              ! [nbr] Number of mobilization candidates
      INTEGER SFC_TYP_SLICE(IIPAR) ! [idx] LSM surface type lat slice (0..28)
      REAL*8 DNS_MDP(IIPAR)             ! [kg/m3] Midlayer density
      REAL*8 RAT_DNS(IIPAR)             ! [unitless] particle-air density ratio (DNS_SLT/DNS_MDP)

      REAL*8 LND_FRC_DRY_SLICE(IIPAR)   ! [frc] Dry land fraction
      REAL*8 MSS_FRC_CLY_SLICE(IIPAR)   ! [frc] Mass fraction of clay
      REAL*8 MSS_FRC_SND_SLICE(IIPAR)   ! [frc] Mass fraction of sand
      REAL*8 MSS_FRC_LOA_SLICE(IIPAR)   ! [frc] Mass fraction of loam 
      REAL*8 MSS_FRC_SLO_SLICE(IIPAR)   ! [frc] Mass fraction of sandy loam

      REAL*8 MPL_AIR(IIPAR)             ! [kg/m2] Air mass path in layer

      REAL*8 TM_DLT                ! [s] Mobilization timestep
      REAL*8 TPT_GND_SLICE(IIPAR)  ! [K] Ground temperature
      REAL*8 TPT_SOI_SLICE(IIPAR)  ! [K] Soil temperature
      REAL*8 TPT_SOI_FRZ           ! [K] Temperature of frozen soil
      REAL*8 TPT_VRT_MDP           ! [K] Midlayer virtual temperature
      REAL*8 VAI_DST_SLICE(IIPAR)  ! [m2/m2] Vegetation area index,
                                   !         one-sided

      REAL*8 FRC_THR_NCR_WTR(IIPAR)     ! [frc] Threshold friction velocity
                                        !       increase from moisture
      REAL*8 VWC_DRY(IIPAR)        ! [m3/s] Dry volumetric water content
                                   !        (no E-T)
      REAL*8 VWC_OPT(IIPAR)        ! [m3/m3] E-T optimal volumetric water
                                   !         content
      REAL*8 VWC_SAT(IIPAR)        ! [m3/m3] Saturated volumetric water
                                   !         content (sand-dependent)
      REAL*8 GWC_SFC(IIPAR)        ! [kg/kg] Gravimetric water content


!****************************byk(090519)for shao scheme - start ******************

      ! constant.inc
        ! PI, g0                                         : include "CMN_GCTM" in QF_SHAO
        ! rhoa (air density)                 [kg/m3]     : DNS_MDP in QF_SHAO
        ! rhop (particle density)            [kg/m3]     : DNS_SLT in QF_SHAO
        ! sigma (particle-air density ratio) [unitless]  : sigma in QF_SHAO

      ! param_in.rdr

       REAL*8,  PARAMETER   :: PHL       = 30000.0d0     !plastic pressure [N/m2]  
       REAL*8,  PARAMETER   :: RHOS      = 1000.0d0      !bulk density of soil [kg/m3]  
       REAL*8,  PARAMETER   :: CYS       = 1.0d-5        !a parameter  
       REAL*8,  PARAMETER   :: CALPHA = 5.d0, CBETA = 1.37d0

!       INTEGER :: I
       INTEGER :: J, N, S, II
       INTEGER :: IDST, IDST1, IDST2, IDST3, IDST4, IDST5
!
       INTEGER, PARAMETER :: IMAX = 100                ! No. of particle size intervals for psd
       INTEGER, PARAMETER :: JMAX = 4                  ! No. of log-normal distributions for constructing psd
       INTEGER, PARAMETER :: IMOD = 5                  ! Options: imod=1,2,3,4,5, 5 is recommended
       INTEGER, PARAMETER :: S_TYPE = 4                ! soil type (1:sand, 2:sloa, 3:loam, 4:clay)
!
       REAL*8, dimension(0:imax) :: D0
       REAL*8, dimension(imax)   :: DD, PSDM, DPSDM, PPSDM
       REAL*8, dimension(imax)   ::     PSDF, DPSDF, PPSDF
       REAL*8, dimension(imax)   ::     PSDS, DPSDS, PPSDS
!
       REAL*8, dimension(S_TYPE)   :: FRC_SOI      ! soil type (s_type=1:sand, 2:sloa, 3:loam, 4:clay)
!
       REAL*8, dimension(imax)         :: USTART, Q, FFQ, FF, FFF
       REAL*8, dimension(imax,imax)    :: F
       REAL*8, dimension(S_TYPE)       :: QTOTAL, FTOTAL  ! horizontal, vertical dust flux each soil
       REAL*8, dimension(S_TYPE,NDSTBIN) :: F_BIN_SLICE   ! vertical dust flux each soil each bin
       REAL*8, dimension(IIPAR,NDSTBIN)  :: F_BIN         ! vertical dust flux each bin each grid
!
       REAL*8, parameter :: C_LAMBDA=0.35d0
       REAL*8 :: R, H, LAMBDA, UST_MIN
       REAL*8 :: CF, SM, U_S, SIGMA, RHOA
       REAL*8 :: GHL, FC
       REAL*8, parameter :: DCUT0=0.2d0              ! dust1 startup particle size
       REAL*8, parameter :: DCUT1=2.1d0               ! dust1 cutoff particle size
       REAL*8, parameter :: DCUT2=3.6d0               ! dust2 cutoff particle size
       REAL*8, parameter :: DCUT3=6.0d0               ! dust3 cutoff particle size
       REAL*8, parameter :: DCUT4=1.2d1              ! dust4 cutoff particle size
       REAL*8 :: RYS
       REAL*8 :: A1, A2, A3
       REAL*8 :: ZETA, SIGMA_M                      ! u*sqrt(rhos/p), bombardment coefficient
!
!       REAL*8 :: R_C, H_C
!       REAL*8 :: QWHITE, USTART0
!       external R_C, H_C, QWHITE, USTART0

!**********FROM SOILPSD.INC*******  start**************************************************
       INTEGER, parameter :: MMAX=4

       REAL*8, dimension(3, MMAX) :: CSANDM         ! Coefs for sand minimally dispersed
       data CSANDM /0.,     0.,     0.,
     &             0.0329, 4.3733, 0.8590,
     &             0.9671, 5.7689, 0.2526,
     &             0.,     0.,     0.      /
!
       REAL*8, dimension(3, MMAX) :: CLOAMM         ! Coefs for loam minimally dispersed
!       data cloamm /0.1114, 4.3565, 0.0257, & 
       data CLOAMM /0.1114, 4.3565, 0.4257,
     &             0.4554, 5.1674, 0.3824,
     &             0.4331, 5.4092, 1.0000,
     &             0.,     0.,     0.      /
!
       REAL*8, dimension(3, MMAX) :: CSLOAMM        ! Coefs for sandy clay loam minimally dispersed, very dusty
       data CSLOAMM/0.0722, 2.2675, 1.0000,
     &             0.6266, 4.9654, 0.3496,
     &             0.3012, 5.5819, 0.5893,
     &             0.,     0.,     0.      /
!
       REAL*8, dimension(3, MMAX) :: CCLAYM         ! Coefs for clay minimally dispersed
       data CCLAYM /0.3902, 3.5542, 1.0000,
     &             0.2813, 4.2239, 0.2507,
     &             0.3286, 5.1638, 0.4632,
     &             0.,     0.,     0.      /
!
       REAL*8, dimension(3, MMAX) :: CSANDF         ! Coefs for sand fully dispersed
       data CSANDF /0.,     0.,     0.,
     &             0.0338, 0.6931, 1.0000,
     &             0.9662, 5.6300, 0.2542,
     &             0.,     0.,     0.      /
!
       REAL*8, dimension(3, MMAX) :: CLOAMF         ! Coefs for loam fully dispersed
       data CLOAMF /0.5844, 4.6079, 0.6141,
     &             0.3304, 5.2050, 0.2897,
     &             0.0522, 7.0553, 1.0000,
     &             0.0330, 0.6931, 1.0000  /
!
       REAL*8, dimension(3, MMAX) :: CSLOAMF         ! Coefs for sandy clay loam fully dispersed
       data CSLOAMF/0.2344, 1.8079, 0.6141,
     &             0.3634, 4.2050, 0.2897,
     &             0.4022, 5.6553, 1.0000,
     &             0.,     0.,     0.      /
!
       REAL*8, dimension(3, MMAX) :: CCLAYF         ! Coefs for clay fully dispersed
       data CCLAYF /0.0872, 0.6931, 1.0000,
     &             0.4464, 3.9323, 0.9181,
     &             0.4665, 5.4486, 0.3916,
     &             0.,     0.,     0.      /
!
!**********FROM SOILPSD.INC*******  end ******************************************
!
!       WRITE(6,*) '***************************************************'
!       WRITE(6,*) CSANDM(1,2),CSANDM(2,2),CSANDM(3,2)
!       WRITE(6,*) '***************************************************'
!          pause

!****************************byk(090519)for shao scheme - end   ******************

      !=================================================================
      ! QF_SHAO begins here!
      !=================================================================

      ! Time step [s]
      TM_DLT                 = TM_ADJ

      ! Freezing pt of soil [K] -- assume it's 0C
      TPT_SOI_FRZ            = TPT_FRZ_PNT

      ! Initialize output fluxes and tendencies
      FRC_THR_NCR_WTR(:)     = 0.0D0       ! [frc]

      F_BIN(:,:)             = 0.0D0
      DSRC(:,:)              = 0.0D0

      !=================================================================
      ! Compute necessary derived fields
      !=================================================================
      DO I = 1, LONLOOP

         ! Stop occasional haywire model runs
         IF ( TPT_MDP(I) > 350.0d0 ) THEN
            CALL ERROR_STOP( 'TPT_MDP(i) > 350.0',
     &                       'QF_SHAO ("dust_shao_mod.f")' )
         ENDIF

         ! Midlayer virtual temperature [K]
         TPT_VRT_MDP = TPT_MDP(I)
     &               * (1.0d0 + EPS_H2O_RCP_M1 * Q_H2O_VPR(I))

         ! Density at center of gridbox [kg/m3]
         DNS_MDP(I) = PRS_MDP(I)
     &              / (TPT_VRT_MDP * GAS_CST_DRY_AIR)

!****************************byk(090519)for shao scheme - start ******************

         ! particle-air density ratio [unitless]          
         RAT_DNS(I)   = DNS_SLT / DNS_MDP(I)

!****************************byk(090519)for shao scheme - end   ******************

         ! Commented out
         !cApproximate surface virtual temperature (uses midlayer moisture)
         !c tpt_vrt_sfc=tpt_sfc(i)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i)) ! [K]
         !c
         !c Surface density
         !c dns_sfc(i)=prs_sfc(i)/(tpt_vrt_sfc*gas_cst_dry_air) ! [kg m-3]

         ! Mass of air currently in gridbox [kg/m2]
         MPL_AIR(I) = PRS_DLT(I) * GRV_SFC_RCP

         ! Mean surface layer horizontal wind speed
         WND_MDP(I) = SQRT( WND_ZNL_MDP(I)*WND_ZNL_MDP(I)
     &              +       WND_MRD_MDP(I)*WND_MRD_MDP(I) )

      ENDDO

      !=================================================================
      ! Gather input variables from GEOS-CHEM modules etc.
      !=================================================================

      ! Get LSM Surface type (0..28)
      CALL SFC_TYP_GET( LAT_IDX, SFC_TYP_SLICE )

      ! Get erodability and mass fractions
      CALL SOI_TXT_GET(
     &    LAT_IDX,             ! I [idx] Latitude index
     &    LND_FRC_DRY_SLICE,   ! O [frc] Dry land fraction
     &    MSS_FRC_CLY_SLICE,   ! O [frc] Mass fraction of clay
     &    MSS_FRC_SND_SLICE,   ! O [frc] Mass fraction of sand
     &    MSS_FRC_LOA_SLICE,   ! O [frc] Mass fraction of loam 
     &    MSS_FRC_SLO_SLICE )  ! O [frc] Mass fraction of sandy loam

      ! Get surface and soil temperature
      CALL TPT_GND_SOI_GET(
     &     LAT_IDX,            ! I [idx] Latitude index!
     &     TPT_GND_SLICE,      ! O [K] Ground temperature
     &     TPT_SOI_SLICE )     ! O [K] Soil temperature

      ! Get time-varying vegetation area index
      CALL DST_TVBDS_GET(
     &    LAT_IDX,             ! I [idx] Latitude index
     &    VAI_DST_SLICE)       ! O [m2 m-2] Vegetation area index, one-sided

      ! Get fraction of surface covered by snow
      CALL SNW_FRC_GET(
     &    SNW_HGT_LQD,         ! I [m] Equivalent liquid water snow depth
     &    SNW_FRC )            ! O [frc] Fraction of surface covered by snow

      !=================================================================
      ! Use the variables retrieved above to compute the fraction
      ! of each gridcell suitable for dust mobilization
      !=================================================================
      CALL LND_FRC_MBL_GET(
     &    DOY,                 ! I [day] Day of year [1.0..366.0)
     &    FLG_MBL_SLICE,       ! O [flg] Mobilization candidate flag
     &    LAT_RDN,             ! I [rdn] Latitude
     &    LND_FRC_DRY_SLICE,   ! I [frc] Dry land fraction
     &    LND_FRC_MBL_SLICE,   ! O [frc] Bare ground fraction
     &    MBL_NBR,             ! O [flg] Number of mobilization candidates
     &    ORO,                 ! I [frc] Orography
     &    SFC_TYP_SLICE,       ! I [idx] LSM surface type (0..28)
     &    SNW_FRC,             ! I [frc] Fraction of surface covered by snow
     &    TPT_SOI_SLICE,       ! I [K] Soil temperature
     &    TPT_SOI_FRZ,         ! I [K] Temperature of frozen soil
     &    VAI_DST_SLICE)       ! I [m2 m-2] Vegetation area index, one-sided

      ! Much ado about nothing
      if (mbl_nbr == 0) then
ctdf        print *,' no mobilisation candidates'
        goto 737
      endif

      !=================================================================
      ! Compute time-invariant hydrologic properties
      ! NB flg_mbl IS time-dependent, so keep this in time loop.
      !=================================================================
      CALL HYD_PRP_GET(        ! NB: These properties are time-invariant
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    MSS_FRC_CLY_SLICE,   ! I [frc] Mass fraction clay
     &    MSS_FRC_SND_SLICE,   ! I [frc] Mass fraction sand
     &    VWC_DRY,             ! O [m3/m3] Dry vol'mtric water content (no E-T)
     &    VWC_OPT,             ! O [m3/m3] E-T optimal volumetric water content
     &    VWC_SAT)             ! O [m3/m3] Saturated volumetric water content

      !=================================================================
      ! Convert volumetric water content to gravimetric water content
      ! NB: Owen effect included in wnd_frc_slt_get
      !=================================================================
      CALL VWC2GWC(
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    GWC_SFC,             ! O [kg kg-1] Gravimetric water content
     &    VWC_SAT,             ! I [m3 m-3] Saturated VWC (sand-dependent)
     &    SOI_MST )            ! I [m3 m-3] Volumetric water content

      !=================================================================
      ! Factor by which soil moisture increases threshold friction
      ! velocity -- i.e. the inhibition of saltation by soil mositure,
      ! Zender et al., exp(5).
      !=================================================================
      CALL FRC_THR_NCR_WTR_GET(
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    FRC_THR_NCR_WTR,     ! O [frc] Factor by which moisture increases
                               !         threshold friction velocity
     &    MSS_FRC_CLY_SLICE,   ! I [frc] Mass fraction of clay
     &    GWC_SFC)             ! I [kg kg-1] Gravimetric water content

!****************************byk(090519)for shao scheme - start ******************

      DO II = 1, LONLOOP        ! do loop for longitude

         CF    = dmin1 (VAI_DST_SLICE(II), 1.0d0)
         SM    = SOI_MST(II)
         U_S = FRC_VEL(II)
         RHOA  = DNS_MDP(II)
         SIGMA = RAT_DNS(II)

!-------------------
! frontal area index
!-------------------
         LAMBDA = - C_LAMBDA*dlog( 1.d0 - CF )
         R = R_C(LAMBDA)
         H = FRC_THR_NCR_WTR(II)

      Do S  = 1, S_TYPE         ! do loop for soil type (1:sand, 2:sloa, 3:loam, 4:clay)

       IF     (S .eq. 1) then
        call psd_create(D0, DD, PSDM, DPSDM, PPSDM, IMAX, CSANDM, JMAX)      ! psd for minimally dispersed
        call psd_create(D0, DD, PSDF, DPSDF, PPSDF, IMAX, CSANDF, JMAX)      ! psd for fully dispersed
        FRC_SOI(S) = MSS_FRC_SND_SLICE(II)

!--------------
! CUT SIZE BIN
!--------------
      IDST  =0
      IDST1 =0
      IDST2 =0
      IDST3 =0
      IDST4 =0
      IDST5 =0
      do I = 1,IMAX
        if     ( DD(I) .le. DCUT1 ) then
              IDST1 = IDST1 +1
        elseif ( DD(I) .gt. DCUT1 .and. DD(I) .le. DCUT2 ) then
              IDST2 = IDST2 + 1
        elseif ( DD(I) .gt. DCUT2 .and. DD(I) .le. DCUT3 ) then
              IDST3 = IDST3 + 1
        elseif ( DD(I) .gt. DCUT3 .and. DD(I) .le. DCUT4 ) then
              IDST4 = IDST4 + 1
        else
              IDST5 = IDST5 + 1
        endif
      enddo
        IDST = IDST1 + IDST2 +IDST3 +IDST4

!----------------------------------------------
! for each particle size group, estimate ustart
!----------------------------------------------
!
         UST_MIN = 999.d0
         do I = 1,IMAX
            USTART(I) = USTART0( DD(I),SIGMA, RHOA )
            USTART(I) = R * H * USTART(I)
            UST_MIN = dmin1(UST_MIN, USTART(I))
            Q(I) = QWHITE( USTART(I), U_S, RHOA )
         enddo

         if ( U_S .le. UST_MIN ) then
            Q                = 0.d0
            FFQ              = 0.d0
            FFF              = 0.d0
            F_BIN_SLICE(S,:) = 0.d0
            QTOTAL(S)        = 0.d0
            FTOTAL(S)        = 0.d0

         else
            GHL = dexp( -(U_S - UST_MIN)**3.d0 )
            DPSDS = GHL*DPSDM + (1.-GHL)*DPSDF
            PSDS  = GHL*PSDM  + (1.-GHL)*PSDF
            PPSDS = GHL*PPSDM + (1.-GHL)*PPSDF

            FFQ = Q * DPSDS
            QTOTAL(S) = sum(FFQ)

!--------------------------------
! Shao (2004) dust emission model
!--------------------------------
        IF (IMOD .eq. 5) THEN
          ZETA = U_S*dsqrt( RHOS/PHL )
          SIGMA_M = 12.d0*ZETA*ZETA*(1.d0+14.d0*ZETA)

          do I = IDST+1, IMAX
             do J = 1, IDST
                RYS = PSDM(J)/PSDF(J)
                A1 = CYS * DPSDF(J) * ( (1.-GHL) + GHL*RYS )
                A2 = (1. + SIGMA_M)
                A3 = FFQ(I)*g0/U_S**2
                F(I,J) = A1*A2*A3
             enddo
          enddo

          FTOTAL(S) = 0.d0
          FFF       = 0.d0
          do J = 1, IDST
             FF(J) = 0.d0
             do I = IDST+1, IMAX
                FF(J) = FF(J) + F(I,J)
             enddo
             FFF(J)    = FF(J)
             FTOTAL(S) = FTOTAL(S) + FFF(J)
          enddo

          F_BIN_SLICE(S,:) = 0.d0
          DO M =1, NDSTBIN
            if (M .eq. 1) then
              do J =1, IDST1
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo

            elseif (M .eq. 2) then
              do J =IDST1+1, IDST1+IDST2
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo

            elseif (M .eq. 3) then
              do J =IDST1+IDST2+1, IDST1+IDST2+IDST3
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo

            elseif (M .eq. 4) then
              do J =IDST1+IDST2+IDST3+1, IDST1+IDST2+IDST3+IDST4
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            endif

          ENDDO
        ENDIF       ! for IMOD = 5
      endif         ! for u_s .le. ust_min

       elseif (S .eq. 2) then                                        ! sandy loam
        call psd_create(D0, DD, PSDM, DPSDM, PPSDM, IMAX, CSLOAMM, JMAX)      
        call psd_create(D0, DD, PSDF, DPSDF, PPSDF, IMAX, CSLOAMF, JMAX)      
        FRC_SOI(S) = MSS_FRC_SLO_SLICE(II)

!--------------
! CUT SIZE BIN
!--------------
      IDST  =0
      IDST1 =0
      IDST2 =0
      IDST3 =0
      IDST4 =0
      IDST5 =0
      do I = 1,IMAX
        if     ( DD(I) .le. DCUT1 ) then
              IDST1 = IDST1 +1
        elseif ( DD(I) .gt. DCUT1 .and. DD(I) .le. DCUT2 ) then
              IDST2 = IDST2 + 1
        elseif ( DD(I) .gt. DCUT2 .and. DD(I) .le. DCUT3 ) then
              IDST3 = IDST3 + 1
        elseif ( DD(I) .gt. DCUT3 .and. DD(I) .le. DCUT4 ) then
              IDST4 = IDST4 + 1
        else
              IDST5 = IDST5 + 1
        endif
      enddo
        IDST = IDST1 + IDST2 +IDST3 +IDST4

!----------------------------------------------
! for each particle size group, estimate ustart
!----------------------------------------------
!
         UST_MIN = 999.d0
         do I = 1,IMAX
            USTART(I) = USTART0( DD(I),SIGMA, RHOA )
            USTART(I) = R * H * USTART(I)
            UST_MIN = dmin1(UST_MIN, USTART(I))
            Q(I) = QWHITE( USTART(I), U_S, RHOA )
         enddo

         if ( U_S .le. UST_MIN ) then
            Q                = 0.d0
            FFQ              = 0.d0
            FFF              = 0.d0
            F_BIN_SLICE(S,:) = 0.d0
            QTOTAL(S)        = 0.d0
            FTOTAL(S)        = 0.d0

         else
            GHL = dexp( -(U_S - UST_MIN)**3.d0 )
            DPSDS = GHL*DPSDM + (1.-GHL)*DPSDF
            PSDS  = GHL*PSDM  + (1.-GHL)*PSDF
            PPSDS = GHL*PPSDM + (1.-GHL)*PPSDF

            FFQ = Q * DPSDS
            QTOTAL(S) = sum(FFQ)

!--------------------------------
! Shao (2004) dust emission model
!--------------------------------
        IF (IMOD .eq. 5) THEN
          ZETA = U_S*dsqrt( RHOS/PHL )
          SIGMA_M = 12.d0*ZETA*ZETA*(1.d0+14.d0*ZETA)

          do I = IDST+1, IMAX
             do J = 1, IDST
                RYS = PSDM(J)/PSDF(J)
                A1 = CYS * DPSDF(J) * ( (1.-GHL) + GHL*RYS )
                A2 = (1. + SIGMA_M)
                A3 = FFQ(I)*g0/U_S**2
                F(I,J) = A1*A2*A3
             enddo
          enddo

          FTOTAL(S) = 0.d0
          FFF       = 0.d0
          do J = 1, IDST
             FF(J) = 0.d0
             do I = IDST+1, IMAX
                FF(J) = FF(J) + F(I,J)
             enddo
             FFF(J)    = FF(J)
             FTOTAL(S) = FTOTAL(S) + FFF(J)
          enddo

          F_BIN_SLICE(S,:) = 0.d0
          DO M =1, NDSTBIN
            if (M .eq. 1) then
              do J =1, IDST1
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 2) then
              do J =IDST1+1, IDST1+IDST2
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 3) then
              do J =IDST1+IDST2+1, IDST1+IDST2+IDST3
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 4) then
              do J =IDST1+IDST2+IDST3+1, IDST1+IDST2+IDST3+IDST4
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            endif
          ENDDO
        ENDIF       ! for IMOD = 5
      endif         ! for u_s .le. ust_min

       elseif (S .eq. 3) then
        call psd_create(D0, DD, PSDM, DPSDM, PPSDM, IMAX, CLOAMM, JMAX)      
        call psd_create(D0, DD, PSDF, DPSDF, PPSDF, IMAX, CLOAMF, JMAX)      
        FRC_SOI(S) = MSS_FRC_LOA_SLICE(II)

!--------------
! CUT SIZE BIN
!--------------
      IDST  =0
      IDST1 =0
      IDST2 =0
      IDST3 =0
      IDST4 =0
      IDST5 =0
      do I = 1,IMAX
        if     ( DD(I) .le. DCUT1 ) then
              IDST1 = IDST1 +1
        elseif ( DD(I) .gt. DCUT1 .and. DD(I) .le. DCUT2 ) then
              IDST2 = IDST2 + 1
        elseif ( DD(I) .gt. DCUT2 .and. DD(I) .le. DCUT3 ) then
              IDST3 = IDST3 + 1
        elseif ( DD(I) .gt. DCUT3 .and. DD(I) .le. DCUT4 ) then
              IDST4 = IDST4 + 1
        else
              IDST5 = IDST5 + 1
        endif
      enddo
        IDST = IDST1 + IDST2 +IDST3 +IDST4

!----------------------------------------------
! for each particle size group, estimate ustart
!----------------------------------------------
!
         UST_MIN = 999.d0
         do I = 1,IMAX
            USTART(I) = USTART0( DD(I),SIGMA, RHOA )
            USTART(I) = R * H * USTART(I)
            UST_MIN = dmin1(UST_MIN, USTART(I))
            Q(I) = QWHITE( USTART(I), U_S, RHOA )
         enddo

         if ( U_S .le. UST_MIN ) then
            Q                = 0.d0
            FFQ              = 0.d0
            FFF              = 0.d0
            F_BIN_SLICE(S,:) = 0.d0
            QTOTAL(S)        = 0.d0
            FTOTAL(S)        = 0.d0

         else
            GHL = dexp( -(U_S - UST_MIN)**3.d0 )
            DPSDS = GHL*DPSDM + (1.-GHL)*DPSDF
            PSDS  = GHL*PSDM  + (1.-GHL)*PSDF
            PPSDS = GHL*PPSDM + (1.-GHL)*PPSDF

            FFQ = Q * DPSDS
            QTOTAL(S) = sum(FFQ)

!--------------------------------
! Shao (2004) dust emission model
!--------------------------------
        IF (IMOD .eq. 5) THEN
          ZETA = U_S*dsqrt( RHOS/PHL )
          SIGMA_M = 12.d0*ZETA*ZETA*(1.d0+14.d0*ZETA)

          do I = IDST+1, IMAX
             do J = 1, IDST
                RYS = PSDM(J)/PSDF(J)
                A1 = CYS * DPSDF(J) * ( (1.-GHL) + GHL*RYS )
                A2 = (1. + SIGMA_M)
                A3 = FFQ(I)*g0/U_S**2
                F(I,J) = A1*A2*A3
             enddo
          enddo

          FTOTAL(S) = 0.d0
          FFF       = 0.d0
          do J = 1, IDST
             FF(J) = 0.d0
             do I = IDST+1, IMAX
                FF(J) = FF(J) + F(I,J)
             enddo
             FFF(J)    = FF(J)
             FTOTAL(S) = FTOTAL(S) + FFF(J)
          enddo

          F_BIN_SLICE(S,:) = 0.d0
          DO M =1, NDSTBIN
            if (M .eq. 1) then
              do J =1, IDST1
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 2) then
              do J =IDST1+1, IDST1+IDST2
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 3) then
              do J =IDST1+IDST2+1, IDST1+IDST2+IDST3
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 4) then
              do J =IDST1+IDST2+IDST3+1, IDST1+IDST2+IDST3+IDST4
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            endif
          ENDDO
        ENDIF       ! for IMOD = 5
      endif         ! for u_s .le. ust_min

       elseif (S .eq. 4) then
        call psd_create(D0, DD, PSDM, DPSDM, PPSDM, IMAX, CCLAYM, JMAX)      
        call psd_create(D0, DD, PSDF, DPSDF, PPSDF, IMAX, CCLAYF, JMAX)      
        FRC_SOI(S) = MSS_FRC_CLY_SLICE(II)

!--------------
! CUT SIZE BIN
!--------------
      IDST  =0
      IDST1 =0
      IDST2 =0
      IDST3 =0
      IDST4 =0
      IDST5 =0
      do I = 1,IMAX
        if     ( DD(I) .le. DCUT1 ) then
              IDST1 = IDST1 +1
        elseif ( DD(I) .gt. DCUT1 .and. DD(I) .le. DCUT2 ) then
              IDST2 = IDST2 + 1
        elseif ( DD(I) .gt. DCUT2 .and. DD(I) .le. DCUT3 ) then
              IDST3 = IDST3 + 1
        elseif ( DD(I) .gt. DCUT3 .and. DD(I) .le. DCUT4 ) then
              IDST4 = IDST4 + 1
        else
              IDST5 = IDST5 + 1
        endif
      enddo
        IDST = IDST1 + IDST2 +IDST3 +IDST4

!----------------------------------------------
! for each particle size group, estimate ustart
!----------------------------------------------
!
         UST_MIN = 999.d0
         do I = 1,IMAX
            USTART(I) = USTART0( DD(I),SIGMA, RHOA )
            USTART(I) = R * H * USTART(I)
            UST_MIN = dmin1(UST_MIN, USTART(I))
            Q(I) = QWHITE( USTART(I), U_S, RHOA )
         enddo

         if ( U_S .le. UST_MIN ) then
            Q                = 0.d0
            FFQ              = 0.d0
            FFF              = 0.d0
            F_BIN_SLICE(S,:) = 0.d0
            QTOTAL(S)        = 0.d0
            FTOTAL(S)        = 0.d0

         else
            GHL = dexp( -(U_S - UST_MIN)**3.d0 )
            DPSDS = GHL*DPSDM + (1.-GHL)*DPSDF
            PSDS  = GHL*PSDM  + (1.-GHL)*PSDF
            PPSDS = GHL*PPSDM + (1.-GHL)*PPSDF

            FFQ = Q * DPSDS
            QTOTAL(S) = sum(FFQ)

!--------------------------------
! Shao (2004) dust emission model
!--------------------------------
        IF (IMOD .eq. 5) THEN
          ZETA = U_S*dsqrt( RHOS/PHL )
          SIGMA_M = 12.d0*ZETA*ZETA*(1.d0+14.d0*ZETA)

          do I = IDST+1, IMAX
             do J = 1, IDST
                RYS = PSDM(J)/PSDF(J)
                A1 = CYS * DPSDF(J) * ( (1.-GHL) + GHL*RYS )
                A2 = (1. + SIGMA_M)
                A3 = FFQ(I)*g0/U_S**2
                F(I,J) = A1*A2*A3
             enddo
          enddo

          FTOTAL(S) = 0.d0
          FFF       = 0.d0
          do J = 1, IDST
             FF(J) = 0.d0
             do I = IDST+1, IMAX
                FF(J) = FF(J) + F(I,J)
             enddo
             FFF(J)    = FF(J)
             FTOTAL(S) = FTOTAL(S) + FFF(J)
          enddo

          F_BIN_SLICE(S,:) = 0.d0
          DO M =1, NDSTBIN
            if (M .eq. 1) then
              do J =1, IDST1
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 2) then
              do J =IDST1+1, IDST1+IDST2
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 3) then
              do J =IDST1+IDST2+1, IDST1+IDST2+IDST3
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            elseif (M .eq. 4) then
              do J =IDST1+IDST2+IDST3+1, IDST1+IDST2+IDST3+IDST4
                F_BIN_SLICE(S,M) = F_BIN_SLICE(S,M) + FFF(J)
              enddo
            endif
          ENDDO
        ENDIF       ! for IMOD = 5
      endif         ! for u_s .le. ust_min

       ENDIF                     ! FOR SOIL TEXTURE (S=1,2,3,4)
      ENDDO                      ! for S (soil type)

         DO M= 1, NDSTBIN
         DO S= 1, S_TYPE
         F_BIN(II,M) = F_BIN(II,M) +            ! verticla dust flux each bin each grid [kg/m2/s]
     &         F_BIN_SLICE(S,M)*FRC_SOI(S)
         ENDDO
         ENDDO

      ENDDO                    ! for II (longitude)

!DO LOOP OUTPUT 
      DO  I=1, LONLOOP          ! NB: Inefficient loop order
         IF (FLG_MBL_SLICE(I)) THEN

            ! Loop over dust bins
            DO M = 1, NDSTBIN

               ! Introduce DSRC: dust mixing ratio increment   12/9/2K3
               DSRC(I,M) =      ! [kg]
     &              TM_ADJ * AREA_M2 * F_BIN(I,M)

           ENDDO
         ENDIF
      ENDDO

      ! Jump to here when no points are mobilization candidates
  737 CONTINUE

      ! Return to calling program
      END SUBROUTINE QF_SHAO

!****************************byk(090519)for shao scheme - start  ****************

!********************************************************************
      REAL*8 function r_c(x)
!
!   Y. Shao 17-07-92
!   CORRECTION FUNCTION FOR UST(D) BASED ON Raupach et al. (1992)
!   x = frontal area index
!
!   R_C = (1 - sig m x)^{1/2} (1 + m beta x)^{1/2}
!         Note I deife R_C = u_{*tR}/u_{*tS}
!         While Raupach et al. defined 
!                      R_C = u_{*tS}/u_{*tR} and their R function is
!   R_C = (1 - sig m x)^{-1/2} (1 + m beta x)^{-1/2}
!                      
!   sig   : basal to frontal area; about 1
!   m     : parameter less than 1; about 0.5
!   beta  : a ratio of drag coef.; about 90.
!
      REAL*8 :: X, XC
      REAL*8, parameter :: SIG=1., M=0.5, BETA=90.

      XC = 1./(SIG*M)
      IF (X .ge. XC) THEN
        R_C = 999.           ! Full covered surface
      ELSE
        R_C = dsqrt(1.-SIG*M*X)*dsqrt(1.+M*BETA*X)
      ENDIF
!
      END FUNCTION r_c

!******************************************************************
      real*8 function h_c(sm)
!
!   Y. Shao 17-06-2000
!   Correction of ust for soil water content
!   w: water content in m^3/m^3, NOT %
!   H_C is an empirical function based on the observations made by 
!   Shao during the 1993 wind tunnel experiments. See notes
!   H_C = exp(22.7 w) for w < 0.03
!   H_C = exp(95.3 w - 2.029) w > 0.03  ???
!
      REAL*8 :: SM
      IF (SM .lt. 0.) THEN
        WRITE(*,*) 'h_c, SM = ', SM, ' less than 0 at dust_shao_mod.f'
        STOP
      ENDIF

      H_C = dexp(22.7*SM)
      if (H_C .gt. 20.) H_C=20.     ! h_c=20 corresponding to w=0.132

      END FUNCTION h_c

!******************************************************************
      real*8 function ustart0(dum,sigma,rhoa)
!
! Y. Shao, 13 June 2000
!
! Calculate ustart0(d) using Shao and Lu (2000) for uncovered
! dry surface
!
! dum:    particle diameter                     [um]
! ustart0: threshold friction velocity           [m/s]
!
      USE CMN_GCTM_MOD        ! g0

      real(8) :: DM, DUM, SIGMA, RHOA
      real(8), parameter :: GAMMA = 1.65d-4      ! a constant
      real(8), parameter :: F = 0.0123d0

      DM = DUM*1d-6

      USTART0 = F*(SIGMA*g0*DM + GAMMA/(RHOA*DM) )
      USTART0 = dsqrt( USTART0 )

      END FUNCTION ustart0

!********************************************************************
      real(8) function qwhite(ust, u_s, rhoa)
!
! Yaping Shao 17-07-99!
!
! White (1979) Sand Flux Equation
! Q = c*rho*u_*^3 over g (1 - u_*t over u_*)(1 + u_*t^2/u_*^2)
! qwhite: Streamwise Sand Flux;       [kg m-1 s-1]
! c     : 2.6
! ust   : threhold friction velocity  [m/s]
! u_s : friction velocity           [m/s]
!C
      USE CMN_GCTM_MOD        ! g0

      real(8) :: C, UST, U_S, RHOA
      real(8) :: A, B

      C = 0.5d0
      A = RHOA/g0
      IF (U_S .lt. UST) THEN
        QWHITE = 0.d0
      ELSE
        B = UST/U_S
        QWHITE = C*A*U_S**3.*(1.-B)*(1.+B*B)
      ENDIF

      END FUNCTION qwhite

!****************************************************************************
      SUBROUTINE psd_create(D, DM, PSD, DPSD, PPSD, IMAX, CMTRIX, JMAX)
!
!----------------------------------------------------------------------------
! Yaping Shao, 13 June 2000
!     
! - Generate particle size distribution density function 
!   (both minimally-dispersed and fully-dispersed as the 
!   sum of four log-normal distributions.
!     
! d(0,imax):    output, particle size at 0, 1, 2, ..., imax points      [um]
! dm(imax):     output, particle size at middle of 0-1, 1-2, etc        [um]
! psd(imax):    output, particle size distribution density at dm        [um^-1]
! dpsd(imax):   output, Delta P for sections 0-1, 1-2, etc.             [ ]
! ppsd(imax):   output, P for sections 0-1, 1-2, etc.                   [ ]
! imax:         input, length dm, psd, dpsdm, ppsd, etc.
! cmtrix:       jmaxx coefficient matrix 
!               e.g.
!               w1   = cmtrix(1, 1): weight for first log-normal distribution
!               dln1 = cmtrix(2, 1): mean log-particle size of first log-normal distribution    
!               sig1 = cmtrix(3, 1): sigma of log-particle size for first log-normal distribution
!               etc. 
!               careful with the dimension of dln and sig
!----------------------------------------------------------------------------

      USE CMN_GCTM_MOD        ! PI, g0

      INTEGER :: I, J, IMAX, JMAX
      REAL*8, dimension(3, JMAX) :: CMTRIX
      REAL*8 :: D(0:IMAX), DM(IMAX)
      REAL*8 :: PSD(IMAX), DPSD(IMAX), PPSD(IMAX)   ! for p(d), Delta P(d) and P(d)
      REAL*8 :: P, PP, W, DLN, SIG
      REAL*8 :: CN
      REAL*8, parameter :: EPS   =   1.d-7
      REAL*8, parameter :: DREF  = 1000.d0
      REAL*8 :: FU, FD, PHI
!     
      CN = 1.d0/dsqrt(2.d0*PI)
!
! initialise psd, dpsd, ppsd
!     
      PSD = 0.d0
      DPSD = 0.d0
      PPSD = 0.d0
!
! Estimate d using phi scale. phi varies between from 9 to -1
! with increment 0.1. Reference particle size d0 = 1000 um
!
      FU =  9.d0
      FD = -1.d0
      do I = 0, IMAX
        PHI = FU - I*(FU-FD)/IMAX
        D(I) = DREF/2.d0**PHI
      enddo
!
      do I = 1, IMAX
        DM(I) = dexp( (dlog(D(I))+dlog(D(I-1)) )/2.d0 )
        PP = 0.d0
        do J = 1, JMAX
          W   = CMTRIX(1, J)
          DLN = CMTRIX(2, J)
          SIG = CMTRIX(3, J)
          if ( (W .gt. EPS) .and. (SIG .ne. 0.) ) then
            P = W*CN/SIG*dexp(-(dlog(DM(I))-DLN)**2/(2*SIG**2) )
          else
            P=0.d0
          endif
          PP = PP + P
        enddo
!     
        DPSD(I) = PP * ( dlog(D(I)) - dlog(D(I-1)) )        ! Delta P over i
        if (I .eq. 1) then
          PPSD(I) =      0.d0 + DPSD(I)                        ! P(d), with P(0) = 0
        else
          PPSD(I) = PPSD(I-1) + DPSD(I)
        endif
        PSD(I) = PP/DM(I)                                 ! p(d), particle size distribution density

      enddo
!
! Renormalisation, in case ppsd(imax) is not 1
!
      DPSD = DPSD/PPSD(IMAX)
       PSD =  PSD/PPSD(IMAX)
      PPSD = PPSD/PPSD(IMAX)
!
      END SUBROUTINE PSD_CREATE

!****************************byk(090519)for shao scheme - end  *******************

!------------------------------------------------------------------------------

      SUBROUTINE SOI_TXT_GET( J,               LND_FRC_DRY_OUT,
     &                        MSS_FRC_CLY_OUT, MSS_FRC_SND_OUT,
     &                        MSS_FRC_LOA_OUT, MSS_FRC_SLO_OUT )
!
!******************************************************************************
!  Subroutine SOI_GET_TXT returns a latitude slice of soil texture to the
!  calling program DST_MBL.  (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J                 (INTEGER) : Grid box latitude index
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) lnd_frc_dry_out   (REAL*8 ) : Dry land fraction      [fraction]
!  (5 ) mss_frc_cly_out   (REAL*8 ) : Mass fraction of clay  [fraction]
!  (6 ) mss_frc_snd_out   (REAL*8 ) : Mass fraction of sand  [fraction]
!  (7 ) mss_frc_loa_out   (REAL*8 ) : Mass fraction of loam  [fraction]
!  (8 ) mss_frc_slo_out   (REAL*8 ) : Mass fraction of sandy loam  [fraction]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!

      USE CMN_SIZE_MOD

      ! Arguments
      INTEGER, INTENT(IN)  :: J
      REAL*8,  INTENT(OUT) :: LND_FRC_DRY_OUT(IIPAR)
      REAL*8,  INTENT(OUT) :: MSS_FRC_CLY_OUT(IIPAR)
      REAL*8,  INTENT(OUT) :: MSS_FRC_SND_OUT(IIPAR)
      REAL*8,  INTENT(OUT) :: MSS_FRC_LOA_OUT(IIPAR)
      REAL*8,  INTENT(OUT) :: MSS_FRC_SLO_OUT(IIPAR)

      ! Local variables
      INTEGER              :: I

      ! Ad hoc globally uniform clay mass fraction [kg/kg]
      REAL*8,  PARAMETER   :: MSS_FRC_CLY_GLB = 0.20d0

      !=================================================================
      ! SOI_GET_TXT begins here!
      !=================================================================
      DO I = 1, LONLOOP

         ! Save dry land fraction slice
         LND_FRC_DRY_OUT(I) = LND_FRC_DRY(I,J)

         ! fxm Temporarily set mss_frc_cly used in mobilization to globally
         !     uniform SGS value of 0.20, and put excess mass fraction
         !     into sand
!         MSS_FRC_CLY_OUT(I) = MSS_FRC_CLY_GLB
!         MSS_FRC_SND_OUT(I) = MSS_FRC_SND(I,J) +
!     &                        MSS_FRC_CLY(I,J) - MSS_FRC_CLY_GLB
         MSS_FRC_SND_OUT(I) = FRC_S01(I,J) + FRC_S02(I,J)      ! sand       (stype=1)
         MSS_FRC_LOA_OUT(I) = FRC_S06(I,J) + FRC_S08(I,J) +    ! loam       (stype=3)  
     &                        FRC_S09(I,J) + FRC_S10(I,J) +
     &                        FRC_S04(I,J)
         MSS_FRC_SLO_OUT(I) = FRC_S03(I,J) + FRC_S07(I,J)      ! sandy loam (stype=2)   
         MSS_FRC_CLY_OUT(I) = FRC_S11(I,J) + FRC_S12(I,J)      ! clay       (stype=4)  

      ENDDO

      ! Return to calling program
      END SUBROUTINE SOI_TXT_GET

!------------------------------------------------------------------------------

      SUBROUTINE SFC_TYP_GET( J, SFC_TYP_OUT )
!
!******************************************************************************
!  Subroutine SFC_TYP_GET returns a latitude slice of LSM surface type
!  to the calling programs DST_MBL & DST_DPS_DRY. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J           (INTEGER) : Grid box latitude index
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) sfc_typ_out (REAL*8 ) : LSM surface type (0..28) [unitless]
!
!  NOTES
!  (1 ) Updated comments & cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!
      USE CMN_SIZE_MOD    ! IIPAR

      ! Arguments
      INTEGER, INTENT(IN)  :: J
      INTEGER, INTENT(OUT) :: SFC_TYP_OUT(IIPAR)

      ! Local variables
      INTEGER              :: I

      !=================================================================
      ! SFC_TYP_GET begins here!
      !=================================================================
      DO I = 1, LONLOOP
         SFC_TYP_OUT(I) = SFC_TYP(I,J)
      ENDDO

      ! Return to calling program
      END SUBROUTINE SFC_TYP_GET                       ! end sfc_typ_get()

!------------------------------------------------------------------------------

      SUBROUTINE TPT_GND_SOI_GET( J, TPT_GND_OUT, TPT_SOI_OUT )
!
!******************************************************************************
!  Subroutine TPT_GND_SOI_GET returns a latitude slice of soil temperature and
!  ground temperature to the calling program DST_MBL. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J           (INTEGER) : Grid box latitude index
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) TPT_GND_OUT (REAL*8 ) : Ground temperature array slice [K]
!  (3 ) tpt_soi_out (REAL*8 ) : Soil temperature array slice   [K]
!
!  NOTES
!  (1 ) Updated comments & cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : TS

      USE CMN_SIZE_MOD   ! IIPAR

      ! Arguments
      INTEGER, INTENT(IN)  :: J
      REAL*8,  INTENT(OUT) :: TPT_GND_OUT(IIPAR)
      REAL*8,  INTENT(OUT) :: TPT_SOI_OUT(IIPAR)

      ! Local variables
      INTEGER              :: I

      !=================================================================
      ! TPT_GND_SOI_GET begins here!
      !=================================================================

      ! Use TS from GEOS-CHEM (tdf, 3/30/04)
      DO I = 1, LONLOOP
         TPT_GND_OUT(I) = TS(I,J)
         TPT_SOI_OUT(I) = TS(I,J)
      ENDDO

      ! Return to calling program
      END SUBROUTINE TPT_GND_SOI_GET

!------------------------------------------------------------------------------

      SUBROUTINE GET_ORO( OROGRAPHY )
!
!******************************************************************************
!  Subroutine GET_ORO creates a 2D orography array, OROGRAPHY, from the
!  GMAO LWI fields.  Ocean= 0; Land=1; ice=2. (tdf, bmy, 3/30/04, 8/17/05)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) OROGRAPHY (INTEGER) : Array for orography flags
!
!  NOTES:
!  (1 ) Added parallel DO-loop (bmy, 4/14/04)
!  (2 ) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 6/9/05)
!  (3 ) Now use IS_LAND, IS_WATER, IS_ICE functions from "dao_mod.f"
!        (bmy, 8/17/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : IS_LAND, IS_WATER, IS_ICE

      USE CMN_SIZE_MOD

      ! Arguments
      INTEGER, INTENT(OUT)  :: OROGRAPHY(IIPAR,JJPAR)

      ! Local variables
      INTEGER :: I, J, TEMP

      !=================================================================
      ! GET_ORO begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, LATLOOP
      DO I = 1, LONLOOP

         ! Ocean
         IF ( IS_WATER( I, J ) ) OROGRAPHY(I,J) = 0

         ! Land
         IF ( IS_LAND(  I, J ) ) OROGRAPHY(I,J) = 1

         ! Ice
         IF ( IS_ICE (  I, J ) ) OROGRAPHY(I,J) = 2

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE GET_ORO

!------------------------------------------------------------------------------

      SUBROUTINE HYD_PRP_GET( FLG_MBL, MSS_FRC_CLY, MSS_FRC_SND,
     &                        VWC_DRY, VWC_OPT,     VWC_SAT )
!
!******************************************************************************
!  Subroutine HYD_PRP_GET determines hydrologic properties from soil texture.
!  (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FLG_MBL     (LOGICAL) : Mobilization candidate flag [unitless]
!  (2 ) MSS_FRC_CLY (REAL*8 ) : Mass fraction clay          [fraction]
!  (3 ) MSS_FRC_SND (REAL*8 ) : Mass fraction sand          [fraction]
!
!  Arguments as Output:
!  ============================================================================
!  (4 ) VWC_DRY     (REAL*8 ) : Dry volumetric water content (no E-T) [m3/m3]
!  (5 ) VWC_OPT     (REAL*8 ) : E-T optimal volumetric water content  [m3/m3]
!  (6 ) VWC_SAT     (REAL*8 ) : Saturated volumetric water content    [m3/m3]
!
!  NOTES:
!  (1 ) All I/O for this routine is time-invariant, thus, the hydrologic
!        properties could be computed once at initialization.  However,
!        FLG_MBL is time-dependent, so we should keep this as-is.
!        (tdf, 10/27/03)
!******************************************************************************
!
      USE CMN_SIZE_MOD   ! IIPAR

      ! Arguments
      LOGICAL, INTENT(IN)  :: FLG_MBL(IIPAR)
      REAL*8,  INTENT(IN)  :: MSS_FRC_CLY(IIPAR)
      REAL*8,  INTENT(IN)  :: MSS_FRC_SND(IIPAR)
      REAL*8,  INTENT(OUT) :: VWC_DRY(IIPAR)
      REAL*8,  INTENT(OUT) :: VWC_OPT(IIPAR)
      REAL*8,  INTENT(OUT) :: VWC_SAT(IIPAR)

      ! Local variables
      INTEGER              :: LON_IDX

      ! [frc] Exponent "b" for smp (clay-dependent)
      REAL*8               :: SMP_XPN_B(IIPAR)

      ! [mm H2O] Saturated soil matric potential (sand-dependent)
      REAL*8               :: SMP_SAT(IIPAR)

      !=================================================================
      ! HYD_PRP_GET begins here
      !=================================================================

      ! Initialize output values
      VWC_DRY(:) = 0.0D0
      VWC_OPT(:) = 0.0D0
      VWC_SAT(:) = 0.0D0

      ! Time-invariant soil hydraulic properties
      ! See Bon96 p. 98, implemented in CCM:lsm/lsmtci()
      DO LON_IDX = 1, LONLOOP

         IF ( FLG_MBL(LON_IDX) ) THEN

           ! Exponent "b" for smp (clay-dependent) [fraction]
           SMP_XPN_B(LON_IDX) =
!     &         2.91D0 +0.159D0 * MSS_FRC_CLY(LON_IDX) * 100.0D0
     &         2.91D0 +0.159D0 * 0.20d0 * 100.0D0

           ! NB: Adopt convention that matric potential is positive definite
           ! Saturated soil matric potential (sand-dependent) [mm H2O]
           SMP_SAT(LON_IDX) =
     &         10.0D0 * (10.0D0**(1.88D0-0.0131D0
     &                          * MSS_FRC_SND(LON_IDX)*100.0D0))

           ! Saturated volumetric water content (sand-dependent) ! [m3 m-3]
           VWC_SAT(LON_IDX)=
     &         0.489D0 - 0.00126D0 * MSS_FRC_SND(LON_IDX)*100.0D0

           ! [m3 m-3]
           VWC_DRY(LON_IDX) =

                ! Dry volumetric water content (no E-T)
     &          VWC_SAT(LON_IDX)*(316230.0D0/SMP_SAT(LON_IDX))
     &                       **(-1.0D0/SMP_XPN_B(LON_IDX))

           ! E-T optimal volumetric water content! [m3 m-3]
           VWC_OPT(LON_IDX) =
     &         VWC_SAT(LON_IDX)*(158490.0D0/SMP_SAT(LON_IDX))
     &                        **(-1.0D0/SMP_XPN_B(LON_IDX))
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE HYD_PRP_GET

!------------------------------------------------------------------------------

      LOGICAL FUNCTION ORO_IS_OCN( ORO_VAL )
!
!******************************************************************************
!  Function ORO_IS_OCN returns TRUE if a grid box contains more than 50%
!  ocean. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ORO_VAL (REAL*8) : Orography at a grid box (0=ocean; 1=land; 2=ice)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: ORO_VAL

      !=================================================================
      ! ORO_IS_OCN begins here!
      !=================================================================
      ORO_IS_OCN = ( NINT( ORO_VAL ) == 0 )

      ! Return to calling program
      END FUNCTION ORO_IS_OCN

!------------------------------------------------------------------------------

      LOGICAL FUNCTION ORO_IS_LND( ORO_VAL )
!
!******************************************************************************
!  Function ORO_IS_LND returns TRUE if a grid box contains more than 50%
!  land. (tdf, bmy, 3/30/04, 3/1/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ORO_VAL (REAL*8) : Orography at a grid box (0=ocean; 1=land; 2=ice)
!
!  NOTES:
!  (1 ) Bug fix: Replaced ": :" with "::" in order to prevent compile error
!        on Linux w/ PGI compiler.  (bmy, 3/1/05)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: ORO_VAL

      !=================================================================
      ! ORO_IS_OCN begins here!
      !=================================================================
      ORO_IS_LND = ( NINT( ORO_VAL ) == 1 )

      ! Return to calling program
      END FUNCTION ORO_IS_LND

!------------------------------------------------------------------------------

      LOGICAL FUNCTION ORO_IS_ICE( ORO_VAL )
!
!******************************************************************************
!  Function ORO_IS_LND returns TRUE if a grid box contains more than 50%
!  ice. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ORO_VAL (REAL*8) : Orography at a grid box (0=ocean; 1=land; 2=ice)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: ORO_VAL

      !=================================================================
      ! ORO_IS_ICE begins here!
      !=================================================================
      ORO_IS_ICE = ( NINT( ORO_VAL ) == 2 )

      ! Return to calling program
      END FUNCTION ORO_IS_ICE

!------------------------------------------------------------------------------

      SUBROUTINE SNW_FRC_GET( SNW_HGT_LQD, SNW_FRC )
!
!******************************************************************************
!  Subroutine SNW_FRC_GET converts equivalent liquid water snow depth to
!  fractional snow cover.  Uses the snow thickness -> fraction algorithm of
!  Bon96.  (tdf bmy, 3/30/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) snw_hgt_lqd (REAL*8) : Equivalent liquid water snow depth [m]
!
!  Arguments as Output:
!  ===========================================================================
!  (2 ) snw_frc     (REAL*8 ) : Fraction of surface covered by snow
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents. (bmy, 3/30/04)
!******************************************************************************
!
      USE CMN_SIZE_MOD

      !----------------
      ! Arguments
      !----------------
      REAL*8, INTENT(IN)  :: SNW_HGT_LQD(IIPAR)
      REAL*8, INTENT(OUT) :: SNW_FRC(IIPAR)

      !----------------
      ! Parameters
      !----------------

      ! Note disparity in bulk snow density between CCM and LSM
      ! WiW80 p. 2724, 2725 has some discussion of bulk snow density
      !
      ! Bulk density of snow [kg m-3]
      REAL*8,  PARAMETER  :: DNS_H2O_SNW_GND_LSM = 250.0D0

      ! Standard bulk density of snow on ground [kg m-3]
      REAL*8,  PARAMETER  :: DNS_H2O_SNW_GND_STD = 100.0D0

      ! Geometric snow thickness for 100% coverage ! [m]
      REAL*8,  PARAMETER  :: SNW_HGT_THR         = 0.05D0

      ! Liquid water density! [kg/m3]
      REAL*8,  PARAMETER  :: DNS_H2O_LQD_STD     = 1000.0D0

      !-----------------
      ! Local variables
      !-----------------

      ! [idx] Counting index for lon
      INTEGER             :: LON_IDX

      ! [m] Geometric bulk thickness of snow
      REAL*8              :: SNW_HGT(IIPAR)

      ! Conversion factor from liquid water depth
      ! to geometric snow thickness [fraction]
      REAL*8              :: HGT_LQD_SNW_CNV

      !=================================================================
      ! SNW_FRC_GET begins here!
      !=================================================================

      ! Conversion factor from liquid water depth to
      ! geometric snow thickness [fraction]
      HGT_LQD_SNW_CNV = DNS_H2O_LQD_STD
     &                / DNS_H2O_SNW_GND_STD

      ! Fractional snow cover
      DO LON_IDX = 1, LONLOOP

         ! Snow height [m]
         SNW_HGT(LON_IDX) = SNW_HGT_LQD(LON_IDX)
     &                    * HGT_LQD_SNW_CNV

         ! Snow fraction
         ! NB: CCM and LSM seem to disagree on this
         SNW_FRC(LON_IDX) = MIN(SNW_HGT(LON_IDX)/SNW_HGT_THR, 1.0D0)
      ENDDO

      ! Return to calling program
      END SUBROUTINE SNW_FRC_GET

!------------------------------------------------------------------------------

      SUBROUTINE VWC2GWC( FLG_MBL, GWC_SFC, VWC_SAT, SOI_MST )
!
!******************************************************************************
!  Subroutine VWC2GWC converts volumetric water content to gravimetric water
!  content -- assigned only for mobilisation candidates. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) FLG_MBL (LOGICAL) : Mobilization candidate flag     [flag]
!  (3 ) VWC_SAT (REAL*8 ) : Saturated VWC (sand-dependent)  [m3/m3]
!  (4 ) SOI_MST (REAL*8 ) : Volumetric water content!       [m3/m3
!
!  Arguments as Output:
!  ===========================================================================
!  (2 ) gwc_sfc (REAL*8 ) : Gravimetric water content       [kg/kg]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 3/30/04)
!******************************************************************************
!
      USE CMN_SIZE_MOD

      !----------------
      ! Arguments
      !----------------
      LOGICAL, INTENT(IN)  :: FLG_MBL(IIPAR)
      REAL*8,  INTENT(IN)  :: VWC_SAT(IIPAR)
      REAL*8,  INTENT(IN)  :: SOI_MST(IIPAR)
      REAL*8,  INTENT(OUT) :: GWC_SFC(IIPAR)

      !----------------
      ! Parameters
      !----------------

      ! Dry density of soil ! particles (excluding pores) [kg/m3]
      REAL*8,  PARAMETER   :: DNS_PRT_SFC     = 2650.0d0

      ! liq. H2O density [kg/m3]
      REAL*8,  PARAMETER   :: DNS_H2O_LQD_STD = 1000.0d0

      !-----------------
      ! Local variables
      !-----------------

      ! Longitude index
      INTEGER              :: LON_IDX

      ! [kg m-3] Bulk density of dry surface soil
      REAL*8               :: DNS_BLK_DRY(IIPAR)

      !=================================================================
      ! VWC2GWC begins here!
      !=================================================================
      GWC_SFC(:)     = 0.0D0
      DNS_BLK_DRY(:) = 0.0D0

      ! Loop over longitudes
      DO LON_IDX = 1, LONLOOP

         ! If this is a mobilization candidate then...
         IF ( FLG_MBL(LON_IDX) ) THEN

            ! Assume volume of air pores when dry equals saturated VWC
            ! This implies air pores are completely filled by water in
            ! saturated soil

            ! Bulk density of dry surface soil  [kg m-3]
            DNS_BLK_DRY(LON_IDX) = DNS_PRT_SFC
     &                           * ( 1.0d0 - VWC_SAT(LON_IDX) )

            ! Gravimetric water content [ kg kg-1]
            GWC_SFC(LON_IDX) = SOI_MST(LON_IDX)
     &                       * DNS_H2O_LQD_STD
     &                       / DNS_BLK_DRY(LON_IDX)

         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE VWC2GWC

!------------------------------------------------------------------------------

      SUBROUTINE FRC_THR_NCR_WTR_GET( FLG_MBL,     FRC_THR_NCR_WTR,
     &                                MSS_FRC_CLY, GWC_SFC )
!
!******************************************************************************
!  Subroutine FRC_THR_NCR_WTR_GET computes the factor by which soil moisture
!  increases threshold friction velocity. This parameterization is based on
!  FMB99. Zender et al., exp. (5). (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) FLG_MBL         (LOGICAL) : Mobilization candidate flag  [flags   ]
!  (3 ) MSS_FRC_CLY     (REAL*8 ) : Mass fraction of clay        [fraction]
!  (4 ) GWC_SFC         (REAL*8 ) : Gravimetric water content    [kg/kg   ]
!
!  Arguments as Output:
!  ===========================================================================
!  (2 ) FRC_THR_NCR_WTR (REAL*8 ) : Factor by which moisture increases
!                                    threshold friction velocity [fraction]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!
      USE CMN_SIZE_MOD

      ! Arguments
      LOGICAL, INTENT(IN)  :: FLG_MBL(IIPAR)
      REAL*8,  INTENT(IN)  :: MSS_FRC_CLY(IIPAR)
      REAL*8,  INTENT(IN)  :: GWC_SFC(IIPAR)
      REAL*8,  INTENT(OUT) :: FRC_THR_NCR_WTR(IIPAR)

      ! local variables
      INTEGER              :: LON_IDX        ! [idx] Counting index
      REAL*8               :: GWC_THR(IIPAR) ! [kg/kg] Threshold GWC

      !=================================================================
      ! FRC_THR_NCR_WTR_GET begins here!
      !=================================================================

      ! Initialize
      frc_thr_ncr_wtr(:) = 1.0D0
      gwc_thr(:)         = 0.0D0

      ! Loop over longitudes
      DO LON_IDX = 1, LONLOOP

         ! If this is a candidate for mobilization...
         IF ( FLG_MBL(LON_IDX) ) THEN

            !===========================================================
            ! Adjust threshold velocity for inhibition by moisture
            ! frc_thr_ncr_wtr(lon_idx)=exp(22.7D0*vwc_sfc(lon_idx))
            ! [frc] SRL96
            !
            ! Compute threshold soil moisture based on clay content
            ! GWC_THR=MSS_FRC_CLY*(0.17D0+0.14D0*MSS_FRC_CLY) [m3/m3]
            ! FMB99 p. 155 (14)
            !
            ! 19991105 remove factor of mss_frc_cly from gwc_thr to
            ! improve large scale behavior.
            !===========================================================

            ! [m3 m-3]
!            GWC_THR(LON_IDX) = 0.17D0 + 0.14D0 * MSS_FRC_CLY(LON_IDX)
            GWC_THR(LON_IDX) = 0.17D0 + 0.14D0 * 0.20d0

            IF ( GWC_SFC(LON_IDX) > GWC_THR(LON_IDX) )
     &           FRC_THR_NCR_WTR(LON_IDX) = SQRT(1.0D0+1.21D0
     &           * (100.0D0 * (GWC_SFC(LON_IDX)-GWC_THR(LON_IDX)))
     &           ** 0.68D0)     ! [frc] FMB99 p. 155 (15)
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE FRC_THR_NCR_WTR_GET

!------------------------------------------------------------------------------

      SUBROUTINE LND_FRC_MBL_GET( DOY,         FLG_MBL,     LAT_RDN,
     &                            LND_FRC_DRY, LND_FRC_MBL, MBL_NBR,
     &                            ORO,         SFC_TYP,     SNW_FRC,
     &                            TPT_SOI,     TPT_SOI_FRZ, VAI_DST )
!
!******************************************************************************
!  Subroutine LND_FRC_MBL_GET returns the fraction of each GEOS-CHEM grid
!  box which is suitable for dust mobilization.  This routine is called
!  by DST_MBL. (tdf, bmy, 4/5/04, 1/13/10)
!
!  The DATE is used to obtain the time-varying vegetation cover.
!  Routine currently uses latitude slice of VAI from time-dependent surface
!  boundary dataset (tdf, 10/27/03).  LAI/VAI algorithm is from CCM:lsm/phenol
!  () Bon96.  The LSM data are mid-month values, i.e., valid on the 15th of !
!  the month.!
!
!  Criterion for mobilisation candidate (tdf, 4/5/04):
!  (1) first, must be a land point, not ocean, not ice
!  (2) second, it cannot be an inland lake, wetland or ice
!  (3) modulated by vegetation type
!  (4) modulated by subgridscale wetness
!  (5) cannot be snow covered
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DOY         (REAL*8 ) : Day of year                         [1.0-366.0]
!  (3 ) LAT_RDN     (REAL*8 ) : Latitude                            [radians  ]
!  (4 ) LND_FRC_DRY (REAL*8 ) : Dry land fraction                   [fraction ]
!  (7 ) ORO         (REAL*8 ) : Orography: land/ocean/ice           [flags    ]
!  (8 ) SFC_TYP     (INTEGER) : LSM surface type (0..28)            [unitless ]
!  (9 ) SNW_FRC     (REAL*8 ) : Fraction of surface covered by snow [fraction ]
!  (10) TPT_SOI     (REAL*8 ) : Soil temperature                    [K        ]
!  (11) TPT_SOI_FRZ (REAL*8 ) : Temperature of frozen soil          [K        ]
!  (12) VAI_DST     (REAL*8 ) : Vegetation area index, one-sided    [m2/m2    ]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) FLG_MBL     (LOGICAL) : Mobilization candidate flag         [flag     ]
!  (5 ) LND_FRC_MBL (REAL*8 ) : Bare ground fraction                [fraction ]
!  (6 ) MBL_NBR     (INTEGER) : Number of mobilization candidates   [unitless ]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!  (2 ) For the GOCART source function, we don't use VAI, so set FLG_VAI_TVBDS
!         = .FALSE. and disable calls to ERROR_STOP (tdf, bmy, 1/25/07)
!  (3 ) Modification for GEOS-4 1 x 1.25 grids (lok, bmy, 1/13/10)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      USE CMN_SIZE_MOD    ! Size parameters
      USE CMN_GCTM_MOD    ! PI

      !------------------
      ! Arguments
      !------------------
      INTEGER, INTENT(IN)  :: SFC_TYP(IIPAR)
      REAL*8,  INTENT(IN)  :: DOY
      REAL*8,  INTENT(IN)  :: LAT_RDN
      REAL*8,  INTENT(IN)  :: LND_FRC_DRY(IIPAR)
      REAL*8,  INTENT(IN)  :: ORO(IIPAR)
      REAL*8,  INTENT(IN)  :: SNW_FRC(IIPAR)
      REAL*8,  INTENT(IN)  :: TPT_SOI(IIPAR)
      REAL*8,  INTENT(IN)  :: TPT_SOI_FRZ
      REAL*8,  INTENT(IN)  :: VAI_DST(IIPAR)
      INTEGER, INTENT(OUT) :: MBL_NBR
      LOGICAL, INTENT(OUT) :: FLG_MBL(IIPAR)
      REAL*8,  INTENT(OUT) :: LND_FRC_MBL(IIPAR)

      !------------------
      ! Parameters
      !------------------

      ! VAI threshold quench [m2/m2]
      REAL*8,  PARAMETER   :: VAI_MBL_THR = 0.30D0

      !------------------
      ! Local variables
      !------------------

      ! [idx] Counting index
      INTEGER              :: IDX_IDX

      ! [idx] Interpolation month, future
      INTEGER              :: IDX_MTH_GLB

      ! [idx] Interpolation month, past
      INTEGER              :: IDX_MTH_LUB

      ! [idx] Longitude index array (land)
      INTEGER              :: LND_IDX(IIPAR)

      ! [nbr] Number of land points
      INTEGER              :: LND_NBR

      ! [idx] Counting index for longitude
      INTEGER              :: LON_IDX

      ! [idx] Surface type index
      INTEGER              :: SFC_TYP_IDX

      ! [idx] Surface sub-gridscale index
      INTEGER              :: SGS_IDX

      !-------------------------------------------------------------------
      ! Prior to 1/25/07:
      ! For GOCART source function, we don't use VAI (tdf, bmy, 1/25/07)
      !
      ! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
      !
      !! [flg] Use VAI data from time-varying boundary dataset
      ! LOGICAL              :: FLG_VAI_TVBDS = .TRUE.
      !-------------------------------------------------------------------

      ! For GOCART source function, we do not use VAI (tdf, bmy, 1/25/07)
      LOGICAL              :: FLG_VAI_TVBDS = .FALSE.

      ! [flg] Add 182 days in southern hemisphere
      LOGICAL              :: FLG_SH_ADJ = .TRUE.

      ! [dgr] Latitude
      REAL*8               :: LAT_DGR

      ! [m2 m-2] Leaf + stem area index, one-sided
      REAL*8               :: VAI_SGS

      !=================================================================
      ! LND_FRC_MBL_GET begins here!
      !=================================================================

      ! Error check
      IF ( VAI_MBL_THR <= 0.0d0 ) THEN
         CALL ERROR_STOP( 'VAI_MBL_THR <= 0.0!',
     &                    'LND_FRC_MBL_GET ("dust_shao_mod.f")' )
      ENDIF

      ! Latitude (degrees)
      LAT_DGR = 180.0D0 * LAT_RDN/PI

      ! Initialize outputs
      MBL_NBR = 0

      DO LON_IDX = 1, LONLOOP
         FLG_MBL(LON_IDX) = .FALSE.
      ENDDO

      LND_FRC_MBL(:) = 0.0D0

      !=================================================================
      ! For dust mobilisation, we need to have land!  tdf 10/27/2K3
      ! Set up lnd_idx to hold the longitude indices for land
      ! Land ahoy!
      !=================================================================
      LND_NBR = 0
      DO LON_IDX = 1, LONLOOP
         IF ( ORO_IS_LND( ORO(LON_IDX)) ) THEN
            LND_NBR          = LND_NBR + 1
            LND_IDX(LND_NBR) = LON_IDX
         ENDIF
      ENDDO

      ! Much ado about nothing (no land points)
      IF ( LND_NBR == 0 ) RETURN

!-----------------------------------------------------------------------------
! Prior to 1/25/07:
! When GOCART source function is used, VAI flag is NOT used, so
! we need to disable the ERROR_STOP call (tdf, bmy, 1/25/07)
!
! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
!
!      ! Introduce error message for flg_vai_tvbds=F (VAI not used!)
!      IF ( .not. FLG_VAI_TVBDS ) THEN
!c         print *,' FLG_VAI_TVBDS is false: GOCART source function used'
!         CALL ERROR_STOP( 'FLG_VAI_TVBDS=F',
!     &                    'LND_FRC_MBL_GET ("dust_dead_mod.f")' )
!      ENDIF
!-----------------------------------------------------------------------------

      !=================================================================
      ! Only land points are possible candidates for dust mobilization
      !=================================================================

      ! Loop over land points
      DO IDX_IDX = 1, LND_NBR
         LON_IDX = LND_IDX(IDX_IDX)

         ! Store surface blend of current gridpoint
         SFC_TYP_IDX = SFC_TYP(LON_IDX)

         ! Check for wet or frozen conditions - no mobilisation allowed
         ! Surface type 1  = inland lakes & land ice
         ! Surface type 27 = wetlands
         IF ( SFC_TYP_IDX <= 1  .OR. SFC_TYP_IDX >= 27 .OR.
     &        TPT_SOI(LON_IDX) < TPT_SOI_FRZ )          THEN

              ! SET bare ground fraction to zero
              LND_FRC_MBL(LON_IDX) = 0.0D0

         ELSE

              LND_FRC_MBL(LON_IDX) = 1.0D0

         ENDIF                 ! endif normal land

         !==============================================================
         ! We have now filled "lnd_frc_mbl" the land fraction suitable
         ! for mobilisation.  Adjust for factors which constrain entire
         ! gridcell  LND_FRC_MBL modulated by LND_FRC_DRY and SNW_FRC.
         ! (tdf, 4/5/04)
         !==============================================================

         ! Take the bare ground fraction, multiply by the fraction
         ! that is dry and that is NOT covered by snow
         LND_FRC_MBL(LON_IDX) = LND_FRC_MBL(LON_IDX)
     &                        * LND_FRC_DRY(LON_IDX)
     &                        * ( 1.0D0 - SNW_FRC(LON_IDX) )

         ! Error check
         IF ( LND_FRC_MBL(lon_idx) > 1.0D0 ) THEN
            CALL ERROR_STOP( 'LND_FRC_MBL > 1!',
     &                       'LND_FRC_MBL_GET ("dust_shao_mod.f")' )
         ENDIF

         IF ( LND_FRC_MBL(LON_IDX) < 0.0D0 )   then
            CALL ERROR_STOP( 'LND_FRC_MBL < 0!',
     &                       'LND_FRC_MBL_GET ("dust_shao_mod.f")' )
         ENDIF

         ! If there is dry land in this longitude
         if ( LND_FRC_MBL(LON_IDX) > 0.0D0 ) then

            ! Set flag, we have a candidate!
            FLG_MBL(LON_IDX) = .TRUE.

            ! Increment # of candidates
            MBL_NBR          = MBL_NBR + 1
         ENDIF

      ENDDO

      ! Return to calling program
      END SUBROUTINE LND_FRC_MBL_GET

!------------------------------------------------------------------------------

      SUBROUTINE DST_TVBDS_GET( LAT_IDX, VAI_DST_OUT )
!
!******************************************************************************
!  Subroutine DST_TVBDS_GET returns a specifed latitude slice of VAI data.
!  (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LAT_IDX     (INTEGER) : Latitude index
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) VAI_DST_OUT (REAL*8 ) : Vegetation area index, 1-sided, current [m2/m2]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!
      USE CMN_SIZE_MOD  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: LAT_IDX
      REAL*8,  INTENT(OUT) :: VAI_DST_OUT(:)

      ! Local variables
      INTEGER              :: LON_IDX

      !=================================================================
      ! DST_TVBDS_GET begins here!
      !=================================================================

      ! Return lat slice of VAI [m2/m2]
      DO LON_IDX = 1, LONLOOP
         VAI_DST_OUT(LON_IDX) = VAI_DST(LON_IDX,LAT_IDX)
      ENDDO

      ! Return to calling program
      END SUBROUTINE DST_TVBDS_GET

!------------------------------------------------------------------------------

      SUBROUTINE GET_TIME_INVARIANT_DATA
!
!******************************************************************************
!  Subroutine GET_TIME_INVARIANT_DATA gets data for the SHAO model which
!  does not vary w/ time.  This routine is called from SRC_DUST_SHAO in
!  "dust_mod.f" only on the first timestep. (bmy, 4/5/04, 1/25/07)
!
!  NOTES:
!  (1 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now can read data for both GEOS & GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references "file_mod.f", "transfer_mod.f".  Also now read from
!        dust_200605 directory.  Now reads GOCART source function from a
!        separate file. (tdf, bmy, 1/25/07)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
!      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD     ! Size parameters

      ! Local variables
      INTEGER               :: I, J, IOS
      REAL*4                :: ARRAY(I2x25,J2x25,1)
      REAL*8                :: XTAU
      CHARACTER(LEN=255)    :: FILENAME

      INTEGER               :: F_SFC_TYP(I2x25,J2x25)
      INTEGER, EXTERNAL     :: CLX, CLY

      !=================================================================
      ! GET_TIME_INVARIANT_DATA begins here!
      !=================================================================

      ! Initialize data arrays
      CALL INIT_DUST_SHAO

      !=================================================================
      ! Need also to provide surface boundary information here
      ! read time-invariant boundary fields data set (labelled 1,1,1985)
      !
      ! The following time-invariant fields are read in
      !
      ! Erodibility field should be copied onto mbl_bsn_fct
      ! which is the one used by the DEAD code   Duncan 8/1/2003
      !
      ! MSS_FRC_CLY    ; mass fraction of clay:        IIPAR JJPAR
      ! MSS_FRC_SND    ; mass fraction of sand:        IIPAR JJPAR
      ! SFC_TYP        ; surface type:                 IIPAR JJPAR
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )         //
     &           'dust_200605/dst_tibds.' // GET_NAME_EXT_2D() //
     &           '.'                      // GET_RES_EXT()

      ! TAU value for reading the bpch files
      XTAU     = GET_TAU0( 1, 1, 1985 )

      ! Echo info
      IF (MASTER) THEN
      WRITE( 6, 100 ) TRIM( FILENAME )
      ENDIF
 100  FORMAT( '     - GET_TIME_INVARIANT_DATA: Reading ', a )

      !-----------------
      ! LND_FRC_DRY
      !-----------------
      CALL READ_BPCH2( FILENAME, 'DEAD-2D', 6,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), LND_FRC_DRY, 2, 1, 2 )

      !-----------------
      ! MSS_FRC_CLY
      !-----------------
      CALL READ_BPCH2( FILENAME, 'DEAD-2D', 8,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), MSS_FRC_CLY, 2, 1, 2 )

      !-----------------
      ! MSS_FRC_SND
      !-----------------
      CALL READ_BPCH2( FILENAME, 'DEAD-2D', 9,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), MSS_FRC_SND, 2, 1, 2 )

      !-----------------
      ! SFC_TYP
      !-----------------
      CALL READ_BPCH2( FILENAME, 'DEAD-2D', 10,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! NINT is not defined for REAL*8
      !CALL TRANSFER_2D( ARRAY(:,:,1), SFC_TYP )

      ! Also round off
      F_SFC_TYP = NINT( ARRAY(:,:,1) )

#ifdef MP
      CALL MPBCASTI(F_SFC_TYP,I2x25*J2x25)
#endif

      DO J = 1, LATLOOP
      DO I = 1, LONLOOP
         SFC_TYP(I,J) = F_SFC_TYP(CLX(I,J),CLY(I,J))
      ENDDO
      ENDDO

!****************   SOIL TEXTURE INFO    **********
      ! Filename
      FILENAME = TRIM( DATA_DIR )         //
     &           'dust_200605/dst_fsoil.' // GET_NAME_EXT_2D() //
     &           '.'                      // GET_RES_EXT()

      ! Echo info
      IF (MASTER) THEN
      WRITE( 6, 101 ) TRIM( FILENAME )
      ENDIF
 101  FORMAT( '  - GET_TIME_INVARIANT_DATA: Soil Texture Reading ', a )

      !-----------------
      ! FRC_01  SAND
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 21,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S01, 2, 1, 2 )

      !-----------------
      ! FRC_02  LOAMY SAND
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 22,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S02, 2, 1, 2 )

      !-----------------
      ! FRC_03  SANDY LOAM
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 23,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S03, 2, 1, 2 )

      !-----------------
      ! FRC_04  SILT LOAM
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 24,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S04, 2, 1, 2 )

      !-----------------
      ! FRC_05  SILT
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 25,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S05, 2, 1, 2 )

      !-----------------
      ! FRC_06  LOAM
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 26,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S06, 2, 1, 2 )

      !-----------------
      ! FRC_07  SANDY CLAY LOAM
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 27,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S07, 2, 1, 2 )

      !-----------------
      ! FRC_08  SILTY CLAY LOAM
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 28,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S08, 2, 1, 2 )

      !-----------------
      ! FRC_09  CLAY LOAM
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 29,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S09, 2, 1, 2 )

      !-----------------
      ! FRC_10  SANDY CLAY
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 30,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S10, 2, 1, 2 )

      !-----------------
      ! FRC_11  SILTY CLAY
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 31,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S11, 2, 1, 2 )

      !-----------------
      ! FRC_12  CLAY
      !-----------------
      CALL READ_BPCH2( FILENAME, 'SHAO-2D', 32,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      CALL TRANS_REGRID( DBLE(ARRAY), FRC_S12, 2, 1, 2 )

      ! Return to calling program
      END SUBROUTINE GET_TIME_INVARIANT_DATA

!------------------------------------------------------------------------------

      SUBROUTINE GET_MONTHLY_DATA
!
!******************************************************************************
!  Subroutine GET_MONTHLY_DATA gets data for the SHAO model which varies by
!  month.  This routine is called from SRC_DUST_DEAD in "dust_mod.f".
!  (tdf, bmy, 4/5/04, 1/25/07)
!
!  NOTES:
!  (1 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now can read data for both GEOS & GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now read from dust_200605 directory (tdf, bmy, 1/25/07)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_MONTH,       ITS_A_NEW_MONTH
!      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD        ! Size parameters

      ! Local variables
      INTEGER               :: THISMONTH
      REAL*4                :: ARRAY(I2x25,J2x25,1)
      REAL*8                :: XTAU
      CHARACTER(LEN=255)    :: FILENAME

      !=================================================================
      ! GET_MONTHLY_DATA begins here!
      !=================================================================

      ! Filename and time
      FILENAME  = TRIM( DATA_DIR )         //
     &            'dust_200605/dst_tvbds.' // GET_NAME_EXT_2D() //
     &            '.'                      // GET_RES_EXT()

      ! TAU for reading the bpch files
      THISMONTH = GET_MONTH()
      XTAU      = GET_TAU0( THISMONTH, 1, 1985 )

      ! Echo info
      IF (MASTER) THEN
      WRITE( 6, 100 ) TRIM( FILENAME )
      ENDIF
 100  FORMAT( '     - GET_MONTHLY_DATA: Reading ', a )

      !-----------------------
      ! Veg. Area Index (VAI)
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'DEAD-2D', 13,
     &                 XTAU,      I2x25,    J2x25,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANS_REGRID( DBLE(ARRAY), VAI_DST, 2, 1, 2 )

      ! Return to calling program
      END SUBROUTINE GET_MONTHLY_DATA

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DUST_SHAO
!
!******************************************************************************
!  Subroutine INIT_DUST_SHAO initializes all allocatable module arrays.
!  (tdf, bmy, 3/30/04, 1/25/07)
!
!  NOTES:
!  (1 ) Now allocate SRCE_FUNC (tdf, bmy, 1/25/07)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_DUST_SHAO begins here!
      !=================================================================
      ALLOCATE( LND_FRC_DRY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LND_FRC_DRY' )
      LND_FRC_DRY = 0d0

      ALLOCATE( MSS_FRC_CLY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MSS_FRC_CLY' )
      MSS_FRC_CLY = 0d0

      ALLOCATE( MSS_FRC_SND( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MSS_FRC_SND' )
      MSS_FRC_SND = 0d0

      ALLOCATE( SFC_TYP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SFC_TYP' )
      SFC_TYP = 0d0

      ALLOCATE( FRC_S01( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S01' )
      FRC_S01 = 0d0

      ALLOCATE( FRC_S02( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S02' )
      FRC_S02 = 0d0

      ALLOCATE( FRC_S03( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S03' )
      FRC_S03 = 0d0

      ALLOCATE( FRC_S04( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S04' )
      FRC_S04 = 0d0

      ALLOCATE( FRC_S05( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S05' )
      FRC_S05 = 0d0

      ALLOCATE( FRC_S06( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S06' )
      FRC_S06 = 0d0

      ALLOCATE( FRC_S07( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S07' )
      FRC_S07 = 0d0

      ALLOCATE( FRC_S08( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S08' )
      FRC_S08 = 0d0

      ALLOCATE( FRC_S09( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S09' )
      FRC_S09 = 0d0

      ALLOCATE( FRC_S10( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S10' )
      FRC_S10 = 0d0

      ALLOCATE( FRC_S11( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S11' )
      FRC_S11 = 0d0

      ALLOCATE( FRC_S12( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRC_S12' )
      FRC_S12 = 0d0

      ALLOCATE( TPT_GND( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TPT_GND' )
      TPT_GND = 0d0

      ALLOCATE( TPT_SOI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TPT_SOI' )
      TPT_SOI = 0d0

      ALLOCATE( VAI_DST( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VAI_DST' )
      VAI_DST = 0d0

      ALLOCATE( OROGRAPHY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OROGRAPHY' )
      OROGRAPHY = 0

      ! Bin size min diameter [m]
      ALLOCATE( DMT_MIN( NDSTBIN ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DMT_MIN' )
      DMT_MIN(1) = 0.2d-6
      DMT_MIN(2) = 2.0d-6
      DMT_MIN(3) = 3.6d-6
      DMT_MIN(4) = 6.0d-6

      ! Bin size max diameter [m]
      ALLOCATE( DMT_MAX( NDSTBIN ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DMT_MAX' )
      DMT_MAX(1) = 2.0d-6
      DMT_MAX(2) = 3.6d-6
      DMT_MAX(3) = 6.0d-6
      DMT_MAX(4) = 1.2d-5

      ! Return to calling program
      END SUBROUTINE INIT_DUST_SHAO

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DUST_SHAO
!
!******************************************************************************
!  Subroutine CLEANUP_DUST_SHAO deallocates all module variables.
!  (tdf, bmy, 3/30/04, 1/25/07)
!
!  NOTES:
!  (1 ) Now deallocate SRCE_FUNC (tdf, bmy, 1/25/07)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DUST_SHAO begins here!
      !=================================================================
      IF ( ALLOCATED( LND_FRC_DRY     ) ) DEALLOCATE( LND_FRC_DRY     )
      IF ( ALLOCATED( MSS_FRC_CLY     ) ) DEALLOCATE( MSS_FRC_CLY     )
      IF ( ALLOCATED( MSS_FRC_SND     ) ) DEALLOCATE( MSS_FRC_SND     )
      IF ( ALLOCATED( SFC_TYP         ) ) DEALLOCATE( SFC_TYP         )
      IF ( ALLOCATED( FRC_S01         ) ) DEALLOCATE( FRC_S01         )
      IF ( ALLOCATED( FRC_S02         ) ) DEALLOCATE( FRC_S02         )
      IF ( ALLOCATED( FRC_S03         ) ) DEALLOCATE( FRC_S03         )
      IF ( ALLOCATED( FRC_S04         ) ) DEALLOCATE( FRC_S04         )
      IF ( ALLOCATED( FRC_S05         ) ) DEALLOCATE( FRC_S05         )
      IF ( ALLOCATED( FRC_S06         ) ) DEALLOCATE( FRC_S06         )
      IF ( ALLOCATED( FRC_S07         ) ) DEALLOCATE( FRC_S07         )
      IF ( ALLOCATED( FRC_S08         ) ) DEALLOCATE( FRC_S08         )
      IF ( ALLOCATED( FRC_S09         ) ) DEALLOCATE( FRC_S09         )
      IF ( ALLOCATED( FRC_S10         ) ) DEALLOCATE( FRC_S10         )
      IF ( ALLOCATED( FRC_S11         ) ) DEALLOCATE( FRC_S11         )
      IF ( ALLOCATED( FRC_S12         ) ) DEALLOCATE( FRC_S12         )
      IF ( ALLOCATED( TPT_GND         ) ) DEALLOCATE( TPT_GND         )
      IF ( ALLOCATED( TPT_SOI         ) ) DEALLOCATE( TPT_SOI         )
      IF ( ALLOCATED( VAI_DST         ) ) DEALLOCATE( VAI_DST         )
      IF ( ALLOCATED( OROGRAPHY       ) ) DEALLOCATE( OROGRAPHY       )
      IF ( ALLOCATED( DMT_MIN         ) ) DEALLOCATE( DMT_MIN         )
      IF ( ALLOCATED( DMT_MAX         ) ) DEALLOCATE( DMT_MAX         )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DUST_SHAO

!------------------------------------------------------------------------------

      END MODULE DUST_SHAO_MOD
