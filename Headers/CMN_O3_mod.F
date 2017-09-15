!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_O3_mod
!
! !DESCRIPTION: Common blocks for anthro emissions (via SMVGEAR!)
!
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_O3_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR,     JJPAR,     LLPAR
      USE CMN_SIZE_MOD, ONLY : NEMPARA,   NEMPARB
      USE CMN_SIZE_MOD, ONLY : NOXLEVELS, NOXEXTENT
      
      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! Rural Emissions: EMISRN = NOx (1:NOXLEVELS), EMISR = all other tracers
      ! Total Emissions: EMISTN = NOx (1:NOXLEVELS), EMIST = all other tracers
      REAL*8, ALLOCATABLE :: EMISR(:,:,:)
      REAL*8, ALLOCATABLE :: EMISRN(:,:,:)
      REAL*8, ALLOCATABLE :: EMIST(:,:,:)
      REAL*8, ALLOCATABLE :: EMISTN(:,:,:)

      ! Rural Emissions:
      ! EMISRRN = NOx emissions into sigma levels L=1,NOXEXTENT
      ! EMISRR  = All other tracer emissions into sigma level L=1
      REAL*8, ALLOCATABLE :: EMISRR(:,:,:)
      REAL*8, ALLOCATABLE :: EMISRRN(:,:,:)


      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !
      ! New biogenic VOC emissions (mpb,2009) 
      !
      !     ---------------------------
      !     Species            | Order 
      !     ---------------------------
      !     Isoprene           = 1
      !     Total Monoterpenes = 2 
      !     MBO                = 3 
      !     Alpha-Pinene       = 4
      !     Beta-Pinene        = 5
      !     Limonene           = 6
      !     Sabinene           = 7
      !     Mycrene            = 8
      !     3-Carene           = 9
      !     Ocimene            = 10
      !     ---------------------------	

      ! Define common block	
      REAL*8, ALLOCATABLE :: EMISS_BVOC(:,:,:)


      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Arrays to read emissions from updated merge file : 
      !        NOx, CO, PRPE, C3H8, ALK4, C2H6, ACET, MEK
      ! NOTE: ALD2 is not emitted in GEIA so we don't need an array for
      !       it below...but it is emitted in EMEP.  It will be saved
      !       into the EMISRR array for SMVGEAR. (bdf, bmy, 11/1/05)
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Arrays to read emissions from updated merge file : 
      !        NOx, CO, PRPE, C3H8, ALK4, C2H6, ACET, MEK
      ! NOTE: ALD2 is not emitted in GEIA so we don't need an array for
      !       it below...but it is emitted in EMEP.  It will be saved

      REAL*4, ALLOCATABLE :: EMISTNOX(:,:,:,:)
      REAL*4, ALLOCATABLE :: EMISTETHE(:,:)
      REAL*4, ALLOCATABLE :: EMISTCO(:,:)
      REAL*4, ALLOCATABLE :: EMISTPRPE(:,:)
      REAL*4, ALLOCATABLE :: EMISTC3H8(:,:)
      REAL*4, ALLOCATABLE :: EMISTALK4(:,:)
      REAL*4, ALLOCATABLE :: EMISTC2H6(:,:)
      REAL*4, ALLOCATABLE :: EMISTSOX(:,:,:,:)
      REAL*4, ALLOCATABLE :: EMISTACET(:,:)
      REAL*4, ALLOCATABLE :: EMISTMEK(:,:)
      REAL*4, ALLOCATABLE :: EMISTBENZ(:,:)
      REAL*4, ALLOCATABLE :: EMISTTOLU(:,:)
      REAL*4, ALLOCATABLE :: EMISTXYLE(:,:)
      REAL*4, ALLOCATABLE :: EMISTC2H4(:,:)
      REAL*4, ALLOCATABLE :: EMISTC2H2(:,:)


      ! Time of day and weekday/weekend scale factors
      ! NOTE: Now SCNR89 is (3,3) because of the weekday scale factor!!!
      REAL*8  :: TODH(6)
      REAL*8  :: TODN(6)
      REAL*8  :: TODB(6)
      REAL*8  :: SCNR89(3,3)


      ! IFSCLYR = Year to use for scaling fossil fuel emissions 
      ! (1985 = no scaling        !)
      INTEGER :: FSCALYR


      ! FTOTCO2 = yearly scale factors based on Total  Fuel CO2 emissions 
      ! FLIQCO2 = yearly scale factors based on Liquid Fuel CO2 emissions
      REAL*4, ALLOCATABLE :: FTOTCO2(:,:)
      REAL*4, ALLOCATABLE :: FLIQCO2(:,:)

      ! FRACO3, FRACNO          = fractions of O3, NO
      ! SAVEOH, SAVENO, SAVENO3 = array to save OH, NO, & NO3 fields 
      ! SAVENO2                 = array to save NO2 fields (rvm, 5/9/00)
      ! FRACNO2                 = fraction of NO2 (rvm, bmy, 2/27/02)
      ! SAVEHO2                 = array to save HO2 fields (rvm, bmy, 2/27/02)
      REAL*8, ALLOCATABLE :: FRACO3(:,:,:)
      REAL*8, ALLOCATABLE :: SAVEOH(:,:,:)
      REAL*8, ALLOCATABLE :: FRACNO(:,:,:)
      REAL*8, ALLOCATABLE :: SAVENO(:,:,:)
      REAL*8, ALLOCATABLE :: SAVENO2(:,:,:)
      REAL*8, ALLOCATABLE :: SAVENO3(:,:,:)
      REAL*8, ALLOCATABLE :: FRACNO2(:,:,:)
      REAL*8, ALLOCATABLE :: SAVEHO2(:,:,:)

      CONTAINS

      SUBROUTINE INIT_CMN_O3

      ALLOCATE ( EMISR(IIPAR,JJPAR,NEMPARA) )
      ALLOCATE ( EMISRN(IIPAR,JJPAR,NOXLEVELS) )
      ALLOCATE ( EMIST(IIPAR,JJPAR,NEMPARA) )
      ALLOCATE ( EMISTN(IIPAR,JJPAR,NOXLEVELS) )
      ALLOCATE ( EMISRR(IIPAR,JJPAR,NEMPARA+NEMPARB) )
      ALLOCATE ( EMISRRN(IIPAR,JJPAR,NOXEXTENT) )
      ALLOCATE ( EMISS_BVOC(IIPAR,JJPAR,10) )
      ALLOCATE ( EMISTNOX(IIPAR,JJPAR,4,2) )
      ALLOCATE ( EMISTETHE(IIPAR,JJPAR) )
      ALLOCATE ( EMISTCO(IIPAR,JJPAR) )
      ALLOCATE ( EMISTPRPE(IIPAR,JJPAR) )
      ALLOCATE ( EMISTC3H8(IIPAR,JJPAR) )
      ALLOCATE ( EMISTALK4(IIPAR,JJPAR) )
      ALLOCATE ( EMISTC2H6(IIPAR,JJPAR) )
      ALLOCATE ( EMISTSOX(IIPAR,JJPAR,4,2) )
      ALLOCATE ( EMISTACET(IIPAR,JJPAR) )
      ALLOCATE ( EMISTMEK(IIPAR,JJPAR) )
      ALLOCATE ( EMISTBENZ(IIPAR,JJPAR) )
      ALLOCATE ( EMISTTOLU(IIPAR,JJPAR) )
      ALLOCATE ( EMISTXYLE(IIPAR,JJPAR) )
      ALLOCATE ( EMISTC2H4(IIPAR,JJPAR) )
      ALLOCATE ( EMISTC2H2(IIPAR,JJPAR) )
      ALLOCATE ( FTOTCO2(IIPAR,JJPAR) )
      ALLOCATE ( FLIQCO2(IIPAR,JJPAR) )
      ALLOCATE ( FRACO3(IIPAR,JJPAR,LLPAR) )
      ALLOCATE ( SAVEOH(IIPAR,JJPAR,LLPAR) )
      ALLOCATE ( FRACNO(IIPAR,JJPAR,LLPAR) )
      ALLOCATE ( SAVENO(IIPAR,JJPAR,LLPAR) )
      ALLOCATE ( SAVENO2(IIPAR,JJPAR,LLPAR) )
      ALLOCATE ( SAVENO3(IIPAR,JJPAR,LLPAR) )
      ALLOCATE ( FRACNO2(IIPAR,JJPAR,LLPAR) )
      ALLOCATE ( SAVEHO2(IIPAR,JJPAR,LLPAR) )

      EMISR      = 0d0
      EMISRN     = 0d0
      EMIST      = 0d0
      EMISTN     = 0d0
      EMISRR     = 0d0
      EMISRRN    = 0d0
      EMISS_BVOC = 0d0
      EMISTNOX   = 0e0
      EMISTETHE  = 0e0
      EMISTCO    = 0e0
      EMISTPRPE  = 0e0
      EMISTC3H8  = 0e0
      EMISTALK4  = 0e0
      EMISTC2H6  = 0e0
      EMISTSOX   = 0e0
      EMISTACET  = 0e0
      EMISTMEK   = 0e0
      EMISTBENZ  = 0e0
      EMISTTOLU  = 0e0
      EMISTXYLE  = 0e0
      EMISTC2H4  = 0e0
      EMISTC2H2  = 0e0
      FTOTCO2    = 0e0
      FLIQCO2    = 0e0
      FRACO3     = 0d0
      SAVEOH     = 0d0
      FRACNO     = 0d0
      SAVENO     = 0d0
      SAVENO2    = 0d0
      SAVENO3    = 0d0
      FRACNO2    = 0d0
      SAVEHO2    = 0d0

      END SUBROUTINE INIT_CMN_O3
!
! !REMARKS:
!  NOTE: Now NEMPARA = max no. of anthropogenic emissions
!            NEMPARB = max no. of biogenic emissions
!                                                                             .
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS NOW DEPRECATED.  IT WILL BE REMOVED WHEN  %%%
!  %%%  THE GEOS-CHEM EMISSIONS MODULE IS UPDATED (SOMETIME SOON)      %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE CMN_O3_MOD
!EOC
