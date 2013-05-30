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
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EMISR!(IIPAR,JJPAR,NEMPARA)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EMISRN!(IIPAR,JJPAR,NOXLEVELS)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EMIST!(IIPAR,JJPAR,NEMPARA)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EMISTN!(IIPAR,JJPAR,NOXLEVELS)

      ! Rural Emissions:
      ! EMISRRN = NOx emissions into sigma levels L=1,NOXEXTENT
      ! EMISRR  = All other tracer emissions into sigma level L=1
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EMISRR! (IIPAR,JJPAR,NEMPARA+NEMPARB)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EMISRRN!(IIPAR,JJPAR,NOXEXTENT      )


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
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EMISS_BVOC!(IIPAR,JJPAR,10)


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

      REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:) :: EMISTNOX! (IIPAR,JJPAR,4,2)
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTETHE!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTCO!  (IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTPRPE!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTC3H8!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTALK4!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTC2H6!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:) :: EMISTSOX! (IIPAR,JJPAR,4,2)
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTACET!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTMEK! (IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTBENZ!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTTOLU!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTXYLE!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTC2H4!(IIPAR,JJPAR    )
      REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: EMISTC2H2!(IIPAR,JJPAR    )


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
      REAL*4, ALLOCATABLE, DIMENSION(:,:) :: FTOTCO2!(IIPAR,JJPAR)
      REAL*4, ALLOCATABLE, DIMENSION(:,:) :: FLIQCO2!(IIPAR,JJPAR)

      ! FRACO3, FRACNO          = fractions of O3, NO
      ! SAVEOH, SAVENO, SAVENO3 = array to save OH, NO, & NO3 fields 
      ! SAVENO2                 = array to save NO2 fields (rvm, 5/9/00)
      ! FRACNO2                 = fraction of NO2 (rvm, bmy, 2/27/02)
      ! SAVEHO2                 = array to save HO2 fields (rvm, bmy, 2/27/02)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: FRACO3!(IIPAR,JJPAR,LLPAR)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SAVEOH!(IIPAR,JJPAR,LLPAR)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: FRACNO!(IIPAR,JJPAR,LLPAR)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SAVENO!(IIPAR,JJPAR,LLPAR)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SAVENO2!(IIPAR,JJPAR,LLPAR)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SAVENO3!(IIPAR,JJPAR,LLPAR)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: FRACNO2!(IIPAR,JJPAR,LLPAR)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SAVEHO2!(IIPAR,JJPAR,LLPAR)

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
