!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_DEP_mod
!
! !DESCRIPTION: Common blocks for dry deposition.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_DEP_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : I2x25, J2x25, IIPAR, JJPAR, NTYPE

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! IREG   - Number of landtypes in grid square (I,J)
      INTEGER, ALLOCATABLE :: IREG_2x25(:,:)                                

      ! ILAND  - Land type ID in grid square (I,J) for IREG landtypes
      INTEGER, ALLOCATABLE :: ILAND_2x25(:,:,:)

      ! IUSE   - Fraction ((per mil) of gridbox area occupied by land type
      INTEGER, ALLOCATABLE :: IUSE_2x25(:,:,:)                          

      ! XLAI    - Leaf Area Index of land type for current MONTH
      REAL*8,  ALLOCATABLE :: XLAI_2x25(:,:,:)

      ! XLAI2   - Leaf Area Index of land type for following MONTH
      REAL*8,  ALLOCATABLE :: XLAI2_2x25(:,:,:)

      ! IREG   - Number of landtypes in grid square (I,J)
      INTEGER, ALLOCATABLE :: IREG(:,:)                                

      ! ILAND  - Land type ID in grid square (I,J) for IREG landtypes
      INTEGER, ALLOCATABLE :: ILAND(:,:,:)

      ! IUSE   - Fraction ((per mil) of gridbox area occupied by land type
      INTEGER, ALLOCATABLE :: IUSE(:,:,:)                          

      ! Fraction of land in grid box
      REAL*8,  ALLOCATABLE :: FRCLND(:,:)

      ! XLAI    - Leaf Area Index of land type for current MONTH
      REAL*8,  ALLOCATABLE :: XLAI(:,:,:)

      ! XLAI2   - Leaf Area Index of land type for following MONTH
      REAL*8,  ALLOCATABLE :: XLAI2(:,:,:)

      CONTAINS

      SUBROUTINE INIT_CMN_DEP

      ALLOCATE ( IREG_2x25(I2x25,J2x25) )
      ALLOCATE ( ILAND_2x25(I2x25,J2x25,NTYPE) )
      ALLOCATE ( IUSE_2x25(I2x25,J2x25,NTYPE) )
      ALLOCATE ( XLAI_2x25(I2x25,J2x25,NTYPE) )
      ALLOCATE ( XLAI2_2x25(I2x25,J2x25,NTYPE) )
      ALLOCATE ( IREG(IIPAR,JJPAR) )
      ALLOCATE ( ILAND(IIPAR,JJPAR,NTYPE) )
      ALLOCATE ( IUSE(IIPAR,JJPAR,NTYPE) )
      ALLOCATE ( FRCLND(IIPAR,JJPAR) )
      ALLOCATE ( XLAI(IIPAR,JJPAR,NTYPE) )
      ALLOCATE ( XLAI2(IIPAR,JJPAR,NTYPE) )

      IREG_2x25  = 0
      ILAND_2x25 = 0
      IUSE_2x25  = 0
      XLAI_2x25  = 0d0
      XLAI2_2x25 = 0d0
      IREG   = 0
      ILAND  = 0
      IUSE   = 0
      FRCLND = 0d0
      XLAI   = 0d0
      XLAI2  = 0d0

      END SUBROUTINE INIT_CMN_DEP
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS DEPRECATED.  IT WILL BE REMOVED WHEN THE  %%%
!  %%%  GEIA BIOGENIC EMISSIONS ARE REMOVED FROM GEOS-CHEM (v9-01-03)  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long     - Converted to Module from Header file
!  25 Aug 2011 - R. Yantosca - Change IGLOB,JGLOB to IIPAR,JJPAR
!EOP
!-----------------------------------------------------------------------------
!BOC
	END MODULE CMN_DEP_MOD
!EOC
