!=====================================================================
!
!     MODULE GCMLINK_MOD
!
!        This module converts and saves the variables which
!        are needed in the chemistry processes.
!
!                                            written by Seungun Lee
!
!=====================================================================

      MODULE GCMLINK_MOD

      USE BPCH2_MOD,   ONLY : GET_TAU0, READ_BPCH2
      USE DAO_MOD,     ONLY : DTRAIN, GWETTOP, T, TROPP
      USE TIME_MOD,    ONLY : ITS_A_NEW_MONTH, GET_MONTH

      USE TRANS_REGRID_MOD
      USE CMN_SIZE_MOD

      IMPLICIT NONE

      REAL*8, DIMENSION(:,:),   ALLOCATABLE :: QDT
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: _DETRN

      CONTAINS

!==========================================================

      SUBROUTINE UPDATE_GCM

      REAL*4  :: ARRAY(I1x1,J1x1,1)
      REAL*8  :: XTAU
      CHARACTER(LEN=255) :: FILENAME

      ! Detrainment rate
      DTRAIN = _DETRN

      ! Top soil layer wetness [unitless]
      IF ( ITS_A_NEW_MONTH() ) THEN

         FILENAME = '/home/sul/MERRA/GWETTOP.clim.merra.1x1.bpch'

         XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )

         CALL READ_BPCH2( FILENAME, 'DAO-FLDS',    22,
     &                    XTAU,      I1x1,         J1x1,
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         CALL TRANS_REGRID( DBLE(ARRAY), GWETTOP, 1, 1, 2 )

      ENDIF

      ! Tropopause [hPa]
      CALL TPP_REICHLER( T, TROPP )
      TROPP = TROPP * 1d-2

      END SUBROUTINE UPDATE_GCM

!=====================================================================

      SUBROUTINE INIT_GCMLINK

      ALLOCATE( QDT      (IIPAR,LLPAR)       )
      ALLOCATE( _DETRN   (IIPAR,JJPAR,LLPAR) )

      QDT        = 0d0
      _DETRN     = 0d0

      END SUBROUTINE INIT_GCMLINK

!==========================================================

      END MODULE GCMLINK_MOD
