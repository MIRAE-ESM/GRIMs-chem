!=====================================================================
!
!     MODULE DIAG_AOD_MOD
!
!        This module saves and prints diagnostic AOD, SSA, ASM
!
!                                            written by Seungun Lee
!
!=====================================================================

      MODULE DIAG_AOD_MOD

      USE FLEXAOD, ONLY : OUTOD, OUTSSA, OUTG, NBND
      USE CMN_SIZE_MOD

      IMPLICIT NONE

      INTEGER, DIMENSION(:,:,:),   ALLOCATABLE :: SCALEAOD
      REAL*8,  DIMENSION(:,:,:,:), ALLOCATABLE :: DIAGAOD
      REAL*8,  DIMENSION(:,:,:,:), ALLOCATABLE :: DIAGSSA
      REAL*8,  DIMENSION(:,:,:,:), ALLOCATABLE :: DIAGASM

      CONTAINS

!=====================================================================

      SUBROUTINE WRITE_DIAG_AOD

      USE BPCH2_MOD, ONLY : BPCH2
      USE BPCH2_MOD, ONLY : GET_HALFPOLAR
      USE BPCH2_MOD, ONLY : GET_MODELNAME
      USE FILE_MOD,  ONLY : IU_BPCH
      USE TIME_MOD,  ONLY : GET_DIAGb
      USE TIME_MOD,  ONLY : GET_DIAGe

      USE CMN_SIZE_MOD  ! Size parameters

      INTEGER            :: I,         J,           L,      N
      INTEGER            :: CENTER180, HALFPOLAR
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4             :: LONRES,    LATRES
      REAL*8             :: DIAGb,     DIAGe
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY,  RESERVED,    UNIT

      ! Initialize
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      LATRES    = DJSIZE
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''

      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Select proper unit string
      UNIT = 'unitless'

      ! AOD
      CATEGORY  = 'DIAGAOD'

      DO N = 1, NBND

      ! Divide by SCALEAOD
      DO J = 1, LATLOOP
      DO I = 1, LONLOOP
         ARRAY(I,J,:) = DIAGAOD(I,J,:,N) / DBLE(SCALEAOD(I,J,N))
         IF (SCALEAOD(I,J,N) .EQ. 0) ARRAY(I,J,:) = 0e0
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ENDDO

      ! SSA
      CATEGORY  = 'DIAGSSA'

      DO N = 1, NBND

      ! Divide by SCALEAOD
      DO J = 1, LATLOOP
      DO I = 1, LONLOOP
         ARRAY(I,J,:) = DIAGSSA(I,J,:,N) / DBLE(SCALEAOD(I,J,N))
         IF (SCALEAOD(I,J,N) .EQ. 0) ARRAY(I,J,:) = 0e0
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ENDDO

      ! ASM
      CATEGORY  = 'DIAGASM'

      DO N = 1, NBND

      ! Divide by SCALEAOD
      DO J = 1, LATLOOP
      DO I = 1, LONLOOP
         ARRAY(I,J,:) = DIAGASM(I,J,:,N) / DBLE(SCALEAOD(I,J,N))
         IF (SCALEAOD(I,J,N) .EQ. 0) ARRAY(I,J,:) = 0e0
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ENDDO

      IF (MASTER) THEN
      WRITE( 6, '(a)' ) '     - DIAG_AOD: Diagnostics written to bpch!'
      ENDIF

      SCALEAOD = 0
      DIAGAOD = 0d0
      DIAGSSA = 0d0
      DIAGASM = 0d0

      END SUBROUTINE WRITE_DIAG_AOD

!=====================================================================

      SUBROUTINE INIT_DIAG_AOD

      ALLOCATE( SCALEAOD(IIPAR,JJPAR,NBND) )
      ALLOCATE( DIAGAOD(IIPAR,JJPAR,LLPAR,NBND) )
      ALLOCATE( DIAGSSA(IIPAR,JJPAR,LLPAR,NBND) )
      ALLOCATE( DIAGASM(IIPAR,JJPAR,LLPAR,NBND) )

      SCALEAOD = 0
      DIAGAOD = 0d0
      DIAGSSA = 0d0
      DIAGASM = 0d0

      END SUBROUTINE INIT_DIAG_AOD

!=====================================================================

      END MODULE DIAG_AOD_MOD
