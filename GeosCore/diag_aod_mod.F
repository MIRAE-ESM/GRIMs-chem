!=====================================================================
!
!     MODULE DIAG_AOD_MOD
!
!        This module makes diagnostic AOD, SSA, ASM
!
!                                                       Seungun Lee
!
!=====================================================================

      MODULE DIAG_AOD_MOD

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      INTEGER, ALLOCATABLE :: SCALEAOD(:,:)
      REAL*8,  ALLOCATABLE :: DIAGAOD(:,:,:)
      REAL*8,  ALLOCATABLE :: DIAGSSA(:,:,:)
      REAL*8,  ALLOCATABLE :: DIAGASM(:,:,:)

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

      INTEGER            :: I,         J,           L
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

      ! Divide by SCALEAOD
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF (SCALEAOD(I,J) .EQ. 0) THEN
            ARRAY(I,J,L) = 0e0
         ELSE
            ARRAY(I,J,L) = DIAGAOD(I,J,L) / DBLE(SCALEAOD(I,J))
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 1,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ! SSA
      CATEGORY  = 'DIAGSSA'

      ! Divide by SCALEAOD
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF (SCALEAOD(I,J) .EQ. 0) THEN
            ARRAY(I,J,L) = 0e0
         ELSE
            ARRAY(I,J,L) = DIAGSSA(I,J,L) / DBLE(SCALEAOD(I,J))
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 1,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ! ASM
      CATEGORY  = 'DIAGASM'

      ! Divide by SCALEAOD
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF (SCALEAOD(I,J) .EQ. 0) THEN
            ARRAY(I,J,L) = 0e0
         ELSE
            ARRAY(I,J,L) = DIAGASM(I,J,L) / DBLE(SCALEAOD(I,J))
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 1,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      IF (LMASTER) THEN
      WRITE( 6, '(a)' ) '     - DIAG_AOD: Diagnostics written to bpch!'
      ENDIF

      ! Initialize
      SCALEAOD = 0
      DIAGAOD  = 0d0
      DIAGSSA  = 0d0
      DIAGASM  = 0d0

      END SUBROUTINE WRITE_DIAG_AOD

!=====================================================================

      SUBROUTINE INIT_DIAG_AOD

      ! Allocate module variables
      ALLOCATE( SCALEAOD(IIPAR,JJPAR) )
      ALLOCATE( DIAGAOD(IIPAR,JJPAR,LLPAR) )
      ALLOCATE( DIAGSSA(IIPAR,JJPAR,LLPAR) )
      ALLOCATE( DIAGASM(IIPAR,JJPAR,LLPAR) )

      ! Initialize module variables
      SCALEAOD = 0
      DIAGAOD  = 0d0
      DIAGSSA  = 0d0
      DIAGASM  = 0d0

      END SUBROUTINE INIT_DIAG_AOD

!=====================================================================

      END MODULE DIAG_AOD_MOD
