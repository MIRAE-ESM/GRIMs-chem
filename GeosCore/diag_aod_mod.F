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

      INTEGER :: SCALESW, SCALELW
      REAL*8,  ALLOCATABLE :: DIAGAOD(:,:,:)
      REAL*8,  ALLOCATABLE :: DIAGSSA(:,:,:)
      REAL*8,  ALLOCATABLE :: DIAGASM(:,:,:)
      REAL*8,  ALLOCATABLE :: AERSWHR(:,:,:,:)
      REAL*8,  ALLOCATABLE :: AERLWHR(:,:,:,:)
      REAL*8,  ALLOCATABLE :: AERSWRF(:,:,:)
      REAL*8,  ALLOCATABLE :: AERLWRF(:,:,:)
      REAL*8,  ALLOCATABLE :: CLDSWRF(:,:,:)
      REAL*8,  ALLOCATABLE :: CLDLWRF(:,:,:)

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

      !=================================================================
      ! CATEGORY = AEROPROP
      !=================================================================
      CATEGORY  = 'AEROPROP'

      ! Select proper unit string
      UNIT = 'unitless'

      !=================================================================
      ! AOD
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,L) = DIAGAOD(I,J,L) / SCALESW
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 1,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! SSA
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,L) = DIAGSSA(I,J,L) / SCALESW
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 2,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! ASM
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,L) = DIAGASM(I,J,L) / SCALESW
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 3,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! CATEGORY = AEROSWRF
      !=================================================================
      CATEGORY  = 'AEROSWRF'

      ! Select proper unit string
      UNIT = 'K/s'

      !=================================================================
      ! SWHR
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,L) = AERSWHR(I,J,L,1) / SCALESW
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 1,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ! Select proper unit string
      UNIT = 'W/m2'

      !=================================================================
      ! DRFTOA
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = AERSWRF(I,J,1) / SCALESW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 2,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! DRFSFC
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = AERSWRF(I,J,2) / SCALESW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 3,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      ! Select proper unit string
      UNIT = 'K/s'

      !=================================================================
      ! SWHRC
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,L) = AERSWHR(I,J,L,2) / SCALESW
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 4,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ! Select proper unit string
      UNIT = 'W/m2'

      !=================================================================
      ! DRFTOAC
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = AERSWRF(I,J,3) / SCALESW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 5,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! DRFSFCC
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = AERSWRF(I,J,4) / SCALESW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 6,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! CATEGORY = AEROLWRF
      !=================================================================
      CATEGORY  = 'AEROLWRF'

      ! Select proper unit string
      UNIT = 'K/s'

      !=================================================================
      ! LWHR
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,L) = AERLWHR(I,J,L,1) / SCALELW
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 1,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ! Select proper unit string
      UNIT = 'W/m2'

      !=================================================================
      ! DRFTOA
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = AERLWRF(I,J,1) / SCALELW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 2,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! DRFSFC
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = AERLWRF(I,J,2) / SCALELW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 3,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      ! Select proper unit string
      UNIT = 'K/s'

      !=================================================================
      ! LWHRC
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,L) = AERLWHR(I,J,L,2) / SCALELW
      ENDDO
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 4,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,    1,
     &            1,         1,         ARRAY(:,:,1) )

      ! Select proper unit string
      UNIT = 'W/m2'

      !=================================================================
      ! DRFTOAC
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = AERLWRF(I,J,3) / SCALELW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 5,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! DRFSFCC
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = AERLWRF(I,J,4) / SCALELW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 6,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! CATEGORY = CLDSWRF
      !=================================================================
      CATEGORY  = 'CLDSWRF'

      ! Select proper unit string
      UNIT = 'W/m2'

      !=================================================================
      ! CRFTOA
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = CLDSWRF(I,J,1) / SCALESW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 1,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! CRFSFC
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = CLDSWRF(I,J,2) / SCALESW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 2,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! CATEGORY = CLDLWRF
      !=================================================================
      CATEGORY  = 'CLDLWRF'

      ! Select proper unit string
      UNIT = 'W/m2'

      !=================================================================
      ! CRFTOA
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = CLDLWRF(I,J,1) / SCALELW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 1,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      !=================================================================
      ! CRFSFC
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ARRAY(I,J,1) = CLDLWRF(I,J,2) / SCALELW
      ENDDO
      ENDDO

      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, 2,
     &            UNIT,      DIAGb,     DIAGe,    RESERVED,
     &            IIPAR,     JJPAR,     1,        1,
     &            1,         1,         ARRAY(:,:,1) )

      IF (LMASTER) THEN
      WRITE( 6, '(a)' ) '     - DIAG_AOD: Diagnostics written to bpch!'
      ENDIF

      ! Initialize
      SCALESW = 0
      SCALELW = 0
      DIAGAOD = 0d0
      DIAGSSA = 0d0
      DIAGASM = 0d0
      AERSWHR = 0d0
      AERLWHR = 0d0
      AERSWRF = 0d0
      AERLWRF = 0d0
      CLDSWRF = 0d0
      CLDLWRF = 0d0

      END SUBROUTINE WRITE_DIAG_AOD

!=====================================================================

      SUBROUTINE INIT_DIAG_AOD

      ! Allocate module variables
      ALLOCATE( DIAGAOD(IIPAR,JJPAR,LLPAR) )
      ALLOCATE( DIAGSSA(IIPAR,JJPAR,LLPAR) )
      ALLOCATE( DIAGASM(IIPAR,JJPAR,LLPAR) )
      ALLOCATE( AERSWHR(IIPAR,JJPAR,LLPAR,2) )
      ALLOCATE( AERLWHR(IIPAR,JJPAR,LLPAR,2) )
      ALLOCATE( AERSWRF(IIPAR,JJPAR,4) )
      ALLOCATE( AERLWRF(IIPAR,JJPAR,4) )
      ALLOCATE( CLDSWRF(IIPAR,JJPAR,2) )
      ALLOCATE( CLDLWRF(IIPAR,JJPAR,2) )

      ! Initialize module variables
      SCALESW = 0
      SCALELW = 0
      DIAGAOD = 0d0
      DIAGSSA = 0d0
      DIAGASM = 0d0
      AERSWHR = 0d0
      AERLWHR = 0d0
      AERSWRF = 0d0
      AERLWRF = 0d0
      CLDSWRF = 0d0
      CLDLWRF = 0d0

      END SUBROUTINE INIT_DIAG_AOD

!=====================================================================

      END MODULE DIAG_AOD_MOD
