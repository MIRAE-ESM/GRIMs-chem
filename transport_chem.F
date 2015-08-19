!=====================================================================
!
!     SUBROUTINE TRANSPORT_CHEM
!
!        This subroutine transports the chemical tracers using
!        semi-Lagrangian mass advection scheme.
!
!                                            written by Seungun Lee
!
!=====================================================================

      SUBROUTINE TRANSPORT_CHEM

      USE COMFVER,      ONLY : DELTIM
      USE DAO_MOD,      ONLY : AD
      USE NISLQ,        ONLY : SLQ_PSFC2, SLQ_U2, SLQ_V2, SLQ_W2
      USE TRACER_MOD,   ONLY : N_TRACERS, STT1, STT

      USE CMN_SIZE_MOD

      IMPLICIT NONE
#     include "define.h"

      INTEGER :: I, J, L, N
      REAL*8 :: WORK1(IIPAR,N_TRACERS*LLPAR,JJPAR)
      REAL*8 :: WORK2(IIPAR,N_TRACERS*LLPAR,JJPAR)
      REAL*8 :: TR1(IIPAR,JJPAR,LLPAR)
      REAL*8 :: TR2(IIPAR,JJPAR,LLPAR)
      REAL*8 :: MASS1, MASS2, FAC
      REAL*8, EXTERNAL :: MAKE_SUM

      !=================================================================
      ! Do the advection
      !=================================================================

      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         WORK1(I,(N_TRACERS+1-N)*LLPAR+1-L,J) = STT1(I,J,L,N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL NISLQ_CHEM_ADVECT
     &     (DELTIM, SLQ_PSFC2, SLQ_U2, SLQ_V2, SLQ_W2, WORK1, WORK2)

      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         STT(I,J,L,N) = WORK2(I,(N-1)*LLPAR+L,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! Ensure mass conservation
      !=================================================================

      DO N = 1, N_TRACERS

         TR1 = STT1(:,:,:,N) * AD(:,:,:)
         TR2 = STT (:,:,:,N) * AD(:,:,:)
         MASS1 = MAKE_SUM( TR1, LLPAR )
         MASS2 = MAKE_SUM( TR2, LLPAR )

         IF ( MASS2 .EQ. 0d0 ) THEN
            FAC = 1d0
         ELSE
            FAC = MASS1 / MASS2
         ENDIF

#ifdef MP
         CALL MPBCASTR(FAC,1)
#endif

         STT(:,:,:,N) = STT(:,:,:,N) * FAC

      ENDDO

      END SUBROUTINE TRANSPORT_CHEM
