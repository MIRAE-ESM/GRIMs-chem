!=====================================================================
!
!     MODULE TRANSPORT_MOD
!
!        This module transports the chemical tracers using
!        semi-Lagrangian mass advection scheme.
!
!                                            written by Seungun Lee
!
!=====================================================================

      MODULE TRANSPORT_MOD

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      CONTAINS

!=====================================================================

      SUBROUTINE DO_TRANSPORT

      USE COMCHEM,    ONLY : TNLEVS, TRANSFORM_3D
      USE COMFVER,    ONLY : DELTIM
      USE NISLQ,      ONLY : SLQ_PSFC2, SLQ_U2, SLQ_V2, SLQ_W2
      USE TRACER_MOD, ONLY : N_TRACERS, STT, STT1, STT2

      REAL*8 :: WORK1(IIPAR,TNLEVS,JJPAR)
      REAL*8 :: WORK2(IIPAR,TNLEVS,JJPAR)
      REAL*8 :: STT_REV(IIPAR,TNLEVS,JJPAR)
      INTEGER :: I, J, K, N

      WORK1 = 0d0
      WORK2 = 0d0
      STT_REV = 0d0

      STT1 = STT2
      STT2 = STT

      DO N=1, N_TRACERS
         CALL TRANSFORM_STT
     &        ( STT1(:,:,:,N), WORK1(:,(N-1)*LLPAR+1:N*LLPAR,:) )
      ENDDO

      FORALL (I=1:IIPAR, K=1:TNLEVS, J=1:JJPAR)
     &       STT_REV(I,K,J) = WORK1(I,TNLEVS+1-K,J)
      CALL TRACER_MASS_ADVECT
     &     (DELTIM, SLQ_PSFC2, SLQ_U2, SLQ_V2, SLQ_W2, STT_REV, WORK2)

      DO N=1, N_TRACERS
         CALL TRANSFORM_3D
     &        ( WORK2(:,(N-1)*LLPAR+1:N*LLPAR,:), STT(:,:,:,N) )
      ENDDO

      END SUBROUTINE DO_TRANSPORT

!=====================================================================

      SUBROUTINE TRANSFORM_STT ( INDATA, OUTDATA )

      REAL*8, INTENT(IN)  :: INDATA(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: OUTDATA(IIPAR,LLPAR,JJPAR)
      INTEGER :: L

      DO L=1, LLPAR
         OUTDATA(:,L,:) = INDATA(:,:,L)
      ENDDO

      END SUBROUTINE TRANSFORM_STT

!=====================================================================

      END MODULE TRANSPORT_MOD
