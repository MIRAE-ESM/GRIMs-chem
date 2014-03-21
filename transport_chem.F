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

      USE COMCHEM,      ONLY : NLEVS, TRANSFORM_3D, DP1, DP2, DP3
      USE COMFSPEC_VR,  ONLY : QM, Q, Z
      USE COMFVER,      ONLY : DELTIM
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE NISLQ,        ONLY : SLQ_PSFC2, SLQ_U2, SLQ_V2, SLQ_W2
      USE TRACER_MOD,   ONLY : N_TRACERS, STT1, STT, TCVV

      USE CMN_SIZE_MOD
      USE CMN_GCTM_MOD

      IMPLICIT NONE

      REAL*8 :: WORK1(IIPAR,NLEVS,JJPAR)
      REAL*8 :: WORK2(IIPAR,NLEVS,JJPAR)
      REAL*8 :: WORK3(IGLOB*2,JGLOB/2,LLPAR)
      REAL*8 :: WORK4(IGLOB*2,JGLOB/2,LLPAR)
      REAL*8 :: STT_REV(IIPAR,NLEVS,JJPAR)
      REAL*8 :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8 :: TR_B(IIPAR,JJPAR,LLPAR,N_TRACERS)
      REAL*8 :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8 :: TR_A(IIPAR,JJPAR,LLPAR)
      REAL*8 :: TR_DIFF, G0_1000, A_M2
      INTEGER :: I, J, K, L, N
      LOGICAL, SAVE :: FIRST = .TRUE.

      G0_1000 = 1000d0 / G0

      WORK1   = 0d0
      WORK2   = 0d0
      STT_REV = 0d0
      AD_B    = 0d0
      TR_B    = 0d0
      AD_A    = 0d0
      TR_A    = 0d0

      IF (FIRST) THEN
         CALL NISLQ_DP(QM,DP1)
         CALL NISLQ_DP(Q,DP2)
         CALL NISLQ_DP(Z,DP3)
         FIRST = .FALSE.
      ELSE
         DP1 = DP2
         DP2 = DP3
         CALL NISLQ_DP(Z,DP3)
      ENDIF

      !=================================================================
      ! Get the air and tracer mass before advection
      !=================================================================

      ! Airmass [kg] before transport
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         A_M2 = GET_AREA_M2(I,J)
         AD_B(I,J,L) = DP1(I,LLPAR+1-L,J) * G0_1000 * A_M2
      ENDDO
      ENDDO
      ENDDO

      ! Tracer mass [kg] before transport
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         TR_B(I,J,L,N) = STT1(I,J,L,N) * AD_B(I,J,L) / TCVV(N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! Do the advection
      !=================================================================

      DO N=1, N_TRACERS
         CALL TRANSFORM_STT
     &        ( STT1(:,:,:,N), WORK1(:,(N-1)*LLPAR+1:N*LLPAR,:) )
      ENDDO

      FORALL (I=1:IIPAR, K=1:NLEVS, J=1:JJPAR)
     &       STT_REV(I,K,J) = WORK1(I,NLEVS+1-K,J)
      CALL NISLQ_CHEM_ADVECT
     &     (DELTIM, SLQ_PSFC2, SLQ_U2, SLQ_V2, SLQ_W2, STT_REV, WORK2)

      DO N=1, N_TRACERS
         CALL TRANSFORM_3D
     &        ( WORK2(:,(N-1)*LLPAR+1:N*LLPAR,:), STT(:,:,:,N) )
      ENDDO

      !=================================================================
      ! Ensure mass conservation
      !=================================================================

      DO N = 1, N_TRACERS

         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Air mass [kg] after transport
            A_M2 = GET_AREA_M2(I,J)
            AD_A(I,J,L) = DP3(I,LLPAR+1-L,J) * G0_1000 * A_M2

            ! Tracer mass [kg] after transport
            TR_A(I,J,L) = STT(I,J,L,N) * AD_A(I,J,L) / TCVV(N)

         ENDDO
         ENDDO
         ENDDO

         CALL MPGP2F(TR_B(:,:,:,N),IIPAR,JJPAR,
     &               WORK3,IGLOB*2,JGLOB/2,LLPAR)
         CALL MPGP2F(TR_A,IIPAR,JJPAR,WORK4,IGLOB*2,JGLOB/2,LLPAR)

         IF (MASTER) THEN

         ! Residual mass difference: before / after
         TR_DIFF = SUM( WORK3 ) / SUM( WORK4 )

         ENDIF

         CALL MPBCASTR(TR_DIFF,1)

         ! Add mass difference back to STT
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            TR_A(I,J,L) = TR_A(I,J,L) * TR_DIFF

            STT(I,J,L,N) = TR_A(I,J,L) / AD_A(I,J,L) * TCVV(N)

            STT(I,J,L,N) = MAX( STT(I,J,L,N), 0d0 )

         ENDDO
         ENDDO
         ENDDO

      ENDDO

      END SUBROUTINE TRANSPORT_CHEM

!=====================================================================

      SUBROUTINE TRANSFORM_STT ( INDATA, OUTDATA )

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      REAL*8, INTENT(IN)  :: INDATA(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: OUTDATA(IIPAR,LLPAR,JJPAR)
      INTEGER :: L

      DO L=1, LLPAR
         OUTDATA(:,L,:) = INDATA(:,:,L)
      ENDDO

      END SUBROUTINE TRANSFORM_STT
