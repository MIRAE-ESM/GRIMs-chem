!=====================================================================
!
!     SUBROUTINE INTERP_VER_F
!
!        This subroutine is for interpolating GEOS-Chem full
!        vertical pressure grid to GRIMs vertical pressure grid.
!
!                                            written by Seungun Lee
!
!=====================================================================

      SUBROUTINE INTERP_VER_F ( INDATA, OUTDATA, NLEV, FLAG )

      USE COMFVER,   ONLY : AK5, BK5
      USE CMN_SIZE_MOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NLEV
      INTEGER, INTENT(IN) :: FLAG
      REAL*8, INTENT(IN)  :: INDATA(IGLOB,JGLOB,NLEV)
      REAL*8, INTENT(OUT) :: OUTDATA(IGLOB,JGLOB,LLPAR)
      REAL*8              :: TEMPARRAY(IGLOB,JGLOB,73)
      REAL*8              :: AP(73), BP(73)
      REAL*8              :: PEDGE1(73), PEDGE2(LLPAR+1)
      REAL*8              :: FSUM(IGLOB,JGLOB,LLPAR)
      REAL*8              :: FR1, FR2
      INTEGER             :: L, ID, ID1, ID2

      AP = (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,
     &        1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,
     &        4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,
     &        7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,
     &        1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,
     &        1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,
     &        2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,
     &        2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02,
     &        1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01,
     &        7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01,
     &        4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01,
     &        1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01,
     &        9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00,
     &        4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00,
     &        1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01,
     &        6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01,
     &        2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02,
     &        6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02,
     &        1.000000d-02 /)
      
      BP = (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01,
     &        9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01,
     &        8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01,
     &        7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01,
     &        6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01,
     &        4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01,
     &        2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01,
     &        6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00 /)

      !================
      ! Initialization
      !================
      PEDGE1(:) = 0d0
      PEDGE2(:) = 0d0
      FSUM(:,:,:) = 0d0
      OUTDATA(:,:,:) = 0d0
      TEMPARRAY(:,:,:) = 0d0
      TEMPARRAY(:,:,1:NLEV) = INDATA(:,:,:)

      !===============================================
      ! Set up GEOS-5 reduced vertical pressure level
      !===============================================
      DO L=1, 73
         PEDGE1(L) = AP(L) + BP(L) * 1013.25
      ENDDO

      !======================================
      ! Set up GRIMs vertical pressure level
      !======================================
      DO L=1, LLPAR+1
         PEDGE2(L) = AK5(LLPAR+2-L) * 1d1 + BK5(LLPAR+2-L) * 1013.25
         IF (PEDGE2(L) .LE. 0) PEDGE2(L) = 1d-2
      ENDDO

      DO L=1, LLPAR

         ID1=0
         ID2=0

         !=====================================================
         ! Finding GEOS-5 reduced pressure level indices
         ! which correspond to the tip of GRIMs pressure level
         !=====================================================
         DO ID=1, 73
         IF (PEDGE2(L) .LE. PEDGE1(ID) .AND.
     &       PEDGE2(L) .GE. PEDGE1(ID+1)) THEN
            ID1=ID
         ENDIF
         IF (PEDGE2(L+1) .LE. PEDGE1(ID) .AND.
     &       PEDGE2(L+1) .GE. PEDGE1(ID+1)) THEN
            ID2=ID
         ENDIF
         ENDDO

         !=======================================================
         ! Calculate coverage at the tip of GRIMs pressure level
         !=======================================================
         IF (FLAG .EQ. 3) THEN

            ! Volume weighted
            IF (ID1 .EQ. ID2) THEN
               FR1 = (LOG(PEDGE2(L))-LOG(PEDGE2(L+1)))
     &             / (LOG(PEDGE1(ID1))-LOG(PEDGE1(ID1+1)))
               FR2 = 0d0
            ELSE
               FR1 = (LOG(PEDGE2(L))-LOG(PEDGE1(ID1+1)))
     &             / (LOG(PEDGE1(ID1))-LOG(PEDGE1(ID1+1)))
               FR2 = (LOG(PEDGE2(L+1))-LOG(PEDGE1(ID2)))
     &             / (LOG(PEDGE1(ID2+1))-LOG(PEDGE1(ID2)))
            ENDIF

         ELSEIF (FLAG .EQ. 1 .OR. FLAG .EQ. 2) THEN

            ! Mass weighted
            IF (ID1 .EQ. ID2) THEN
               FR1 = (PEDGE2(L)-PEDGE2(L+1))
     &             / (PEDGE1(ID1)-PEDGE1(ID1+1))
               FR2 = 0d0
            ELSE
               FR1 = (PEDGE2(L)-PEDGE1(ID1+1))
     &             / (PEDGE1(ID1)-PEDGE1(ID1+1))
               FR2 = (PEDGE2(L+1)-PEDGE1(ID2))
     &             / (PEDGE1(ID2+1)-PEDGE1(ID2))
            ENDIF

         ENDIF

         !============================
         ! Interpolate pressure level
         !============================
         OUTDATA(:,:,L) = OUTDATA(:,:,L) + TEMPARRAY(:,:,ID1) * FR1
         FSUM(:,:,L) = FSUM(:,:,L) + FR1

         DO ID=ID1+1, ID2-1
            OUTDATA(:,:,L) = OUTDATA(:,:,L) + TEMPARRAY(:,:,ID)
            FSUM(:,:,L) = FSUM(:,:,L) + 1d0
         ENDDO

         OUTDATA(:,:,L) = OUTDATA(:,:,L) + TEMPARRAY(:,:,ID2) * FR2
         FSUM(:,:,L) = FSUM(:,:,L) + FR2

      ENDDO

      !================================
      ! Average if it is concentration
      !================================
      IF (FLAG .NE. 1) OUTDATA(:,:,:) = OUTDATA(:,:,:) / FSUM(:,:,:)

      END SUBROUTINE INTERP_VER_F
