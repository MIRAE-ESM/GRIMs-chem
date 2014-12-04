#include "define.h"
!=====================================================================
!
!     MODULE RESTART
!
!        This module makes restart file for chemistry module.
!
!                                            written by Seungun Lee
!
!=====================================================================

      MODULE RESTART_MOD

      USE DAO_MOD,      ONLY : AIRDEN, BXHEIGHT, RH
      USE MODULE_TRANS, ONLY : DYN_TRANS2MODEL_GRID
      USE MODULE_TRANS, ONLY : DYN_TRANS2OUTPUT_GRID
      USE TRACER_MOD,   ONLY : STT1, STT2, N_TRACERS

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      CONTAINS

!======================================================================

      SUBROUTINE READ_RESTART

      REAL*8 :: WORK(IIPAR,JJPAR,LLPAR)
      INTEGER :: I, J, L, N, IOS

      IF (MASTER) THEN
      OPEN(9999, FILE='chemi', FORM='UNFORMATTED', STATUS='OLD',
     &     ACCESS='SEQUENTIAL', IOSTAT=IOS)
      ENDIF

#ifdef MP
      CALL MPBCASTI(IOS,1)
#endif

      IF (IOS .NE. 0) RETURN

      DO N=1,N_TRACERS
         CALL GET_RESTART_DATA( STT1(:,:,:,N) )
      ENDDO

      IF (MASTER) THEN
      CLOSE(9999)
      WRITE(6,*) '    - RESTART_MOD: Reading chemi'
      ENDIF

      IF (MASTER) THEN
      OPEN(9999, FILE='chemidt', FORM='UNFORMATTED', STATUS='OLD',
     &     ACCESS='SEQUENTIAL', IOSTAT=IOS)
      ENDIF

#ifdef MP
      CALL MPBCASTI(IOS,1)
#endif

      IF (IOS .NE. 0) THEN
         STT2 = STT1
         RETURN
      ENDIF

      DO N=1,N_TRACERS
         CALL GET_RESTART_DATA( STT2(:,:,:,N) )
      ENDDO

      IF (MASTER) THEN
      CLOSE(9999)
      WRITE(6,*) '    - RESTART_MOD: Reading chemidt'
      ENDIF

      IF (MASTER) THEN
      OPEN(9999, FILE='bxhght', FORM='UNFORMATTED', STATUS='OLD',
     &     ACCESS='SEQUENTIAL', IOSTAT=IOS)
      ENDIF

#ifdef MP
      CALL MPBCASTI(IOS,1)
#endif

      IF (IOS .NE. 0) RETURN

      CALL GET_RESTART_DATA( RH )
      CALL GET_RESTART_DATA( BXHEIGHT )
      CALL GET_RESTART_DATA( WORK )

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AIRDEN(L,I,J) = WORK(I,J,L)
      ENDDO
      ENDDO
      ENDDO

      IF (MASTER) THEN
      CLOSE(9999)
      WRITE(6,*) '    - RESTART_MOD: Reading bxhght'
      ENDIF

      END SUBROUTINE READ_RESTART

!======================================================================

      SUBROUTINE WRITE_RESTART

      REAL*8 :: WORK(IIPAR,JJPAR,LLPAR)
      INTEGER :: I, J, L, N

      IF (MASTER) THEN
      OPEN(9999, FILE='chemi', FORM='UNFORMATTED', STATUS='REPLACE',
     &     ACCESS='SEQUENTIAL')
      ENDIF

      DO N=1,N_TRACERS
         CALL PUT_RESTART_DATA( STT1(:,:,:,N) )
      ENDDO

      IF (MASTER) THEN
      CLOSE(9999)
      WRITE(6,*) '    - RESTART_MOD: Restart chemi written!'
      ENDIF

      IF (MASTER) THEN
      OPEN(9999, FILE='chemidt', FORM='UNFORMATTED', STATUS='REPLACE',
     &     ACCESS='SEQUENTIAL')
      ENDIF

      DO N=1,N_TRACERS
         CALL PUT_RESTART_DATA( STT2(:,:,:,N) )
      ENDDO

      IF (MASTER) THEN
      CLOSE(9999)
      WRITE(6,*) '    - RESTART_MOD: Restart chemidt written!'
      ENDIF

      IF (MASTER) THEN
      OPEN(9999, FILE='bxhght', FORM='UNFORMATTED', STATUS='REPLACE',
     &     ACCESS='SEQUENTIAL')
      ENDIF

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         WORK(I,J,L) = AIRDEN(L,I,J)
      ENDDO
      ENDDO
      ENDDO

      CALL PUT_RESTART_DATA( RH )
      CALL PUT_RESTART_DATA( BXHEIGHT )
      CALL PUT_RESTART_DATA( WORK )

      IF (MASTER) THEN
      CLOSE(9999)
      WRITE(6,*) '    - RESTART_MOD: Restart bxhght written!'
      ENDIF

      END SUBROUTINE WRITE_RESTART

!======================================================================

      SUBROUTINE GET_RESTART_DATA( OUTDATA )

      REAL*8, INTENT(OUT) :: OUTDATA(IIPAR,JJPAR,LLPAR)
      REAL*8 :: WORK(IGLOB*2,JGLOB/2,LLPAR)
      INTEGER :: IOS

      IF (MASTER) THEN
      READ(9999, IOSTAT=IOS) WORK
      CALL DYN_TRANS2MODEL_GRID(WORK,LLPAR)
      ENDIF

#ifdef MP
      CALL MPGF2P(WORK,IGLOB*2,JGLOB/2,OUTDATA,IIPAR,JJPAR,LLPAR)
#else
      OUTDATA = WORK
#endif

      IF (MASTER) THEN
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) ' '
         WRITE(6,*) "CANNOT READ : GET_RESTART_DATA"
         CALL FLUSH(6)
         STOP
      ENDIF
      ENDIF

      END SUBROUTINE GET_RESTART_DATA

!======================================================================

      SUBROUTINE PUT_RESTART_DATA( INDATA )

      REAL*8, INTENT(IN) :: INDATA(IIPAR,JJPAR,LLPAR)
      REAL*8 :: WORK(IGLOB*2,JGLOB/2,LLPAR)
      INTEGER :: IOS

#ifdef MP
      CALL MPGP2F(INDATA,IIPAR,JJPAR,WORK,IGLOB*2,JGLOB/2,LLPAR)
#else
      WORK = INDATA
#endif

      IF (MASTER) THEN
      CALL DYN_TRANS2OUTPUT_GRID(WORK,LLPAR)
      WRITE(9999, IOSTAT=IOS) WORK
      ENDIF

      IF (MASTER) THEN
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) ' '
         WRITE(6,*) "CANNOT WRITE : PUT_RESTART_DATA"
         CALL FLUSH(6)
         STOP
      ENDIF
      ENDIF

      END SUBROUTINE PUT_RESTART_DATA

!======================================================================

      END MODULE RESTART_MOD
