#include "define.h"
!=====================================================================
!
!     MODULE RESTART
!
!        This module makes restart file for chemistry module.
!
!                                                       Seungun Lee
!
!=====================================================================

      MODULE RESTART_MOD

      USE DAO_MOD,      ONLY : AIRDEN, BXHEIGHT, RH
      USE TRACER_MOD,   ONLY : STT, N_TRACERS

      USE CMN_SIZE_MOD

      IMPLICIT NONE

      CONTAINS

!======================================================================

      SUBROUTINE READ_RESTART

      INTEGER :: I, J, L, N, IOS
      INTEGER :: NI, NJ, NL, NT
      REAL*8  :: P_AIRDEN(IIPAR,JJPAR,LLPAR)
      REAL*8  :: F_STT(IGLOB,JGLOB,LLPAR,N_TRACERS)
      REAL*8  :: F_AIRDEN(IGLOB,JGLOB,LLPAR)
      REAL*8  :: F_BXHEIGHT(IGLOB,JGLOB,LLPAR)
      REAL*8  :: F_RH(IGLOB,JGLOB,LLPAR)

      IF (LMASTER) THEN
      ! Open file for read
      OPEN(9999, FILE='chemi', FORM='UNFORMATTED', STATUS='OLD',
     &     IOSTAT=IOS)
      ENDIF

#ifdef MP
#ifndef RMP
      CALL mpbcasti(IOS,1)
#else
      CALL rmpbcasti(IOS,1)
#endif
#endif

      ! IOS check
      IF (IOS .NE. 0) THEN

         ! Set initial condition to be 1d-30
         STT = 1d-30

         ! Message
         IF (LMASTER) THEN
         WRITE(6,*) '    - RESTART_MOD: Zero concentration start!'
         ENDIF

         RETURN

      ENDIF

      IF (LMASTER) THEN
      ! Read
      READ(9999) NI, NJ, NL, NT

      ! Dimension check
      IF ( NI .NE. IGLOB .OR. NJ .NE. JGLOB .OR.
     &     NL .NE. LLPAR .OR. NT .NE. N_TRACERS ) THEN
         WRITE(6,*) '    - RESTART_MOD: Dimension of chemi is bad!'
         STOP
      ENDIF

      ! Read
      DO N = 1, N_TRACERS
         READ(9999) F_STT(:,:,:,N)
      ENDDO

      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart chemi is read!'
      ENDIF

      ! Transform
      DO N = 1, N_TRACERS
         CALL TRANS_REGRID( F_STT(:,:,:,N), STT(:,:,:,N),
     &                      IGLOB*2, JGLOB/2, LLPAR, LLPAR, 0, 0, 0 )
      ENDDO

      IF (LMASTER) THEN
      ! Open file for read
      OPEN(9999, FILE='bxhght', FORM='UNFORMATTED', STATUS='OLD',
     &     IOSTAT=IOS)
      ENDIF

#ifdef MP
#ifndef RMP
      CALL mpbcasti(IOS,1)
#else
      CALL rmpbcasti(IOS,1)
#endif
#endif

      ! IOS check
      IF (IOS .NE. 0) THEN

         ! Message
         IF (LMASTER) THEN
         WRITE(6,*) '    - RESTART_MOD: bxhght not found!'
         ENDIF

         RETURN

      ENDIF

      IF (LMASTER) THEN
      ! Read
      READ(9999) NI, NJ, NL

      ! Dimension check
      IF ( NI .NE. IGLOB .OR. NJ .NE. JGLOB .OR.
     &     NL .NE. LLPAR ) THEN
         WRITE(6,*) '    - RESTART_MOD: Dimension of bxhght is bad!'
         STOP
      ENDIF

      ! Read
      READ(9999) F_AIRDEN
      READ(9999) F_BXHEIGHT
      READ(9999) F_RH

      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart bxhght is read!'
      ENDIF

      ! Transform
      CALL TRANS_REGRID( F_AIRDEN, P_AIRDEN,
     &                   IGLOB*2, JGLOB/2, LLPAR, LLPAR, 0, 0, 0 )
      CALL TRANS_REGRID( F_BXHEIGHT, BXHEIGHT,
     &                   IGLOB*2, JGLOB/2, LLPAR, LLPAR, 0, 0, 0 )
      CALL TRANS_REGRID( F_RH, RH,
     &                   IGLOB*2, JGLOB/2, LLPAR, LLPAR, 0, 0, 0 )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO L = 1, LLPAR
         AIRDEN(L,I,J) = P_AIRDEN(I,J,L)
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE READ_RESTART

!======================================================================

      SUBROUTINE WRITE_RESTART

      INTEGER :: I, J, L, N
      REAL*8  :: P_AIRDEN(IIPAR,JJPAR,LLPAR)
      REAL*8  :: F_STT(IGLOB,JGLOB,LLPAR,N_TRACERS)
      REAL*8  :: F_AIRDEN(IGLOB,JGLOB,LLPAR)
      REAL*8  :: F_BXHEIGHT(IGLOB,JGLOB,LLPAR)
      REAL*8  :: F_RH(IGLOB,JGLOB,LLPAR)

      ! Transform
      DO N = 1, N_TRACERS
         CALL TRANS_OUTPUT( STT(:,:,:,N), F_STT(:,:,:,N), LLPAR )
      ENDDO

      IF (LMASTER) THEN
      ! Open file for write
      OPEN(9999, FILE='chemi', FORM='UNFORMATTED', STATUS='REPLACE')

      ! Write
      WRITE(9999) IGLOB, JGLOB, LLPAR, N_TRACERS
      DO N = 1, N_TRACERS
         WRITE(9999) F_STT(:,:,:,N)
      ENDDO

      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart chemi is written!'
      ENDIF

      ! Transform
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         P_AIRDEN(I,J,L) = AIRDEN(L,I,J)
      ENDDO
      ENDDO
      ENDDO
      CALL TRANS_OUTPUT( P_AIRDEN, F_AIRDEN, LLPAR )
      CALL TRANS_OUTPUT( BXHEIGHT, F_BXHEIGHT, LLPAR )
      CALL TRANS_OUTPUT( RH, F_RH, LLPAR )

      IF (LMASTER) THEN
      ! Open file for write
      OPEN(9999, FILE='bxhght', FORM='UNFORMATTED', STATUS='REPLACE')

      ! Write
      WRITE(9999) IGLOB, JGLOB, LLPAR
      WRITE(9999) F_AIRDEN
      WRITE(9999) F_BXHEIGHT
      WRITE(9999) F_RH

      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart bxhght is written!'
      ENDIF

      END SUBROUTINE WRITE_RESTART

!======================================================================

      END MODULE RESTART_MOD
