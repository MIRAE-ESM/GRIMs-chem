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

      USE COMODE_MOD,   ONLY : CSPEC
      USE DAO_MOD,      ONLY : AIRDEN, BXHEIGHT, RH
      USE LOGICAL_MOD,  ONLY : LCHEM, LEMIS, LMEGAN, LSVCSPEC
      USE MEGAN_MOD,    ONLY : DAY_DIM, NUM_DAYS
      USE MEGAN_MOD,    ONLY : T_DAY, T_15
      USE MEGAN_MOD,    ONLY : PARDR_DAY, PARDR_15
      USE MEGAN_MOD,    ONLY : PARDF_DAY, PARDF_15
      USE TRACER_MOD,   ONLY : STT, N_TRACERS
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM

      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD

      IMPLICIT NONE

      CONTAINS

!======================================================================

      SUBROUTINE READ_RESTART

      USE FILE_MOD, ONLY : FILE_EXISTS

      INTEGER :: I, J, L, N
      INTEGER :: NI, NJ, NL, NT
      REAL*8, ALLOCATABLE :: P_ARRAY(:,:,:)
      REAL*8, ALLOCATABLE :: F_ARRAY(:,:,:)

      ! Check chem initial condition file
      IF ( .not. FILE_EXISTS( 'chemi' ) ) THEN

         ! Set initial condition to be 1d-30
         STT = 1d-30

         ! Message
         IF (LMASTER) THEN
         WRITE(6,*) '    - RESTART_MOD: Zero concentration start!'
         ENDIF

         RETURN

      ENDIF

      !------------!
      ! Read chemi !
      !------------!

      IF (LMASTER) THEN
      ! Open file for read
      OPEN(9999, FILE='chemi', FORM='UNFORMATTED', STATUS='OLD')

      ! Read header
      READ(9999) NI, NJ, NL, NT

      ! Dimension check
      IF ( NI .NE. IGLOB .OR. NJ .NE. JGLOB .OR.
     &     NL .NE. LLPAR .OR. NT .NE. N_TRACERS ) THEN
         WRITE(6,*) '    - RESTART_MOD: Dimension of chemi is bad!'
         STOP
      ENDIF
      ENDIF

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, LLPAR ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, LLPAR ) )

      ! Read STT
      DO N = 1, N_TRACERS

         IF (LMASTER) THEN
         READ(9999) F_ARRAY
         ENDIF

         CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                      IGLOB*2, JGLOB/2, LLPAR, LLPAR, 0, 0, 0 )

         STT(:,:,:,N) = P_ARRAY

      ENDDO

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      IF (LMASTER) THEN
      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart chemi is read!'
      ENDIF

      !--------------!
      ! Read bxhghti !
      !--------------!

      ! Check bxhght initial condition file
      IF ( FILE_EXISTS( 'bxhghti' ) ) THEN

      IF (LMASTER) THEN
      ! Open file for read
      OPEN(9999, FILE='bxhghti', FORM='UNFORMATTED', STATUS='OLD')

      ! Read header
      READ(9999) NI, NJ, NL

      ! Dimension check
      IF ( NI .NE. IGLOB .OR. NJ .NE. JGLOB .OR.
     &     NL .NE. LLPAR ) THEN
         WRITE(6,*) '    - RESTART_MOD: Dimension of bxhghti is bad!'
         STOP
      ENDIF
      ENDIF

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, LLPAR ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, LLPAR ) )

      ! Read AIRDEN
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, LLPAR, LLPAR, 0, 0, 0 )

      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO L = 1, LLPAR
         AIRDEN(L,I,J) = P_ARRAY(I,J,L)
      ENDDO
      ENDDO
      ENDDO

      ! Read BXHEIGHT
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, LLPAR, LLPAR, 0, 0, 0 )

      BXHEIGHT = P_ARRAY

      ! Read RH
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, LLPAR, LLPAR, 0, 0, 0 )

      RH = P_ARRAY

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      IF (LMASTER) THEN
      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart bxhghti is read!'
      ENDIF

      ENDIF

      !-------------!
      ! Read megani !
      !-------------!

      ! Check megan initial condition file
      IF ( LMEGAN .and. FILE_EXISTS( 'megani' ) ) THEN

      IF (LMASTER) THEN
      ! Open file for read
      OPEN(9999, FILE='megani', FORM='UNFORMATTED', STATUS='OLD')

      ! Read header
      READ(9999) NI, NJ, NL, NT

      ! Dimension check
      IF ( NI .NE. IGLOB .OR. NJ .NE. JGLOB .OR.
     &     NL .NE. DAY_DIM .OR. NT .NE. NUM_DAYS ) THEN
         WRITE(6,*) '    - RESTART_MOD: Dimension of megani is bad!'
         STOP
      ENDIF
      ENDIF

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, DAY_DIM ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, DAY_DIM ) )

      ! Read T_DAY
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, DAY_DIM, DAY_DIM, 0, 0, 0 )

      T_DAY = P_ARRAY

      ! Read PARDR_DAY
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, DAY_DIM, DAY_DIM, 0, 0, 0 )

      PARDR_DAY = P_ARRAY

      ! Read PARDF_DAY
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, DAY_DIM, DAY_DIM, 0, 0, 0 )

      PARDF_DAY = P_ARRAY

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, NUM_DAYS ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, NUM_DAYS ) )

      ! Read T_15
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, NUM_DAYS, NUM_DAYS, 0, 0, 0 )

      T_15 = P_ARRAY

      ! Read PARDR_15
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, NUM_DAYS, NUM_DAYS, 0, 0, 0 )

      PARDR_15 = P_ARRAY

      ! Read PARDF_15
      IF (LMASTER) THEN
      READ(9999) F_ARRAY
      ENDIF

      CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                   IGLOB*2, JGLOB/2, NUM_DAYS, NUM_DAYS, 0, 0, 0 )

      PARDF_15 = P_ARRAY

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      IF (LMASTER) THEN
      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart megani is read!'
      ENDIF

      ENDIF

      !-------------!
      ! Read cspeci !
      !-------------!

      ! Check cspec initial condition file
      IF ( ( LEMIS .or. LCHEM ) .and. ITS_A_FULLCHEM_SIM() .and.
     &     FILE_EXISTS( 'cspeci' ) ) THEN

      IF (LMASTER) THEN
      ! Open file for read
      OPEN(9999, FILE='cspeci', FORM='UNFORMATTED', STATUS='OLD')

      ! Read header
      READ(9999) NI, NJ, NL, NT

      ! Dimension check
      IF ( NI .NE. IGLOB .OR. NJ .NE. JGLOB .OR.
     &     NL .NE. LLTROP ) THEN
         WRITE(6,*) '    - RESTART_MOD: Dimension of cspeci is bad!'
         STOP
      ENDIF
      ENDIF

#ifdef MP
#ifndef RMP
      CALL mpbcasti(NT,1)
#else
      CALL rmpbcasti(NT,1)
#endif
#endif

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, LLTROP ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, LLTROP ) )

      ! Read CSPEC
      DO N = 1, NT

         IF (LMASTER) THEN
         READ(9999) F_ARRAY
         ENDIF

         CALL TRANS_REGRID( F_ARRAY, P_ARRAY,
     &                      IGLOB*2, JGLOB/2, LLTROP, LLTROP, 0, 0, 0 )

         CSPEC(:,:,:,N) = P_ARRAY

      ENDDO

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      IF (LMASTER) THEN
      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart cspeci is read!'
      ENDIF

      LSVCSPEC = .TRUE.

      ENDIF

      END SUBROUTINE READ_RESTART

!======================================================================

      SUBROUTINE WRITE_RESTART

      INTEGER :: I, J, L, N
      REAL*8, ALLOCATABLE :: P_ARRAY(:,:,:)
      REAL*8, ALLOCATABLE :: F_ARRAY(:,:,:)

      !-------------!
      ! Write chemi !
      !-------------!

      IF (LMASTER) THEN
      ! Open file for write
      OPEN(9999, FILE='chemi', FORM='UNFORMATTED', STATUS='REPLACE')

      ! Write header
      WRITE(9999) IGLOB, JGLOB, LLPAR, N_TRACERS
      ENDIF

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, LLPAR ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, LLPAR ) )

      ! Write STT
      DO N = 1, N_TRACERS

         P_ARRAY = STT(:,:,:,N)

         CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, LLPAR )

         IF (LMASTER) THEN
         WRITE(9999) F_ARRAY
         ENDIF

      ENDDO

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      IF (LMASTER) THEN
      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart chemi is written!'
      ENDIF

      !---------------!
      ! Write bxhghti !
      !---------------!

      IF (LMASTER) THEN
      ! Open file for write
      OPEN(9999, FILE='bxhghti', FORM='UNFORMATTED', STATUS='REPLACE')

      ! Write header
      WRITE(9999) IGLOB, JGLOB, LLPAR
      ENDIF

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, LLPAR ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, LLPAR ) )

      ! Write AIRDEN
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         P_ARRAY(I,J,L) = AIRDEN(L,I,J)
      ENDDO
      ENDDO
      ENDDO

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, LLPAR )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Write BXHEIGHT
      P_ARRAY = BXHEIGHT

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, LLPAR )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Write RH
      P_ARRAY = RH

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, LLPAR )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      IF (LMASTER) THEN
      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart bxhghti is written!'
      ENDIF

      !--------------!
      ! Write megani !
      !--------------!

      IF ( LMEGAN ) THEN

      IF (LMASTER) THEN
      ! Open file for write
      OPEN(9999, FILE='megani', FORM='UNFORMATTED', STATUS='REPLACE')

      ! Write header
      WRITE(9999) IGLOB, JGLOB, DAY_DIM, NUM_DAYS
      ENDIF

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, DAY_DIM ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, DAY_DIM ) )

      ! Write T_DAY
      P_ARRAY = T_DAY

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, DAY_DIM )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Write PARDR_DAY
      P_ARRAY = PARDR_DAY

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, DAY_DIM )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Write PARDF_DAY
      P_ARRAY = PARDF_DAY

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, DAY_DIM )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, NUM_DAYS ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, NUM_DAYS ) )

      ! Write T_15
      P_ARRAY = T_15

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, NUM_DAYS )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Write PARDR_15
      P_ARRAY = PARDR_15

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, NUM_DAYS )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Write PARDF_15
      P_ARRAY = PARDF_15

      CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, NUM_DAYS )

      IF (LMASTER) THEN
      WRITE(9999) F_ARRAY
      ENDIF

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      IF (LMASTER) THEN
      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart megani is written!'
      ENDIF

      ENDIF

      !--------------!
      ! Write cspeci !
      !--------------!

      IF ( ( LEMIS .or. LCHEM ) .and. ITS_A_FULLCHEM_SIM() ) THEN

      IF (LMASTER) THEN
      ! Open file for write
      OPEN(9999, FILE='cspeci', FORM='UNFORMATTED', STATUS='REPLACE')

      ! Write header
      WRITE(9999) IGLOB, JGLOB, LLTROP, NTSPECGAS
      ENDIF

      ! Allocate
      ALLOCATE( P_ARRAY( IIPAR, JJPAR, LLTROP ) )
      ALLOCATE( F_ARRAY( IGLOB, JGLOB, LLTROP ) )

      ! Write CSPEC
      DO N = 1, NTSPECGAS

         P_ARRAY = CSPEC(:,:,:,N)

         CALL TRANS_OUTPUT( P_ARRAY, F_ARRAY, LLTROP )

         IF (LMASTER) THEN
         WRITE(9999) F_ARRAY
         ENDIF

      ENDDO

      ! Deallocate
      DEALLOCATE( P_ARRAY )
      DEALLOCATE( F_ARRAY )

      IF (LMASTER) THEN
      ! Close file
      CLOSE(9999)

      ! Message
      WRITE(6,*) '    - RESTART_MOD: Restart cspeci is written!'
      ENDIF

      ENDIF

      END SUBROUTINE WRITE_RESTART

!======================================================================

      END MODULE RESTART_MOD
