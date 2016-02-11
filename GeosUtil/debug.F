!=====================================================================
!
!     SUBROUTINE DEBUG
!
!        This subroutine writes an array to the disk for debugging.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE DEBUG( ARRAY, NL, FNAME )

      USE CMN_SIZE_MOD

      INTEGER :: NL
      REAL*8  :: ARRAY(IIPAR,JJPAR,NL)
      REAL*8  :: F_ARRAY(IGLOB,JGLOB,NL)
      CHARACTER(LEN=255) :: FNAME

      ! Transform
      CALL TRANS_OUTPUT( ARRAY, F_ARRAY, NL )

      ! Write
      IF (LMASTER) THEN
      OPEN(101, FILE=TRIM(FNAME), FORM='UNFORMATTED', STATUS='REPLACE')
      WRITE(101) F_ARRAY
      CLOSE(101)
      ENDIF

      END SUBROUTINE DEBUG
