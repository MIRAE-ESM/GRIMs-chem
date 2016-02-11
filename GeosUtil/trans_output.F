#include "define.h"
!=====================================================================
!
!     SUBROUTINE TRANS_OUTPUT
!
!        This subroutine transforms an array from model grid to
!        output grid.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE TRANS_OUTPUT( INDATA, OUTDATA, NL )

      USE CMN_SIZE_MOD
      USE module_trans
      USE paramodel
#ifdef MP
      USE commpi
#endif

      IMPLICIT NONE

      INTEGER :: NL
      INTEGER :: I, J, L
      REAL*8  :: INDATA(IIPAR,JJPAR,NL)
      REAL*8  :: OUTDATA(lonf2_,latg2_,NL)
      REAL*8  :: TMP(LONF2S,LATG2S,NL)

      ! Copy
      DO L = 1, NL
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         TMP(I,J,L) = INDATA(I,J,L)
      ENDDO
      ENDDO
      ENDDO

#ifdef MP
      ! MPI partial grid to full grid
      CALL mpgp2f(TMP,LONF2S,LATG2S,OUTDATA,lonf2_,latg2_,NL)
#else
      OUTDATA = TMP
#endif

      IF (LMASTER) THEN
      ! Transform to output grid
      CALL dyn_trans2output_grid(OUTDATA,NL)
      CALL XSHIFT(OUTDATA,lonf_,latg_,NL)
      CALL YREV(OUTDATA,lonf_,latg_,NL)
      ENDIF

      END SUBROUTINE TRANS_OUTPUT
