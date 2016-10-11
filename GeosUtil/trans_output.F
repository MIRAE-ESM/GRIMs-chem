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
      REAL*8  :: OUTDATA(LONF2F,LATG2F,NL)
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
#ifndef RMP
      CALL mpgp2f(TMP,LONF2S,LATG2S,OUTDATA,LONF2F,LATG2F,NL)
#else
      CALL rmpgp2f(TMP,LONF2S,LATG2S,OUTDATA,LONF2F,LATG2F,NL)
#endif
#else
      OUTDATA = TMP
#endif

      IF (LMASTER) THEN
      ! Transform to output grid
#ifndef RMP
      CALL dyn_trans2output_grid(OUTDATA,NL)
      CALL XSHIFT(OUTDATA,LONFD,LATGD,NL)
      CALL YREV(OUTDATA,LONFD,LATGD,NL)
#else
      CALL rmp_trans2output_grid(OUTDATA,NL)
#endif
      ENDIF

      END SUBROUTINE TRANS_OUTPUT
