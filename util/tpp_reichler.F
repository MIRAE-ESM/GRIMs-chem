#include "define.h"
!=====================================================================
!
!     SUBROUTINE TPP_REICHLER
!
!        This subroutine calculates the tropopause pressure
!        based on the algorithm of Reichler et al., 2003.
!
!                                            Written by Seungun Lee
!
!=====================================================================

      SUBROUTINE TPP_REICHLER( TEMP, TPP_PRESS )

      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE CMN_SIZE_MOD

      IMPLICIT NONE

      REAL*8, INTENT(IN)  :: TEMP(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: TPP_PRESS(IIPAR,JJPAR)
      REAL*8, PARAMETER :: LAPSE_RATE=-2d-3
      REAL*8, PARAMETER :: KAP=0.286d0
      REAL*8, PARAMETER :: R=287d0
      REAL*8, PARAMETER :: G=9.81d0
      REAL*8, PARAMETER :: FAKTOR=-1d0*G/R
      REAL*8, PARAMETER :: DELTAZ=2000d0
      REAL*8, PARAMETER :: PLIMU=45000d0
      REAL*8, PARAMETER :: PLIML=7500d0
      REAL*8 :: P, PD, PU, PK
      REAL*8 :: P2, PD2, PU2, PK2
      REAL*8 :: DT, DPK, DTDPK, DTDP, DTDZ
      REAL*8 :: DT2, DPK2, DTDPK2, DTDP2, DTDZ2
      REAL*8 :: DTDZ0, PK0
      REAL*8 :: AG, BG, PTPH, P2KM, ASUM, AQUER
      INTEGER :: I, J, L, L2, ICOUNT

      TPP_PRESS = 0d0

      DO J = 1, LATLOOP
      DO I = 1, LONLOOP

         DO L = 2, LLPAR-1

            P  = GET_PCENTER(I,J,L) * 1d2
            PD = GET_PCENTER(I,J,L-1) * 1d2
            PU = GET_PCENTER(I,J,L+1) * 1d2
            PK = P ** KAP

            DT    = TEMP(I,J,L+1) - TEMP(I,J,L-1)
            DPK   = PU ** KAP - PD ** KAP
            DTDPK = DT / DPK
            DTDP  = DTDPK * KAP * P ** (KAP-1d0)
            DTDZ  = FAKTOR * DTDP * P / TEMP(I,J,L)

            IF ( L .EQ. 2 ) GOTO 100
            IF ( DTDZ .LE. LAPSE_RATE ) GOTO 100
            IF ( P .GT. PLIMU ) GOTO 100
            IF ( DTDZ0 .LT. LAPSE_RATE ) THEN
               AG = ( DTDZ - DTDZ0 ) / ( PK - PK0 )
               BG = DTDZ0 - AG * PK0
               PTPH = EXP( ALOG( ( LAPSE_RATE - BG ) / AG ) / KAP )
            ELSE
               PTPH = P
            ENDIF
            IF ( PTPH .LT. PLIML ) GOTO 300
            IF ( PTPH .GT. PLIMU ) GOTO 100

            P2KM = PTPH + DELTAZ * P / TEMP(I,J,L) * FAKTOR
            ASUM = 0d0
            ICOUNT = 0

            DO L2 = L, LLPAR-1

               P2  = GET_PCENTER(I,J,L2) * 1d2
               PD2 = GET_PCENTER(I,J,L2-1) * 1d2
               PU2 = GET_PCENTER(I,J,L2+1) * 1d2
               PK2 = P2 ** KAP
               IF ( P2 .GT. PTPH ) GOTO 200
               IF ( P2 .LT. P2KM ) GOTO 300

               DT2    = TEMP(I,J,L2+1) - TEMP(I,J,L2-1)
               DPK2   = PU2 ** KAP - PD2 ** KAP
               DTDPK2 = DT2 / DPK2
               DTDP2 = DTDPK2 * KAP * P2 ** (KAP-1d0)
               DTDZ2 = FAKTOR * DTDP2 * P2 / TEMP(I,J,L2)

               ASUM = ASUM + DTDZ2
               ICOUNT = ICOUNT + 1
               AQUER = ASUM / FLOAT(ICOUNT)

               IF ( AQUER .LE. LAPSE_RATE ) GOTO 100

 200           CONTINUE

            ENDDO

 100        PK0 = PK
            DTDZ0 = DTDZ

         ENDDO

         PRINT*, "TPP_REICHLER: Couldn't find the tropopause."
#ifdef MP
         CALL MPABORT
#else
         STOP
#endif

 300     TPP_PRESS(I,J) = MAX( PTPH, PLIML )

      ENDDO
      ENDDO

      END SUBROUTINE TPP_REICHLER
