!=====================================================================
!
!     SUBROUTINE INIT_HEADER
!
!        This subroutine initializes header modules.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE INIT_HEADER

      USE CMN_DEP_MOD,     ONLY : INIT_CMN_DEP
      USE CMN_DIAG_MOD,    ONLY : INIT_CMN_DIAG
      USE CMN_FJ_MOD,      ONLY : INIT_CMN_FJ
      USE CMN_ISOP_MOD,    ONLY : INIT_CMN_ISOP
      USE CMN_MOD,         ONLY : INIT_CMN
      USE CMN_MONOT_MOD,   ONLY : INIT_CMN_MONOT
      USE CMN_NOX_MOD,     ONLY : INIT_CMN_NOX
      USE CMN_O3_MOD,      ONLY : INIT_CMN_O3
      USE CMN_SIZE_MOD,    ONLY : INIT_CMN_SIZE
      USE CMN_VEL_MOD,     ONLY : INIT_CMN_VEL
      USE COMODE_LOOP_MOD, ONLY : INIT_COMODE_LOOP
      USE JV_CMN_MOD,      ONLY : INIT_JV_CMN

      IMPLICIT NONE

      ! Initialize header modules
      CALL INIT_CMN_SIZE
      CALL INIT_CMN_DEP
      CALL INIT_CMN_DIAG
      CALL INIT_CMN_FJ
      CALL INIT_CMN_ISOP
      CALL INIT_CMN
      CALL INIT_CMN_MONOT
      CALL INIT_CMN_NOX
      CALL INIT_CMN_O3
      CALL INIT_CMN_VEL
      CALL INIT_COMODE_LOOP
      CALL INIT_JV_CMN

      END SUBROUTINE INIT_HEADER
