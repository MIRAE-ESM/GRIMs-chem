!=====================================================================
!
!     SUBROUTINE INTEGRATE_CHEM
!
!        This subroutine is a driver subroutine for gmp_integrate.F90.
!
!                                            written by Seungun Lee
!
!=====================================================================

      SUBROUTINE INTEGRATE_CHEM

      USE BPCH2_MOD,      ONLY : OPEN_BPCH2_FOR_WRITE
      USE CHEMISTRY_MOD,  ONLY : DO_CHEMISTRY
      USE COMFVER,        ONLY : STEPONE
      USE CONVECTION_MOD, ONLY : DO_CONVECTION
      USE DAO_MOD,        ONLY : AD, AIRQNT, CONVERT_UNITS, COSSZA
      USE DAO_MOD,        ONLY : MAKE_RH, PS1, SUNCOS, SUNCOS_MID
      USE DAO_MOD,        ONLY : T, TROPP
      USE DIAG41_MOD,     ONLY : DIAG41, ND41
      USE DIAG49_MOD,     ONLY : DIAG49, ITS_TIME_FOR_DIAG49
      USE DIAG_AOD_MOD,   ONLY : WRITE_DIAG_AOD
      USE DRYDEP_MOD,     ONLY : DO_DRYDEP
      USE EMISSIONS_MOD,  ONLY : DO_EMISSIONS
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE FILE_MOD,       ONLY : IU_BPCH, IU_DEBUG, IU_ND48, IU_SMV2LOG
      USE LINOZ_MOD,      ONLY : DO_LINOZ
      USE LOGICAL_MOD,    ONLY : LCHEM, LCONV, LDRYD, LEMIS, LLINOZ
      USE LOGICAL_MOD,    ONLY : LPRT, LTURB, LVARTROP, LWETD
      USE PBL_MIX_MOD,    ONLY : DO_PBL_MIX
      USE PRESSURE_MOD,   ONLY : SET_FLOATING_PRESSURE
      USE TIME_MOD,       ONLY : GET_NYMD, GET_NHMS, GET_TAU
      USE TIME_MOD,       ONLY : GET_YEAR, GET_SEASON, GET_MONTH
      USE TIME_MOD,       ONLY : GET_HOUR, GET_MINUTE
      USE TIME_MOD,       ONLY : GET_ELAPSED_SEC, GET_DAY_OF_YEAR
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_DAY
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_BPCH
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_CHEM
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_CONV
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_DIAG
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_DYN
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_EMIS
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_UNIT
      USE TIME_MOD,       ONLY : SET_CT_CHEM, SET_CT_CONV, SET_CT_DIAG
      USE TIME_MOD,       ONLY : SET_CT_DYN, SET_CT_EMIS
      USE TIME_MOD,       ONLY : SET_DIAGb, SET_DIAGe, TIMESTAMP_DIAG
      USE TIME_MOD,       ONLY : SET_CURRENT_TIME, PRINT_CURRENT_TIME
      USE TIME_MOD,       ONLY : SET_ELAPSED_MIN
      USE TIME_MOD,       ONLY : EXPAND_DATE
      USE TRACER_MOD,     ONLY : CHECK_STT, N_TRACERS
      USE TRACER_MOD,     ONLY : STT1, STT2, STT, TCVV
      USE TRACER_MOD,     ONLY : ITS_A_H2HD_SIM
      USE TROPOPAUSE_MOD, ONLY : CHECK_VAR_TROP, DIAG_TROPOPAUSE
      USE WETSCAV_MOD,    ONLY : DO_WETDEP, INIT_WETSCAV

      USE CMN_SIZE_MOD    ! Size parameters   
      USE COMODE_LOOP_MOD ! SMVGEAR common blocks
      USE CMN_DIAG_MOD    ! Diagnostic switches, NJDAY
      USE CMN_GCTM_MOD    ! Physical constants

      IMPLICIT NONE

      INTEGER :: DAY_OF_YEAR, ELAPSED_SEC, MONTH, NHMS, NYMD, HOUR
      INTEGER :: MINUTE, YEAR, SEASON, ELAPSED_TODAY
      REAL*8  :: TAU
      CHARACTER(LEN=255) :: FILENAME
      LOGICAL, SAVE :: FIRST = .TRUE.

      IF (STEPONE) THEN

         ! Synchronize with shour
         CALL SET_ELAPSED_MIN

      ELSE

         ! 3-D Cloud Voluem Fraction
         CALL MAKE_CLDF

         ! Tropopause [hPa]
         CALL TPP_REICHLER( T, TROPP )
         TROPP = TROPP * 1d-2

         IF (FIRST) THEN

            ! Call AIRQNT to compute air mass quantities from PS1
            CALL SET_FLOATING_PRESSURE( PS1 )
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a SET_FLT_PRS' )

            CALL AIRQNT
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a AIRQNT' )

            ! Read land types and fractions from "vegtype.global"
            CALL RDLAND
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a RDLAND' )

            ! Initialize PBL quantities but do not do mixing
            CALL DO_PBL_MIX( .FALSE. )
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:1' )

            FIRST = .FALSE.

         ENDIF

         ! Compute & print time quantities at start of dyn step
         CALL SET_CURRENT_TIME
         CALL PRINT_CURRENT_TIME

         ! Set time variables for dynamic loop
         DAY_OF_YEAR   = GET_DAY_OF_YEAR()
         ELAPSED_SEC   = GET_ELAPSED_SEC()
         MONTH         = GET_MONTH()
         NHMS          = GET_NHMS()
         NYMD          = GET_NYMD()
         HOUR          = GET_HOUR()
         MINUTE        = GET_MINUTE()
         TAU           = GET_TAU()
         YEAR          = GET_YEAR()
         SEASON        = GET_SEASON()
         ELAPSED_TODAY = ( HOUR * 3600 ) + ( MINUTE * 60 )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a SET_CURRENT_TIME' )

         !==============================================================
         !   ***** W R I T E   D I A G N O S T I C   F I L E S *****
         !==============================================================
         IF ( ITS_TIME_FOR_BPCH() ) THEN

            ! Set time at end of diagnostic timestep
            CALL SET_DIAGe( TAU )

            ! Write bpch file
            CALL DIAG3

            ! Write AOD, SSA, ASM
            CALL WRITE_DIAG_AOD

            ! Flush file units
            CALL FLUSH( IU_ND48    )          
            CALL FLUSH( IU_BPCH    )          
            CALL FLUSH( IU_SMV2LOG )          
            CALL FLUSH( IU_DEBUG   )          

            ! Open a new diagnostic file after saving bpch file
            FILENAME = 'ctm.bpch.YYYYMMDDhh'

            ! Expand YYYYMMDD tokens in the bpch file name
            CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

            ! Open bpch file
            CALL OPEN_BPCH2_FOR_WRITE( IU_BPCH, FILENAME )

            ! Set time at beginning of next diagnostic timestep
            CALL SET_DIAGb( TAU )

            !===========================================================
            !        ***** Z E R O   D I A G N O S T I C S *****
            !===========================================================
            CALL INITIALIZE( 2 ) ! Zero arrays
            CALL INITIALIZE( 3 ) ! Zero counters
         ENDIF

         !==============================================================
         !              ***** D A I L Y   D A T A *****
         !==============================================================
         IF ( ITS_A_NEW_DAY() ) THEN

            ! Read leaf-area index (needed for drydep)
            ! Now include the year (mpb,11/19/09)
            CALL RDLAI( DAY_OF_YEAR, MONTH, YEAR )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG ( '### MAIN: a DAILY DATA' )
         ENDIF

         ! Case of variable tropopause:
         ! Check LLTROP and set LMIN, LMAX, and LPAUSE
         ! since this is not done with READ_TROPOPAUSE anymore.
         ! (Need to double-check that LMIN, Lmax are not used before-phs) 
         IF ( LVARTROP ) CALL CHECK_VAR_TROP

         ! Compute airmass quantities at each grid box 
         CALL SET_FLOATING_PRESSURE( PS1 )
         CALL AIRQNT

         ! Compute the cosine of the solar zenith angle array
         ! SUNCOS     = cos(SZA) at the current time
         ! SUNCOS_MID = cos(SZA) at the midpt of the chemistry timestep
!------------------------------------------------------------------------------
! Prior to 10/7/11:
!         CALL COSSZA( DAY_OF_YEAR, SUNCOS )
!------------------------------------------------------------------------------
         CALL COSSZA( DAY_OF_YEAR, SUNCOS, SUNCOS_MID )

         ! Compute tropopause height for ND55 diagnostic
         IF ( ND55 > 0 ) CALL DIAG_TROPOPAUSE

         ! Update dynamic timestep
         CALL SET_CT_DYN( INCREMENT=.TRUE. )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INTERP, etc' )

         ! Get averaging intervals for local-time diagnostics
         ! (NOTE: maybe improve this later on)
         ! Placed after interpolation to get correct value of TROPP. 
         ! (ccc, 12/9/08)
         CALL DIAG_2PM

         !==============================================================
         !     ***** S T R A T O S P H E R I C   O Z O N E S *****
         !==============================================================
         IF ( LLINOZ ) CALL DO_LINOZ

         !==============================================================
         !              ***** T R A N S P O R T *****
         !==============================================================
         IF ( ITS_TIME_FOR_DYN() ) THEN

            ! Transport of chemical tracers is moved to
            ! gmp_integrate.F90, from this routine.

            ! Get relative humidity (after recomputing pressures)
            ! NOTE: for GEOS-5 we'll read this from disk instead
            CALL MAKE_RH

            ! Initialize wet scavenging and wetdep fields after
            ! the airmass quantities are reset after transport
            IF ( LCONV .or. LWETD ) CALL INIT_WETSCAV
         ENDIF

         !-------------------------------
         ! Test for convection timestep
         !-------------------------------
         IF ( ITS_TIME_FOR_CONV() ) THEN

            ! Increment the convection timestep
            CALL SET_CT_CONV( INCREMENT=.TRUE. )

            !===========================================================
            !      ***** M I X E D   L A Y E R   M I X I N G *****
            !===========================================================
            CALL DO_PBL_MIX( LTURB )
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:1' )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:2' )

            !===========================================================
            !        ***** C L O U D   C O N V E C T I O N *****
            !===========================================================
            IF ( LCONV ) THEN

               CALL DO_CONVECTION

               !### Debug
               IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVECTION' )
            ENDIF
         ENDIF


         !==============================================================
         !    ***** U N I T   C O N V E R S I O N  ( v/v -> kg ) *****
         !==============================================================
         IF ( ITS_TIME_FOR_UNIT() ) THEN
            CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:1' )
         ENDIF

         !-------------------------------
         ! Test for emission timestep
         !-------------------------------
         IF ( ITS_TIME_FOR_EMIS() ) THEN
            
            ! Increment emission counter
            CALL SET_CT_EMIS( INCREMENT=.TRUE. )
            
            !========================================================
            !         ***** D R Y   D E P O S I T I O N *****
            !========================================================
            IF ( LDRYD .and. ( .not. ITS_A_H2HD_SIM() ) ) CALL DO_DRYDEP
            
            !========================================================
            !             ***** E M I S S I O N S *****
            !========================================================
            IF ( LEMIS ) CALL DO_EMISSIONS
         ENDIF

         
         !===========================================================
         !               ***** C H E M I S T R Y *****
         !===========================================================    
         
         ! Every chemistry timestep...
         IF ( ITS_TIME_FOR_CHEM() ) THEN
            
            ! Increment chemistry timestep counter
            CALL SET_CT_CHEM( INCREMENT=.TRUE. )
            
            ! Call the appropriate chemistry routine
            CALL DO_CHEMISTRY
         
         ENDIF
         
         !==============================================================
         ! ***** W E T   D E P O S I T I O N  (rainout + washout) *****
         !==============================================================
         IF ( LWETD .and. ITS_TIME_FOR_DYN() ) THEN
            
            CALL DO_WETDEP
            
         ENDIF

         !==============================================================
         !   ***** I N C R E M E N T   E L A P S E D   T I M E *****
         !============================================================== 
         ! Moved before diagnostics to count the last timestep as done.
         ! Need to save timestamps for filenames.
         ! (ccc, 5/13/09)

         CALL TIMESTAMP_DIAG
         CALL SET_ELAPSED_MIN
         CALL SET_CURRENT_TIME
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after SET_ELAPSED_MIN' )

         !==============================================================
         !       ***** A R C H I V E   D I A G N O S T I C S *****
         !==============================================================
         IF ( ITS_TIME_FOR_DIAG() ) THEN

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: b DIAGNOSTICS' )

            ! Accumulate several diagnostic quantities
            CALL DIAG1
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG1' )

            ! ND41: save PBL height in 1200-1600 LT (amf)
            ! (for comparison w/ Holzworth, 1967)
            IF ( ND41 > 0 ) CALL DIAG41
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG41' )

            ! Increment diagnostic timestep counter. (ccc, 5/13/09)
            CALL SET_CT_DIAG( INCREMENT=.TRUE. )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a DIAGNOSTICS' )
         ENDIF

         !==============================================================
         !   ***** T I M E S E R I E S   D I A G N O S T I C S  *****
         !
         ! NOTE: Since we are saving soluble tracers, we must move
         !       the ND40, ND49, and ND52 timeseries diagnostics
         !       to after the call to DO_WETDEP (bmy, 4/22/04)
         !============================================================== 

         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: before TIMESERIES' )

         ! 3-D timeseries
         IF ( ITS_TIME_FOR_DIAG49() ) CALL DIAG49
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG49' )

         !### After diagnostics
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after TIMESERIES' )

         !==============================================================
         !   ***** U N I T   C O N V E R S I O N  ( kg -> v/v ) *****
         !==============================================================
         IF ( ITS_TIME_FOR_UNIT() ) THEN
            CALL CONVERT_UNITS( 1,  N_TRACERS, TCVV, AD, STT )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:2' )
         ENDIF

         !==============================================================
         !  ***** E N D   O F   D Y N A M I C   T I M E S T E P *****
         !==============================================================
         
         ! Check for NaN, Negatives, Infinities in STT each time diag are
         ! saved. (ccc, 5/13/09)
         IF ( ITS_TIME_FOR_DIAG() ) THEN
            CALL CHECK_STT( 'End of Dynamic Loop' )
         ENDIF

         ! Update STT
         STT1 = STT2
         STT2 = STT

      ENDIF

      END SUBROUTINE INTEGRATE_CHEM
