#include <define.h>
!=====================================================================
!
!     SUBROUTINE INTEGRATE_CHEM
!
!        This subroutine is a driver subroutine for gmp_integrate.F90.
!
!                                                       Seungun Lee
!
!=====================================================================

      SUBROUTINE INTEGRATE_CHEM

      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_WRITE
      USE CHEMISTRY_MOD,     ONLY : DO_CHEMISTRY, RECOMPUTE_OD
      USE CONVECTION_MOD,    ONLY : DO_CONVECTION
      USE DAO_MOD,           ONLY : AD, AIRQNT, CONVERT_UNITS, COSSZA
      USE DAO_MOD,           ONLY : MAKE_RH, PS1, PSC2, SUNCOS
      USE DAO_MOD,           ONLY : SUNCOS_MID, T, TROPP
      USE DIAG41_MOD,        ONLY : DIAG41, ND41
      USE DIAG49_MOD,        ONLY : DIAG49, ITS_TIME_FOR_DIAG49
      USE DIAG_AOD_MOD,      ONLY : WRITE_DIAG_AOD
      USE DRYDEP_MOD,        ONLY : DO_DRYDEP
      USE EMISSIONS_MOD,     ONLY : DO_EMISSIONS
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_BPCH, IU_DEBUG
      USE FILE_MOD,          ONLY : IU_ND48, IU_SMV2LOG
      USE LAI_MOD,           ONLY : RDISOLAI
      USE LIGHTNING_NOX_MOD, ONLY : LIGHTNING
      USE LOGICAL_MOD,       ONLY : LCHEM, LCONV, LDRYD, LEMIS, LNLPBL
      USE LOGICAL_MOD,       ONLY : LPRT, LTRAN, LTURB, LVARTROP, LWETD
      USE LOGICAL_MOD,       ONLY : LSVMET, LSVMETNEST
      USE LOGICAL_MOD,       ONLY : LLIGHTNOX, LMEGAN, LUPBD
      USE LOGICAL_MOD,       ONLY : LAERODIR, LAEROIND
      USE MEGAN_MOD,         ONLY : UPDATE_T_15_AVG
      USE MEGAN_MOD,         ONLY : UPDATE_T_DAY
      USE METEO_MOD,         ONLY : DO_METEO, INIT_A3, INIT_A6
      USE METEO_MOD,         ONLY : SAVE_A3, SAVE_A6, SAVE_I6
      USE METEO_MOD,         ONLY : SAVE_A3_NESTED, SAVE_A6_NESTED
      USE METEO_MOD,         ONLY : SAVE_I6_NESTED
      USE PBL_MIX_MOD,       ONLY : DO_PBL_MIX
      USE PRESSURE_MOD,      ONLY : SET_FLOATING_PRESSURE
      USE TIME_MOD,          ONLY : GET_NYMD, GET_NHMS, GET_TAU
      USE TIME_MOD,          ONLY : GET_YEAR, GET_SEASON, GET_MONTH
      USE TIME_MOD,          ONLY : GET_HOUR, GET_MINUTE
      USE TIME_MOD,          ONLY : GET_ELAPSED_SEC, GET_DAY_OF_YEAR
      USE TIME_MOD,          ONLY : ITS_A_NEW_DAY
      USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,          ONLY : ITS_A_NEW_SEASON
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_BPCH
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_CHEM
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_CONV
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_DIAG
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_DYN
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_EMIS
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_UNIT
      USE TIME_MOD,          ONLY : SET_CT_CHEM, SET_CT_CONV
      USE TIME_MOD,          ONLY : SET_CT_DIAG, SET_CT_DYN, SET_CT_EMIS
      USE TIME_MOD,          ONLY : SET_DIAGb, SET_DIAGe, TIMESTAMP_DIAG
      USE TIME_MOD,          ONLY : SET_CURRENT_TIME, PRINT_CURRENT_TIME
      USE TIME_MOD,          ONLY : SET_ELAPSED_MIN
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : CHECK_STT, N_TRACERS
      USE TRACER_MOD,        ONLY : STT, TCVV
#ifdef CHEM_NEWTP
      USE TRACER_MOD,        ONLY : STT_TP
#endif
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRANSPORT_MOD,     ONLY : DO_TRANSPORT
      USE TROPOPAUSE_MOD,    ONLY : CHECK_VAR_TROP, DIAG_TROPOPAUSE
      USE UPBDFLX_MOD,       ONLY : DO_UPBDFLX, UPBDFLX_NOY
      USE UVALBEDO_MOD,      ONLY : READ_UVALBEDO
      USE VDIFF_MOD,         ONLY : DO_PBL_MIX_2
      USE WETSCAV_MOD,       ONLY : DO_WETDEP, INIT_WETSCAV
      USE comfver,           ONLY : inistp
      USE tprof,             ONLY : tprof_s, tprof_e

      USE CMN_DIAG_MOD    ! Diagnostic switches, NJDAY
      USE CMN_GCTM_MOD    ! Physical constants

      IMPLICIT NONE

      INTEGER :: DAY_OF_YEAR, ELAPSED_SEC, MONTH, NHMS, NYMD, HOUR
      INTEGER :: MINUTE, YEAR, SEASON, ELAPSED_TODAY
      INTEGER :: ierr
      REAL*8  :: TAU
      CHARACTER(LEN=255) :: FILENAME
      LOGICAL, SAVE :: FIRST = .TRUE.

      IF (inistp) THEN
#ifndef CHEM_NEWTP
         ierr = tprof_s('    transport', 2)
         ierr = tprof_e('    transport', 2)
#endif
         IF ( LNLPBL ) THEN
            ierr = tprof_s('    drydep', 2)
            ierr = tprof_e('    drydep', 2)
            ierr = tprof_s('    emission', 2)
            ierr = tprof_e('    emission', 2)
         ENDIF
         ierr = tprof_s('    pblmix', 2)
         ierr = tprof_e('    pblmix', 2)
         ierr = tprof_s('    convection', 2)
         ierr = tprof_e('    convection', 2)
         IF ( .NOT. LNLPBL ) THEN
            ierr = tprof_s('    drydep', 2)
            ierr = tprof_e('    drydep', 2)
            ierr = tprof_s('    emission', 2)
            ierr = tprof_e('    emission', 2)
         ENDIF
         ierr = tprof_s('    chemistry', 2)
         ierr = tprof_e('    chemistry', 2)
         ierr = tprof_s('    wetdep', 2)
         ierr = tprof_e('    wetdep', 2)
         RETURN
      ENDIF

      IF (FIRST) THEN

         ! Read land types and fractions from "vegtype.global"
         CALL RDLAND
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a RDLAND' )

         ! Initialize PBL quantities but do not do mixing
         ! Add option for non-local PBL (Lin, 03/31/09) 
         IF ( .NOT. LNLPBL ) THEN
           CALL DO_PBL_MIX( .FALSE. )
           IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:1' )
         ELSE
           CALL DO_PBL_MIX_2( .FALSE. )
           IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a NLPBL 1' )
         ENDIF

         ! Save I6 fields
         IF ( MOD(GET_NHMS(),010000) .EQ. 0 ) THEN
            IF ( LSVMET     ) CALL SAVE_I6( FIRST )
            IF ( LSVMETNEST ) CALL SAVE_I6_NESTED( FIRST )
         ENDIF

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
      !   ***** U N I T   C O N V E R S I O N  ( v/v -> kg ) *****
      !==============================================================
      IF ( ITS_TIME_FOR_UNIT() ) THEN
         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:0' )
      ENDIF

      !==============================================================
      !   ***** W R I T E   D I A G N O S T I C   F I L E S *****
      !==============================================================
      IF ( ITS_TIME_FOR_BPCH() ) THEN

         ! Set time at end of diagnostic timestep
         CALL SET_DIAGe( TAU )

         ! Write bpch file
         CALL DIAG3

         ! Write AOD, SSA, ASM
         IF ( LAERODIR .or. LAEROIND ) THEN
            CALL WRITE_DIAG_AOD
         ENDIF

         ! Flush file units
         CALL FLUSH( IU_ND48    )          
         CALL FLUSH( IU_BPCH    )          
         CALL FLUSH( IU_SMV2LOG )          
         CALL FLUSH( IU_DEBUG   )          

         ! Open a new diagnostic file after saving bpch file
         FILENAME = 'ctm.bpch.YYYYMMDDhh'

         ! Expand YYYYMMDD tokens in the bpch file name
         CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

         ! Close bpch file
         CLOSE( IU_BPCH )

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

      ! Save GRIMs meteorology for offline GRIMs/GEOS-Chem
      IF ( LSVMET .or. LSVMETNEST ) THEN

         ! Save A3 fields
         IF ( MOD(NHMS,010000) .EQ. 0 ) THEN
            IF ( LSVMET     ) CALL SAVE_A3
            IF ( LSVMETNEST ) CALL SAVE_A3_NESTED
            CALL INIT_A3
         ENDIF

         ! Save A6 fields
         IF ( MOD(NHMS,010000) .EQ. 0 ) THEN
            IF ( LSVMET     ) CALL SAVE_A6
            IF ( LSVMETNEST ) CALL SAVE_A6_NESTED
            CALL INIT_A6
         ENDIF

      ENDIF

      ! Update daily mean temperature archive for MEGAN biogenics
      IF ( MOD(NHMS,010000) .EQ. 0 ) THEN
         IF ( LMEGAN ) CALL UPDATE_T_DAY
      ENDIF

      ! Cloud top level
      CALL MAKE_CLDTOPS

      ! Since CLDTOPS is an A-6 field, update the
      ! lightning NOx emissions [molec/box/6h]
      IF ( LLIGHTNOX ) CALL LIGHTNING

      !==============================================================
      ! ***** M O N T H L Y   O R   S E A S O N A L   D A T A *****
      !==============================================================

      ! UV albedoes
      IF ( LCHEM .and. ITS_A_NEW_MONTH() ) THEN
         CALL READ_UVALBEDO( MONTH )
      ENDIF

      ! Fossil fuel emissions (SMVGEAR)
      ! NOTE: Only call ANTHROEMS if it's time to do emissions;
      ! this will function properly for the central chemistry time
      ! step
      ! (bmy, 10/14/11)
      IF ( ITS_A_FULLCHEM_SIM() ) THEN
         IF ( ITS_TIME_FOR_EMIS() ) THEN
            IF ( LEMIS .and. ITS_A_NEW_SEASON() ) THEN
               CALL ANTHROEMS( SEASON )
            ENDIF
         ENDIF
      ENDIF

      !==============================================================
      !              ***** D A I L Y   D A T A *****
      !==============================================================
      IF ( ITS_A_NEW_DAY() ) THEN

         ! Read leaf-area index (needed for drydep)
         ! Now include the year (mpb,11/19/09)
         CALL RDLAI( DAY_OF_YEAR, MONTH, YEAR )

         ! For MEGAN biogenics ...
         IF ( LMEGAN ) THEN

            ! Read daily leaf-area-index
            ! Now include the year (mpb,11/19/09)
            CALL RDISOLAI( DAY_OF_YEAR, MONTH, YEAR )

            ! Compute 15-day average temperature for MEGAN
            IF ( NHMS == 000000 ) CALL UPDATE_T_15_AVG
         ENDIF

         ! Also read soil-type info for fullchem simulation
         IF ( ITS_A_FULLCHEM_SIM() ) THEN
            CALL RDSOIL
         ENDIF

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG ( '### MAIN: a DAILY DATA' )
      ENDIF

      ! Tropopause [hPa]
      !CALL TPP_REICHLER( T, PSC2, TROPP )

      ! Case of variable tropopause:
      ! Check LLTROP and set LMIN, LMAX, and LPAUSE
      ! since this is not done with READ_TROPOPAUSE anymore.
      ! (Need to double-check that LMIN, Lmax are not used before-phs) 
      !IF ( LVARTROP ) CALL CHECK_VAR_TROP

      ! If we are not doing transport, then make sure that
      ! the floating pressure is set to PSC2 (bdf, bmy, 8/22/02)
      IF ( .not. LTRAN ) CALL SET_FLOATING_PRESSURE( PSC2 )

      ! Compute airmass quantities at each grid box 
      IF ( .not. LTRAN ) CALL AIRQNT

      ! Compute LMIN, LMAX
      IF ( .not. LTRAN ) CALL CHECK_VAR_TROP

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
      !   ***** U N I T   C O N V E R S I O N  ( kg -> v/v ) *****
      !==============================================================
      IF ( ITS_TIME_FOR_UNIT() ) THEN
         CALL CONVERT_UNITS( 1,  N_TRACERS, TCVV, AD, STT )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:1' )
      ENDIF

#ifndef CHEM_NEWTP
      !==============================================================
      !     ***** S T R A T O S P H E R I C   F L U X E S *****
      !==============================================================
      IF ( LUPBD ) CALL DO_UPBDFLX
#endif

      !==============================================================
      !              ***** T R A N S P O R T *****
      !==============================================================
      IF ( ITS_TIME_FOR_DYN() ) THEN

#ifndef CHEM_NEWTP
         ! Call the appropritate version of TPCORE
         ierr = tprof_s('    transport', 2)
         IF ( LTRAN ) CALL DO_TRANSPORT
         ierr = tprof_e('    transport', 2)
#else
         ! We do transport in dyns_slq_driver
         STT = STT_TP
         CALL SET_FLOATING_PRESSURE( PSC2 )
#endif

         ! Reset air mass quantities
         CALL AIRQNT

         ! Compute LMIN, LMAX
         CALL CHECK_VAR_TROP

         ! Repartition [NOy] species after transport
         IF ( LUPBD .and. ITS_A_FULLCHEM_SIM() ) THEN
            CALL UPBDFLX_NOY( 2 )
         ENDIF

         ! Get relative humidity (after recomputing pressures)
         ! NOTE: for GEOS-5 we'll read this from disk instead
         CALL MAKE_RH

         ! Initialize wet scavenging and wetdep fields after
         ! the airmass quantities are reset after transport
         IF ( LCONV .or. LWETD ) CALL INIT_WETSCAV
      ENDIF


      !-------------------------------
      ! Test for emission timestep
      !-------------------------------
      IF ( LNLPBL .AND. ITS_TIME_FOR_EMIS() ) THEN

         IF ( ITS_TIME_FOR_UNIT() )
     &      CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT ) ! v/v -> kg

         ! Increment emission counter
         CALL SET_CT_EMIS( INCREMENT=.TRUE. )

         !========================================================
         !         ***** D R Y   D E P O S I T I O N *****
         !========================================================
         ierr = tprof_s('    drydep', 2)
         IF ( LDRYD ) CALL DO_DRYDEP
         ierr = tprof_e('    drydep', 2)

         !========================================================
         !             ***** E M I S S I O N S *****
         !========================================================
         ierr = tprof_s('    emission', 2)
         IF ( LEMIS ) CALL DO_EMISSIONS
         ierr = tprof_e('    emission', 2)

         IF ( ITS_TIME_FOR_UNIT() )
     &      CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT ) ! kg -> v/v

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
         ! Add option for non-local PBL. (Lin, 03/31/09)
         ierr = tprof_s('    pblmix', 2)
         IF ( .NOT. LNLPBL ) THEN
            CALL DO_PBL_MIX( LTURB )
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:1' )
         ELSE
            CALL DO_PBL_MIX_2( LTURB )
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a NLPBL 2' )
         ENDIF

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:2' )
         ierr = tprof_e('    pblmix', 2)

         !===========================================================
         !        ***** C L O U D   C O N V E C T I O N *****
         !===========================================================
         ierr = tprof_s('    convection', 2)
         IF ( LCONV ) THEN

            CALL DO_CONVECTION

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVECTION' )
         ENDIF
         ierr = tprof_e('    convection', 2)
      ENDIF


      !==============================================================
      !    ***** U N I T   C O N V E R S I O N  ( v/v -> kg ) *****
      !==============================================================
      IF ( ITS_TIME_FOR_UNIT() ) THEN
         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:2' )
      ENDIF

      !-------------------------------
      ! Test for emission timestep
      !-------------------------------
      IF ( ( .NOT. LNLPBL ) .AND. ITS_TIME_FOR_EMIS() ) THEN
            
         ! Increment emission counter
         CALL SET_CT_EMIS( INCREMENT=.TRUE. )
            
         !========================================================
         !         ***** D R Y   D E P O S I T I O N *****
         !========================================================
         ierr = tprof_s('    drydep', 2)
         IF ( LDRYD ) CALL DO_DRYDEP
         ierr = tprof_e('    drydep', 2)
            
         !========================================================
         !             ***** E M I S S I O N S *****
         !========================================================
         ierr = tprof_s('    emission', 2)
         IF ( LEMIS ) CALL DO_EMISSIONS
         ierr = tprof_e('    emission', 2)
      ENDIF

         
      !===========================================================
      !               ***** C H E M I S T R Y *****
      !===========================================================    
         
      ! Every chemistry timestep...
      IF ( ITS_TIME_FOR_CHEM() ) THEN
            
         ! Increment chemistry timestep counter
         CALL SET_CT_CHEM( INCREMENT=.TRUE. )
            
         ! Call the appropriate chemistry routine
         ierr = tprof_s('    chemistry', 2)
         CALL DO_CHEMISTRY
         ierr = tprof_e('    chemistry', 2)
         
      ENDIF
         
      !==============================================================
      ! ***** W E T   D E P O S I T I O N  (rainout + washout) *****
      !==============================================================
      IF ( LWETD .and. ITS_TIME_FOR_DYN() ) THEN
            
         ierr = tprof_s('    wetdep', 2)
         CALL DO_WETDEP
         ierr = tprof_e('    wetdep', 2)
            
      ENDIF

      !==============================================================
      !      ***** U P D A T E  O P T I C A L  D E P T H *****          
      !==============================================================
      ! Recalculate the optical depth at the wavelength specified by
      ! jv_spec_aod.dat This must be done before the call to any
      ! diagnostic and only on a chemistry timestep.
      ! (skim, 02/05/11)
      IF ( ITS_TIME_FOR_CHEM() ) THEN
         CALL RECOMPUTE_OD
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

         ! Accumulate GRIMs meteorology for offline GRIMs/GEOS-Chem
         IF ( LSVMET .or. LSVMETNEST ) CALL DO_METEO

         ! Increment diagnostic timestep counter. (ccc, 5/13/09)
         CALL SET_CT_DIAG( INCREMENT=.TRUE. )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a DIAGNOSTICS' )
      ENDIF

      ! Save I6 fields
      IF ( MOD(GET_NHMS(),010000) .EQ. 0 ) THEN
         IF ( LSVMET     ) CALL SAVE_I6( FIRST )
         IF ( LSVMETNEST ) CALL SAVE_I6_NESTED( FIRST )
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
         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:3' )
      ENDIF

      !==============================================================
      !  ***** E N D   O F   D Y N A M I C   T I M E S T E P *****
      !==============================================================
         
      ! Check for NaN, Negatives, Infinities in STT each time diag are
      ! saved. (ccc, 5/13/09)
      IF ( ITS_TIME_FOR_DIAG() ) THEN
         CALL CHECK_STT( 'End of Dynamic Loop' )
      ENDIF

      END SUBROUTINE INTEGRATE_CHEM
