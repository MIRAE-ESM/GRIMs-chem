!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdsoil
!
! !DESCRIPTION: Subroutine RDSOIL reads in soiltype data, fertilizer data, 
!  and monthly soil precipitation data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDSOIL
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,      ONLY : GET_MONTH
      
      USE CMN_SIZE_MOD   ! Size parameters
      USE COMMSOIL_MOD ! Soil variables

      IMPLICIT NONE
#     include "define.h"

!
! !REMARKS:
!  RDSOIL is one of the original GEOS-CHEM subroutines, and has its origins
!  from the GISS-II model that was used at Harvard in the early 90's.  This
!  was cleaned up and improved error checking was added. (bmy, 4/2/02)
!                                                                             .
!  Variables from "commsoil.h" header file:
!  ============================================================================
!  (1 ) NCONSOIL  (INTEGER) : Olson -> soil type mapping index                
!  (2 ) INDEXSOIL (INTEGER) : Array containing grid box indices (I,J)
!  (3 ) SOILFERT  (REAL*8 ) : Array containing fertilizer NOx [ng N/m2/s]  
!  (4 ) SOILPREP  (REAL*8 ) : Array containing 2 months of observed
!                              soil precipitation [mm/day]    
!                                                                             .
!  Files read in by "rdsoil.f":
!  ============================================================================
!  (1 ) DATA_DIR/soil_NOx_200203/soiltype.dat       : Olson and soil land types
!  (2 ) DATA_DIR/soil_NOx_200203/fert_scale.dat     : NOx from fertilizers
!  (3 ) DATA_DIR/soil_NOx_200203/climatprep4x5.dat  : 1x1   monthly soil precip
!                                climatprep2x25.dat : 2x2.5 monthly soil precip
!                                climatprep1x1.dat  : 4x5   monthly soil precip
! 
! !REVISION HISTORY: 
!  05 Jan 1994 - Y. H. Wang, G. M. Gardner, - Initial version
!  (1 ) Be sure to force double precision with the DBLE function and the "D" 
!        exponent, wherever necessary (bmy, 10/6/99)            *
!  (2 ) Now read soil data files directly from the from 
!        DATA_DIR/soil_NOx_200203/ subdirectory.  Now use IOERROR to trap
!        I/O errors across all platforms.  Added comment header.  Updated
!        comments, cosmetic changes. (bmy, 4/2/02)
!  (3 ) Removed obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE as the file unit number,
!        assign it to IUNIT. (bmy, 6/27/02)
!  (4 ) Now reference GEOS_CHEM_STOP from "error_mod.f".  Bug fix: remove 
!        duplicate declaration of IOS.  This causes compile errors for the 
!        ALPHA platform. (gcc, bmy, 11/6/02)
!  (5 ) Now use function GET_MONTH from "time_mod.f".  Now make MONTH a local
!        variable. (bmy, 2/11/03)
!  (6 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST   = .TRUE.
      INTEGER, SAVE      :: MONSAVE = 0 
      INTEGER            :: I, IUNIT, IOS, J, K, KK, M, M1, MONTH
      REAL*8             :: TMP(12)

      ! Name of file to read in
      CHARACTER(LEN=255) :: FILENAME

      ! Arrays for transform
      INTEGER, PARAMETER :: NLAND2 = 3920
      INTEGER            :: F_INDEX(I2x25,J2x25)
      INTEGER            :: P_INDEX(IIPAR,JJPAR)
      REAL*8             :: SOILPREP2(2,NLAND2)
      REAL*8             :: SOILFERT2(NLAND2)

      !=================================================================
      ! RDSOIL begins here!
      !=================================================================

      ! Define the file unit
      IUNIT = IU_FILE
      
      ! Get the current month
      MONTH = GET_MONTH()

      ! First-time only initialization
      IF ( FIRST ) THEN
         
         ! Reset First-time flag
         FIRST = .FALSE.

         !==============================================================
         ! Read in soil type data (first pass only)
         !==============================================================

         ! Define soiltype file name
         FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_200203/soiltype.dat'

         ! Echo filename
         IF (LMASTER) THEN
         WRITE( 6, 100 ) TRIM( FILENAME )
         ENDIF
 100     FORMAT( '     - RDSOIL: Reading ', a )

         IF (LMASTER) THEN
         ! Open file
         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:1' )

         ! Read header line
         READ( IUNIT, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:2' )

         ! Read data
         DO K = 1, NVEGTYPE
	    READ( IUNIT, *, IOSTAT=IOS ) KK, NCONSOIL(KK+1) 
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:3' )
         ENDDO

         ! Close file
         CLOSE( IUNIT )
         ENDIF

#ifdef MP
#ifndef RMP
         CALL mpbcasti(NCONSOIL,NVEGTYPE)
#else
         CALL rmpbcasti(NCONSOIL,NVEGTYPE)
#endif
#endif

         !==============================================================
         ! Read in fertilizer data (first pass only)
         ! Units are [ng N/m2/s] 
         !==============================================================

         ! Define fertilizer file name
         FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_200203/fert_scale.dat'

         ! Echo filename
         IF (LMASTER) THEN
         WRITE( 6, 100 ) TRIM( FILENAME )
         ENDIF

         IF (LMASTER) THEN
         ! Open file
         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:4' )

         ! Read data -- save (I,J) pairs into INDEXSOIL array
         F_INDEX = 0
         DO M = 1, NLAND2
	    READ( IUNIT,*, IOSTAT=IOS ) I, J, SOILFERT2(M)
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:5' )
            F_INDEX(I,J) = M
         ENDDO

         ! Close file
         CLOSE( IUNIT )
         ENDIF

         CALL TRANS_INDEX( F_INDEX, P_INDEX, I2x25, J2x25 )

#ifdef MP
#ifndef RMP
         CALL mpbcastr(SOILFERT2,NLAND2)
#else
         CALL rmpbcastr(SOILFERT2,NLAND2)
#endif
#endif

         NLAND = 0
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            M = P_INDEX(I,J)
            IF (M .GT. 0) THEN
               NLAND = NLAND + 1
            ENDIF
         ENDDO
         ENDDO

         CALL INIT_COMMSOIL

         K = 0
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            M = P_INDEX(I,J)
            IF (M .GT. 0) THEN
               K = K + 1
               INDEXSOIL(1,K) = I
               INDEXSOIL(2,K) = J
               SOILFERT(K) = SOILFERT2(M)
            ENDIF
         ENDDO
         ENDDO

      ENDIF
         
      !=================================================================
      ! Read in monthly soil precipitation data
      !=================================================================

      ! Only read data when we have entered a new month...
      IF ( MONSAVE /= MONTH ) THEN
      
         ! Save the current month 
         MONSAVE = MONTH
      
         ! M1 is the previous month
         IF ( MONTH == 1 ) THEN
            M1 = 12
         ELSE
            M1 = MONTH - 1
         END IF

         ! Define soil precip file name
         FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_200203/climatprep' //
     &              GET_RES_EXT()    // '.dat'

         ! Echo filename
         IF (LMASTER) THEN
         WRITE( 6, 100 ) TRIM( FILENAME )
         ENDIF

         IF (LMASTER) THEN
         ! Open soil precip file
         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:6' )

         F_INDEX = 0
         ! Loop over Olson land types
         DO M = 1, NLAND2

            ! Read monthly soil precip data for each (I,J) box
            READ( IUNIT, *, IOSTAT=IOS ) I, J, ( TMP(K), K=1,12 )
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:7' )
            F_INDEX(I,J) = M
            SOILPREP2(1,M) = TMP(M1)
            SOILPREP2(2,M) = TMP(MONTH)

         ENDDO

         ! Close file
         CLOSE( IUNIT )
         ENDIF

         CALL TRANS_INDEX( F_INDEX, P_INDEX, I2x25, J2x25 )

#ifdef MP
#ifndef RMP
         CALL mpbcastr(SOILPREP2,2*NLAND2)
#else
         CALL rmpbcastr(SOILPREP2,2*NLAND2)
#endif
#endif

         K = 0
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            M = P_INDEX(I,J)
            IF (M .GT. 0) THEN
               K = K + 1
               ! Error check -- make sure that each (I,J) pair has both
               ! soil precip data and fertilizer data defined
               IF ( INDEXSOIL(1,K) /= I .OR. INDEXSOIL(2,K) /= J ) THEN
                  WRITE(6,*) 'CORRUPTED TEMPCORR OR CLMATPRECIP DATA'
                  WRITE(6,*) 'CHECK (I,J)',I,J
                  CALL GEOS_CHEM_STOP
               ELSE
                  SOILPREP(1,K) = SOILPREP2(1,M)
                  SOILPREP(2,K) = SOILPREP2(2,M)
               ENDIF
            ENDIF
         ENDDO
         ENDDO

      ENDIF

      ! Return to calling program
      END SUBROUTINE RDSOIL






