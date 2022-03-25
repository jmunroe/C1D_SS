MODULE trdtrc
   !!======================================================================
   !!                       ***  MODULE  trdtrc  ***
   !! Ocean diagnostics:  mixed layer passive tracer trends 
   !!======================================================================
   !! History :  3.0  !  2010-07  (C. Ethe)  Original code (from trdtrc.F90)
   !!----------------------------------------------------------------------
#if   defined key_top &&  defined key_trdtrc 
   !!----------------------------------------------------------------------
   !!   'key_trdtrc'                      3D trend diagnostics
   !!----------------------------------------------------------------------
   !!   trdtrc      : passive tracer trends 
   !!----------------------------------------------------------------------
   USE trc               ! tracer definitions (trn, trb, tra, etc.)
   USE trcnam_trp
   USE trd_oce
   USE iom               ! I/O library
   USE sms_pisces

   IMPLICIT NONE
   PRIVATE

   INTERFACE trd_trc
      MODULE PROCEDURE trd_trc_trp, trd_trc_bio
   END INTERFACE

   PUBLIC trd_trc
   PUBLIC trd_trc_novol
   PUBLIC trd_flx ! fluxes: make CHARACTER flag for U,V,W
   PUBLIC trd_flx_2d ! fluxes: make CHARACTER flag for U,V,W

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trdtrc.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trd_trc_trp( ptrtrd, kjn, ktrd, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_trc  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      INTEGER, INTENT( in )  ::   ktrd                                ! tracer trend index
      REAL(wp), DIMENSION(:,:,:), INTENT( in )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER (len=20) :: cltra
      !!----------------------------------------------------------------------

      IF( l_trdtrc .AND. ln_trdtrc( kjn ) ) THEN
         !
         SELECT CASE( ktrd )
         CASE( jptra_xad  )       ;    WRITE (cltra,'("XAD_",4a)')
         CASE( jptra_yad  )       ;    WRITE (cltra,'("YAD_",4a)')
         CASE( jptra_zad  )       ;    WRITE (cltra,'("ZAD_",4a)')
         CASE( jptra_ldf  )       ;    WRITE (cltra,'("LDF_",4a)')
         CASE( jptra_bbl  )       ;    WRITE (cltra,'("BBL_",4a)')
         CASE( jptra_nsr  )       ;    WRITE (cltra,'("FOR_",4a)')
         CASE( jptra_zdf  )       ;    WRITE (cltra,'("ZDF_",4a)')
         CASE( jptra_dmp  )       ;    WRITE (cltra,'("DMP_",4a)')
         CASE( jptra_sms  )       ;    WRITE (cltra,'("SMS_",4a)')
         CASE( jptra_atf  )       ;    WRITE (cltra,'("ATF_",4a)')
         CASE( jptra_radb )       ;    WRITE (cltra,'("RDB_",4a)')
         CASE( jptra_radn )       ;    WRITE (cltra,'("RDN_",4a)')
         END SELECT
         cltra = TRIM(cltra)//TRIM(ctrcnm(kjn))
         IF( iom_use(cltra))                 CALL iom_put( cltra,  ptrtrd(:,:,:)*cvol(:,:,:)*rn_ucf_trc )
         IF (lwp) WRITE(numout,*) "Called trd_trc_trp", cltra
         !
      END IF

   END SUBROUTINE trd_trc_trp

   SUBROUTINE trd_trc_novol(ptrtrd, kjn, cltra0, kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_bio  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      REAL(wp), DIMENSION(:,:,:), INTENT( in )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER (len=20), INTENT( in)  :: cltra0
      CHARACTER (len=20)  :: cltra
      !!----------------------------------------------------------------------
      !! for case where volume already included in diagnostic
         cltra = TRIM(cltra0)//TRIM(ctrcnm(kjn))
         IF( iom_use(cltra))                 CALL iom_put( cltra,  ptrtrd(:,:,:)*tmask(:,:,:)*rn_ucf_trc )
         IF (lwp) WRITE(numout,*) "Called trd_trc_novol", cltra
   END SUBROUTINE trd_trc_novol

   SUBROUTINE trd_trc_bio(ptrtrd, kjn, cltra0, kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_bio  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      REAL(wp), DIMENSION(:,:,:), INTENT( in )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER (len=20), INTENT( in)  :: cltra0
      CHARACTER (len=20)  :: cltra
      !!----------------------------------------------------------------------
         cltra = TRIM(cltra0)//TRIM(ctrcnm(kjn))
         IF( iom_use(cltra))                 CALL iom_put( cltra,  ptrtrd(:,:,:)*cvol(:,:,:)*tmask(:,:,:)*rn_ucf_trc )
         IF (lwp) WRITE(numout,*) "Called trd_trc_bio", cltra
   END SUBROUTINE trd_trc_bio

   SUBROUTINE trd_flx( ptrtrd, kjn, prfx, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_trc  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      REAL(wp), DIMENSION(:,:,:), INTENT( in )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER(len=1), INTENT(in   )           ::   cd_type   ! nature of ptrtrd grid-points ('U','V','W')
      CHARACTER(len=4), INTENT(in   )           ::   prfx   ! prefix for variable name (XXX_)
      CHARACTER (len=20) :: cltra
      !!----------------------------------------------------------------------
      ! trends ptrtrd should already be in units of mmol/s
      IF( l_trdtrc .AND. ln_trdtrc( kjn ) ) THEN
         !
         cltra = TRIM(prfx)//TRIM(ctrcnm(kjn))
         IF( iom_use(cltra)) THEN      
            SELECT CASE( cd_type )
            CASE( "U"  )       ;      CALL iom_put( cltra,  ptrtrd(:,:,:)*umask(:,:,:)*rn_ucf_trc)
            CASE( "V"  )       ;      CALL iom_put( cltra,  ptrtrd(:,:,:)*vmask(:,:,:)*rn_ucf_trc)
            CASE( "W"  )       ;      CALL iom_put( cltra,  ptrtrd(:,:,:)*wmask(:,:,:)*rn_ucf_trc)
            END SELECT
         ENDIF
         IF (lwp) WRITE(numout,*) "Called trd_flux", cltra
         !
      END IF
   END SUBROUTINE trd_flx

   SUBROUTINE trd_flx_2d( ptrtrd, kjn, prfx, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_trc  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      REAL(wp), DIMENSION(:,:), INTENT( in )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER(len=1), INTENT(in   )           ::   cd_type   ! nature of ptrtrd grid-points ('U','V','W')
      CHARACTER(len=4), INTENT(in   )           ::   prfx   ! prefix for variable name (XXX_)
      CHARACTER (len=20) :: cltra
      !!----------------------------------------------------------------------
      ! trends ptrtrd should already be in units of mmol/s
      IF( l_trdtrc .AND. ln_trdtrc( kjn ) ) THEN
         !
         cltra = TRIM(prfx)//TRIM(ctrcnm(kjn))
         IF( iom_use(cltra)) THEN      
            SELECT CASE( cd_type )
            CASE( "U"  )       ;      CALL iom_put( cltra,  ptrtrd(:,:)*umask(:,:,1)*rn_ucf_trc)
            CASE( "V"  )       ;      CALL iom_put( cltra,  ptrtrd(:,:)*vmask(:,:,1)*rn_ucf_trc)
            CASE( "W"  )       ;      CALL iom_put( cltra,  ptrtrd(:,:)*wmask(:,:,1)*rn_ucf_trc)
            END SELECT
         ENDIF
         IF (lwp) WRITE(numout,*) "Called trd_flx_2d", cltra
         !
      END IF
   END SUBROUTINE trd_flx_2d

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------

   INTERFACE trd_trc
      MODULE PROCEDURE trd_trc_trp, trd_trc_bio
   END INTERFACE

   PUBLIC trd_flx ! fluxes: make CHARACTER flag for U,V,W
   PUBLIC trd_trc_novol
   PUBLIC trd_flx_2d 
CONTAINS

   SUBROUTINE trd_trc_trp( ptrtrd, kjn, ktrd, kt )
      INTEGER               , INTENT( in )     ::   kt      ! time step
      INTEGER               , INTENT( in )     ::   kjn     ! tracer index
      INTEGER               , INTENT( in )     ::   ktrd    ! tracer trend index
      REAL, DIMENSION(:,:,:), INTENT( inout )  ::   ptrtrd  ! Temperature or U trend
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrtrd(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kt
   END SUBROUTINE trd_trc_trp

   SUBROUTINE trd_trc_bio(ptrtrd, kjn, cltra0, kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_bio  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      REAL, DIMENSION(:,:,:), INTENT( in )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER (len=20), INTENT( in)  :: cltra0
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrtrd(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', cltra0
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kt
   END SUBROUTINE trd_trc_bio

   SUBROUTINE trd_flx( ptrtrd, kjn, prfx, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_trc  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      REAL, DIMENSION(:,:,:), INTENT( in )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER(len=1), INTENT(in   )           ::   cd_type   ! nature of ptrtrd grid-points ('U','V','W')
      CHARACTER(len=4), INTENT(in   )           ::   prfx   ! prefix for variable name (XXX_)
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrtrd(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', prfx
   END SUBROUTINE trd_flx

   SUBROUTINE trd_trc_novol(ptrtrd, kjn, cltra0, kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_bio  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      REAL, DIMENSION(:,:,:), INTENT( in )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER (len=20), INTENT( in)  :: cltra0
      CHARACTER (len=20)  :: cltra
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrtrd(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', cltra0
   END SUBROUTINE trd_trc_novol

   SUBROUTINE trd_flx_2d( ptrtrd, kjn, prfx, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_trc  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      REAL, DIMENSION(:,:), INTENT( in )  ::   ptrtrd ! Temperature or U trend
      CHARACTER(len=1), INTENT(in   )           ::   cd_type   ! nature of ptrtrd grid-points ('U','V','W')
      CHARACTER(len=4), INTENT(in   )           ::   prfx   ! prefix for variable name (XXX_)
      CHARACTER (len=20) :: cltra
      !!----------------------------------------------------------------------
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrtrd(1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', prfx
   END SUBROUTINE trd_flx_2d

#endif
   !!======================================================================
END MODULE trdtrc
