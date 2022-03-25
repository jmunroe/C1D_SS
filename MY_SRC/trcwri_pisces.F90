MODULE trcwri_pisces
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    SMELT  :   Output of SMELT tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && defined key_iomput && defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                  SMELT / PISCES model
   !!----------------------------------------------------------------------
   !! trc_wri_pisces   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE oce_trc     !  shared variables between ocean and passive tracers
   USE trc         ! passive tracers common variables 
   USE sms_pisces  ! PISCES variables
   USE iom         ! I/O manager
   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_pisces 

#  include "top_substitute.h90"
CONTAINS

   SUBROUTINE trc_wri_pisces
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)           :: cltra, strU, strV, strW
      CHARACTER (len=20)           :: cltraU, cltraV, cltraW
      REAL(wp)                     :: zfact
      INTEGER                      :: ji, jj, jk, jn
      REAL(wp), POINTER, DIMENSION(:,:,:) :: z3d      ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj) :: zdic, zo2min, zdepo2min

      CALL wrk_alloc( jpi , jpj, jpk , z3d )
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      strU="_UT"
      strV="_VT"
      strW="_WT"
      DO jn = jp_pcs0, jp_pcs1
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         cltraU=TRIM(cltra)//TRIM(strU)
         cltraV=TRIM(cltra)//TRIM(strV)
         cltraW=TRIM(cltra)//TRIM(strW)
         IF( iom_use( cltra ) )  CALL iom_put( cltra, trn(:,:,:,jn) )

         IF( iom_use( cltraU )) THEN ! calculate u-dir transport
            z3d(:,:,:) = 0.e0
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d(ji,jj,jk) = un(ji,jj,jk) * e2u(ji,jj) * fse3u_n(ji,jj,jk) * &
                        & umask(ji,jj,jk) * ( trn(ji,jj,jk,jn) + trn(ji+1,jj,jk,jn) )
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z3d, 'U', -1. )
            CALL iom_put( cltraU , 0.5 * z3d )  
         ENDIF
      
         IF( iom_use( cltraV )) THEN ! calculate v-dir transport
            z3d(:,:,:) = 0.e0
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d(ji,jj,jk) = vn(ji,jj,jk) * e1v(ji,jj) * fse3v_n(ji,jj,jk) * &
                        & vmask(ji,jj,jk) * ( trn(ji,jj,jk,jn) + trn(ji,jj+1,jk,jn) )
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z3d, 'V', -1. )
            CALL iom_put( cltraV, 0.5 * z3d )  
         ENDIF
      
         IF( iom_use( cltraW )) THEN ! calculate u-dir transport
            z3d(:,:,:) = 0.e0
            DO jk = 2, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d(ji,jj,jk) = wn(ji,jj,jk) * e1t(ji,jj) * e2t(ji,jj) * &
                        & wmask(ji,jj,jk) * ( trn(ji,jj,jk,jn) + trn(ji,jj,jk-1,jn) )
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z3d, 'W', -1. )
            CALL iom_put( cltraW , 0.5 * z3d )  
         ENDIF
      
      END DO
      
      !
      !IF( iom_use( "O2MIN" ) .OR. iom_use ( "ZO2MIN" ) ) THEN  ! Oxygen minimum concentration and depth 
      !   zo2min   (:,:) = trn(:,:,1,jpriv) * tmask(:,:,1)
      !   zdepo2min(:,:) = fsdepw(:,:,1)    * tmask(:,:,1)
      !   DO jk = 2, jpkm1
      !      DO jj = 1, jpj
      !         DO ji = 1, jpi
      !            IF( tmask(ji,jj,jk) == 1 ) then
      !               IF( trn(ji,jj,jk,jpriv) < zo2min(ji,jj) ) then
      !                  zo2min   (ji,jj) = trn(ji,jj,jk,jpriv)
      !                  zdepo2min(ji,jj) = fsdepw(ji,jj,jk)
      !               ENDIF
      !            ENDIF
      !         END DO
      !      END DO
      !   END DO
         !
      !   CALL iom_put('O2MIN' , zo2min     )                              ! oxygen minimum concentration
      !   CALL iom_put('ZO2MIN', zdepo2min  )                              ! depth of oxygen minimum concentration
          !
      !ENDIF
      !
      CALL wrk_dealloc( jpi , jpj, jpk , z3d )
   END SUBROUTINE trc_wri_pisces

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_pisces
CONTAINS
   SUBROUTINE trc_wri_pisces                     ! Empty routine  
   END SUBROUTINE trc_wri_pisces
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_pisces.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_pisces
