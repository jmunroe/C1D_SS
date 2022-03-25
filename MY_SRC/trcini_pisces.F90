MODULE trcini_pisces
   !!======================================================================
   !!                         ***  MODULE trcini_pisces  ***
   !! TOP :   initialisation of the SMELT biochemical model
   !!======================================================================
   !! History :    -   !  1988-07  (E. Maier-Reiner) Original code
   !!              -   !  1999-10  (O. Aumont, C. Le Quere)
   !!              -   !  2002     (O. Aumont)  PISCES
   !!             1.0  !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.pisces.h90
   !!             3.5  !  2012-05  (C. Ethe) Merge PISCES-LOBSTER
   !!                     2015 (E. Olson) adapt for use with SMELT config
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_ini_pisces   : biochemical model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_pisces   ! called by trcini.F90 module


#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_pisces.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_pisces
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_pisces ***
      !!
      !! ** Purpose :   Initialisation of the biochemical model
      !!----------------------------------------------------------------------

      IF( lk_p4z ) THEN  ;   CALL p4z_ini   !  SMELT
      ELSE ;   CALL ctl_stop( 'STOP', 'trc_ini_pisces: SMELT lk_p4z must be true' )
      ENDIF

   END SUBROUTINE trc_ini_pisces

   SUBROUTINE p4z_ini
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_ini ***
      !!
      !! ** Purpose :   Initialisation of the SMELT biochemical model
      !!----------------------------------------------------------------------
      !
      USE p4zsms          ! Main P4Z routine
      USE p4zsink         !  vertical flux of particulate matter due to sinking
      USE p4zopt          !  optical model
      USE p4zsbc          !  Boundary conditions
      USE p4zrem          !  Remineralisation of organic matter
      USE p4zflx          !  Gas exchange
      USE p4zprod         !  Growth rate of the 2 phyto groups
      USE p4zmicro        !  Sources and sinks of microzooplankton
      USE p4zmeso         !  Sources and sinks of mesozooplankton
      USE p4zmesozoo         !  Sources and sinks of mesozooplankton
      USE p4zmort         !  Mortality terms for phytoplankton
      USE p4zriv       
#if defined key_skog
      USE p4zcar          ! carbon model 	
  !    USE p4zche          ! carbon constants
      USE mocsywrapper    ! allocation
#endif
#if defined key_oxy
      USE p4zoxy          ! oxygen component
#endif
      !
      REAL(wp), SAVE :: sco2   =  2.312e3_wp
      REAL(wp), SAVE :: alka0  =  2.426e3_wp
      REAL(wp), SAVE :: oxyg0  =  160.0_wp 
      REAL(wp), SAVE :: po4    =  2.165_wp 
      REAL(wp), SAVE :: bioma0 =  1.000e-2_wp  
      REAL(wp), SAVE :: silic1 =  91.51_wp  
      REAL(wp), SAVE :: no3    =  30.9_wp 
      !
      INTEGER  ::  ji, jj, jk, ierr
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' p4z_ini :   SMELT biochemical model initialisation'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

                                                 ! Allocate arrays
      ierr =         sms_pisces_alloc()          
      ierr = ierr +  p4z_opt_alloc()
      ierr = ierr +  p4z_prod_alloc()
      ierr = ierr +  p4z_sink_alloc()
#if defined key_skog
      ierr = ierr +  mocsy_alloc()
#endif
!
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'pisces_alloc: unable to allocate PISCES arrays' )
      !
      ryyss    = nyear_len(1) * rday    ! number of seconds per year
      r1_ryyss = 1. / ryyss
      !

      CALL p4z_sms_init       !  Maint routine
      !                                            ! Time-step

      !! Set biological ratios (PISCES version)
      !! ---------------------
      !rno3    =  1._wp/7.6_wp !16._wp / 122._wp redfield ratio
      !po4r    =   1._wp / 122._wp
      !o2nit   =  32._wp / 122._wp
      !rdenit  = 105._wp /  16._wp
      !rdenita =   3._wp /  5._wp
      !o2ut    = 133._wp / 122._wp

      ! Initialization of tracer concentration in case of  no restart 
      !--------------------------------------------------------------
      IF( .NOT. ln_rsttr ) THEN  
         trn(:,:,:,jpno3) = no3 * tmask(:,:,:)
         trn(:,:,:,jpnh4) = bioma0 * tmask(:,:,:)
         trn(:,:,:,jpsil) = silic1 * tmask(:,:,:)
         trn(:,:,:,jpdia) = bioma0 * tmask(:,:,:)
         trn(:,:,:,jpphy) = bioma0 * tmask(:,:,:)
         trn(:,:,:,jpmyr) = bioma0 * tmask(:,:,:)
         trn(:,:,:,jpzoo) = bioma0 * tmask(:,:,:)
         trn(:,:,:,jpdon) = bioma0 * tmask(:,:,:)
         trn(:,:,:,jppon) = bioma0 * tmask(:,:,:)
         trn(:,:,:,jpdsi) = bioma0 * 0.15 * tmask(:,:,:)
         trn(:,:,:,jpriv) = 0.2_wp * tmask(:,:,:)
#if defined key_skog
	 trn(:,:,:,jpdic) = no3 * tmask(:,:,:)
         trn(:,:,:,jpta)  = no3 * tmask(:,:,:)
#endif
#if defined key_oxy
         trn(:,:,:,jpo2)  = oxyg0 * tmask(:,:,:)
#endif
       END IF


      CALL p4z_riv_init       !  river tracer MUST BE BEFORE SINK INIT
      CALL p4z_sink_init      !  vertical flux of particulate organic matter
      CALL p4z_opt_init       !  Optic: PAR in the water column
      CALL p4z_prod_init      !  phytoplankton growth rate over the global ocean.
      CALL p4z_sbc_init       !  boundary conditions
      CALL FLUSH(numout)
      
      CALL p4z_rem_init       !  remineralisation
      CALL p4z_mort_init      !  phytoplankton mortality 
      CALL p4z_micro_init     !  M. rubra
      CALL p4z_meso_init      !  mesozooplankton
      CALL p4z_mesozoo_init   !  mesozooplankton
#if defined key_skog
      CALL p4z_car_init       !  dic
      CALL p4z_flx_init       !  air-sea exchange
#endif
#if defined key_oxy
      CALL p4z_oxy_init       ! read oxygen namelist
#endif
      ndayflxtr = 0

      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'Initialization of SMELT tracers done'
      IF(lwp) WRITE(numout,*) 
      !
   END SUBROUTINE p4z_ini

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_pisces             ! Empty routine
   END SUBROUTINE trc_ini_pisces
#endif

   !!======================================================================
END MODULE trcini_pisces
