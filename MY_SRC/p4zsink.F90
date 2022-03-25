MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  SMELT  vertical fluxes due to gravitational sinking
   !!======================================================================
   !! History :   2014 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp
   USE p4zprod, ONLY: zz_diat_Nlimit
   USE p4zriv, ONLY: wsink

#if defined key_trdtrc
   USE trdtrc
   USE trd_oce
#endif
#if defined key_oxy
   USE p4zoxy
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_sink_alloc

   REAL(wp), PUBLIC  ::   zz_w_sink_Pdiat_min  ! m/d  diatom minimum sinking rate # Alain
   REAL(wp), PUBLIC  ::   zz_w_sink_Pdiat_max  ! m/d diatom maximum sinking rate # Alain
   REAL(wp), PUBLIC  ::   zz_w_sink_D_PON       ! m/s PON detritus sinking rate # Jeffery quoting Dune and Bacon
   REAL(wp), PUBLIC  ::   zz_w_sink_D_bSi       !  m/s  biogenic silicon detritus sinking rate # match NO3 particles
   REAL(wp), PUBLIC  ::   zz_alpha_b_Si         !  fraction of sinking flux reflected back at bottom - Si
   REAL(wp), PUBLIC  ::   zz_alpha_b_D          !  fraction of sinking flux reflected back at bottom - diatoms
   REAL(wp), PUBLIC  ::   zz_alpha_b_N          !  fraction of sinking flux reflected back at bottom - N
   REAL(wp), PUBLIC  ::   zz_alpha_b_T          !  fraction of sinking flux reflected back at bottom - turbidity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  zz_biosink

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsink.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute diatom sinking and execute bottom flux
      !!
      !! ** Method  : - actual sinking now carried out within tvd advection
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk, ikt, jn

      REAL(wp), POINTER, DIMENSION(:,:,:) :: zz_diat_sink
#if defined key_trdtrc
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zz_bflux
#endif
      !!---------------------------------------------------------------------
      !

      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')

      CALL wrk_alloc(jpi, jpj, jpk, zz_diat_sink) 
#if defined key_trdtrc
      CALL wrk_alloc( jpi, jpj, jptra, zz_bflux )
      zz_bflux(:,:,:)=0.0_wp
#endif

      ! sinking defined by nutrient status
      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
              ! Calculate the sinking term for the quantities that sink
              ! to  calculate sinking at the interfaces used the values above
              IF (jk .LE. mbkt(ji,jj)) THEN
                 zz_diat_sink(ji,jj,jk) = -1._wp * (zz_w_sink_Pdiat_min * zz_diat_Nlimit(ji,jj,jk) &
                  &   + zz_w_sink_Pdiat_max * (1._wp - zz_diat_Nlimit(ji,jj,jk)))
              ENDIF 
            END DO
         END DO
      END DO
      zz_biosink(:,:,:,jpdia) = zz_diat_sink(:,:,:)

      ! bottom flux
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt = mbkt(ji,jj)
            tra(ji,jj,ikt,jpdia) = tra(ji,jj,ikt,jpdia) + zz_biosink(ji,jj,ikt,jpdia) * &
                                &        trb(ji,jj,ikt,jpdia) * (1_wp-zz_alpha_b_D) / fse3t(ji,jj,ikt)
            tra(ji,jj,ikt,jppon) = tra(ji,jj,ikt,jppon) + zz_biosink(ji,jj,ikt,jppon) * &
                                &        trb(ji,jj,ikt,jppon) * (1_wp-zz_alpha_b_N) / fse3t(ji,jj,ikt)
            tra(ji,jj,ikt,jpdsi) = tra(ji,jj,ikt,jpdsi) + zz_biosink(ji,jj,ikt,jpdsi) * &
                                &        trb(ji,jj,ikt,jpdsi) * (1_wp-zz_alpha_b_Si) / fse3t(ji,jj,ikt)
            tra(ji,jj,ikt,jpriv) = tra(ji,jj,ikt,jpriv) + zz_biosink(ji,jj,ikt,jpriv) * &
                                &        trb(ji,jj,ikt,jpriv) * (1_wp-zz_alpha_b_T) / fse3t(ji,jj,ikt)
#if defined key_oxy
            ! instantaneous full remineralization to no3 of a fraction of organic matter lost to sediments
            ! note: currently, this is returning oxygen to water column but not nitrate
            ! remember: zz_biosink has negative sign
            tra(ji,jj,ikt,jpo2) = tra(ji,jj,ikt,jpo2) + zz_alpha_SOD * (zz_o2_nitr+zz_o2_proreg) * &
                                &    ( zz_biosink(ji,jj,ikt,jpdia) * trb(ji,jj,ikt,jpdia) * (1_wp-zz_alpha_b_D) + &
                                &      zz_biosink(ji,jj,ikt,jppon) * trb(ji,jj,ikt,jppon) * (1_wp-zz_alpha_b_N) ) &
                                &         / fse3t(ji,jj,ikt)
#endif
#if defined key_trdtrc
            zz_bflux(ji,jj,jpdia) = zz_biosink(ji,jj,ikt,jpdia) * &
                                &        trb(ji,jj,ikt,jpdia) * (1_wp-zz_alpha_b_D) 
            zz_bflux(ji,jj,jppon) = zz_biosink(ji,jj,ikt,jppon) * &
                                &        trb(ji,jj,ikt,jppon) * (1_wp-zz_alpha_b_N) 
            zz_bflux(ji,jj,jpdsi) = zz_biosink(ji,jj,ikt,jpdsi) * &
                                &        trb(ji,jj,ikt,jpdsi) * (1_wp-zz_alpha_b_Si) 
            zz_bflux(ji,jj,jpriv) = zz_biosink(ji,jj,ikt,jpriv) * &
                                &        trb(ji,jj,ikt,jpriv) * (1_wp-zz_alpha_b_T) 
#endif
         END DO
      END DO

#if defined key_trdtrc 
      DO jn = 1, jptra
         CALL trd_flx_2d(zz_bflux(:,:,jn), jn, "BFX_","W")
      END DO
#endif

      CALL wrk_dealloc(jpi, jpj, jpk, zz_diat_sink)
#if defined key_trdtrc
      CALL wrk_dealloc( jpi, jpj, jptra, zz_bflux )
#endif
 
     IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink
   
     subroutine zz_sinking_advection(zz_w_sink, zz_flux, jptracer, knt, rfact) 
    ! Calculate the sinking term of the semi-implicit PDEs for the biology
    ! quantities that sink (some classes of plankton, and some of detritus).
    ! use upwind advection

    ! Arguments:
    real(kind=wp), dimension(jpi,jpj,jpk), intent(in) :: &
         zz_w_sink  ! Sinking velocity [m/s] on interfaces
    real(kind=wp), dimension(jpi,jpj,jpk), intent(out) :: &
         zz_flux  ! sinking flux
    integer, intent (in):: jptracer   ! which tracer
    integer, intent (in):: knt    ! not used
    real(kind=wp), intent (in):: rfact ! not used

    ! local variables
    INTEGER  ::   ji, jj, jk


!    zz_flux(:,:,0) = 0._wp ! no flux in or out from above
    
      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               if (zz_w_sink(ji,jj,jk).lt.0) then
                  write (*,*) "You have got biology sinking upward!"
                  call exit(1)
               else
                  zz_flux(ji,jj,jk) = - zz_w_sink(ji,jj,jk) * trn(ji,jj,jk,jptracer)
               endif
            END DO
         END DO
      END DO

      zz_flux(:,:,jpk)=0._wp
      jk = 1
      DO jj = 1, jpj
      	 DO ji = 1, jpi
                tra(ji,jj,jk,jptracer) = tra(ji,jj,jk,jptracer) - (- zz_flux(ji,jj,jk))/fse3t(ji,jj,jk)
      	END DO
      END DO
      DO jk = 2,jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
                tra(ji,jj,jk,jptracer) = tra(ji,jj,jk,jptracer) - (zz_flux(ji,jj,jk-1) - zz_flux(ji,jj,jk))/fse3t(ji,jj,jk)
            END DO
         END DO
      END DO
  end subroutine zz_sinking_advection


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      INTEGER  ::   ji, jj, jk
      !
      NAMELIST/nampissink/ zz_w_sink_Pdiat_min, zz_w_sink_Pdiat_max, zz_w_sink_D_PON, zz_w_sink_D_bSi, &
                           & zz_alpha_b_Si, zz_alpha_b_D, zz_alpha_b_N, zz_alpha_b_T


      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink_init')
      !

      REWIND( numnatp_ref )              ! Namelist nampissink in reference namelist : SOG sinking
      READ  ( numnatp_ref, nampissink, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissink in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampissink in configuration namelist : SOG sinking
      READ  ( numnatp_cfg, nampissink, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissink in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampissink )

      zz_w_sink_Pdiat_min=zz_w_sink_Pdiat_min/rday
      zz_w_sink_Pdiat_max=zz_w_sink_Pdiat_max/rday

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for sinking of biological material, nampissink'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '  m/d  diatom minimum sinking rate zz_w_sink_Pdiat_min*86400 =', zz_w_sink_Pdiat_min*rday
         WRITE(numout,*) '  m/d diatom maximum sinking rate  zz_w_sink_Pdiat_max*86400 =', zz_w_sink_Pdiat_max*rday
         WRITE(numout,*) '  m/s PON detritus sinking rate        zz_w_sink_D_PON =', zz_w_sink_D_PON
         WRITE(numout,*) '  m/s  bio si detritus sinking rate    zz_w_sink_D_bSi =', zz_w_sink_D_bSi
         WRITE(numout,*) '  BB reflection parameter Si           zz_alpha_b_Si =', zz_alpha_b_Si
         WRITE(numout,*) '  BB reflection parameter diatoms      zz_alpha_b_D =', zz_alpha_b_D
         WRITE(numout,*) '  BB reflection parameter N            zz_alpha_b_N =', zz_alpha_b_N
         WRITE(numout,*) '  BB reflection parameter turbidity    zz_alpha_b_T =', zz_alpha_b_T
      ENDIF

      zz_biosink(:,:,:,:)=0.0_wp

      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                IF (jk .LE. mbkt(ji,jj)) THEN
                    zz_biosink(ji,jj,jk,jppon)=-1._wp*zz_w_sink_D_PON
                    zz_biosink(ji,jj,jk,jpdsi)=-1._wp*zz_w_sink_D_bSi
                    zz_biosink(ji,jj,jk,jpriv)=-1._wp*wsink
                ENDIF
            END DO
         END DO
      END DO
 
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink_init')
      !
  END SUBROUTINE p4z_sink_init

   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( zz_biosink(jpi,jpj,jpk,jptra), STAT = p4z_sink_alloc )
      !
      IF( p4z_sink_alloc /= 0 ) CALL ctl_warn('p4z_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_sink_alloc   

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sink                    ! Empty routine
   END SUBROUTINE p4z_sink
#endif 
   !!======================================================================
END MODULE p4zsink
