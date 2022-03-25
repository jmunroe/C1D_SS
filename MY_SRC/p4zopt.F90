MODULE p4zopt
   !!======================================================================
   !!                         ***  MODULE p4zopt  ***
   !! TOP - SMELT : Compute the light availability in the water column
   !!======================================================================
   !! History :  2014 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined  key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                              SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_opt       : light availability in the water column
   !!----------------------------------------------------------------------
   USE trc            ! tracer variables
   USE oce_trc        ! tracer-ocean share variables
   USE sms_pisces     ! Source Minus Sink
   USE iom            ! I/O manager
   USE fldread         !  time interpolation
   USE prtctl_trc      !  print control for debugging


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_opt        ! called in p4zbio.F90 module
   PUBLIC   p4z_opt_init   ! called in trcsms_pisces.F90 module
   PUBLIC   p4z_opt_alloc

   !! * Shared module variables

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: etot, epar   !: total and  PAR 

   INTEGER  ::   nksrp   ! levels below which the light cannot penetrate ( depth larger than 391 m)
   !REAL(wp) ::   parlux = 0.43_wp / 3._wp

   ! qsr is shortware radiation read in through sbcblk_core.F90 in W/m^2; hourly, interpolated
   REAL(wp) ::   zzialpha
   REAL(wp) ::   zzibeta
   REAL(wp) ::   zzigamma
   REAL(wp) ::   zzisigma
   REAL(wp) ::   zzitheta
   REAL(wp) ::   zzidl
   REAL(wp) ::   zzN2chl


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zopt.F90 3160 2011-11-20 14:27:18Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_opt( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt  ***
      !!
      !! ** Purpose :   Compute the light availability in the water column
      !!              depending on the depth and the chlorophyll concentration
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      real(wp) :: &     ! Values for the Kpar fit
         zzKK, & ! for light atten. calc.
         zzQinter
      real(wp) :: zzIred, zzIblue ! light in red and blue parts of spectrum
      real(wp), POINTER, DIMENSION(:) :: zzabovefactor, zzbelowfactor ! used in interpolation (1:jpkm1)
      real(wp), POINTER, DIMENSION(:,:,:) :: zzI_par_i, zzI_total !(1:jpk)


      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_opt')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpkm1,         zzabovefactor, zzbelowfactor  )
      CALL wrk_alloc( jpi, jpj, jpk, zzI_par_i, zzI_total )


    zzQinter = 0.
    zzI_total(:,:,1) =  qsr(:,:) ! light at w grid points (interfaces of t cells)
    ! 44% of total light is PAR at surface (Jerlov)
    zzI_par_i(:,:,1) = qsr(:,:) * 0.44
    DO jj = 1, jpj
        DO ji = 1, jpi
            zzIblue = zzI_total(ji,jj,1) * 0.7
            zzIred = zzI_total(ji,jj,1) * 0.3
            do jk = 2, jpk
                zzKK = zzialpha +zzibeta * (zzN2chl &
                &    * (trn(ji,jj,jk-1,jpdia) + trn(ji,jj,jk-1,jpmyr) + trn(ji,jj,jk-1,jpphy))) ** 0.665 &
                    + (zzigamma * zzQinter ** zzisigma + zzitheta) * exp(-fsdept(ji,jj,jk-1) / zzidl)

                zzKK = min(2.5d0, zzKK)
                zzI_par_i(ji,jj,jk) = zzI_par_i(ji,jj,jk-1) * exp(-fse3w(ji,jj,jk-1) * zzKK)
                ! Total light for heat budget
                zzIblue = zzIblue * exp(-fse3w(ji,jj,jk-1) * (2.497d0 * zzKK + 1.355d0))
                zzIred = zzIred * exp(-fse3w(ji,jj,jk-1) * (0.8864d0 * zzKK))
                zzI_total(ji,jj,jk) = zzIblue + zzIred
            end do
            ! slight difference from SOG here: SOG has I_par(0) = I_par_i(0)
            zzabovefactor(1:jpkm1) = (fsdept(ji,jj,1:jpkm1)-fsdepw(ji,jj,1:jpkm1))/fse3t(ji,jj,1:jpkm1)
            zzbelowfactor(1:jpkm1) =(fsdepw(ji,jj,2:jpk)-fsdept(ji,jj,1:jpkm1))/fse3t(ji,jj,1:jpkm1)
            epar(ji,jj,1:jpkm1)=zzI_par_i(ji,jj,1:jpkm1)*zzabovefactor(1:jpkm1)+zzI_par_i(ji,jj,2:jpk)*zzbelowfactor(1:jpkm1)
            etot(ji,jj,1:jpkm1)=zzI_total(ji,jj,1:jpkm1)*zzabovefactor(1:jpkm1)+zzI_total(ji,jj,2:jpk)*zzbelowfactor(1:jpkm1)
         END DO
      END DO

      IF( ln_qsr_bio ) etot3 = zzI_total      !* heat flux accros w-level (used in the dynamics)
                                              !  ------------------------


      !                                        !* Euphotic depth and level
      neln(:,:) = 1                            !  ------------------------
      heup(:,:) = 300.
      nksrp=min(jpk,35) ! minimum is for AGRIF reduced vertical domain in case <35 levels

      DO jk = 2, nksrp
         DO jj = 1, jpj
           DO ji = 1, jpi
              IF( etot(ji,jj,jk) * tmask(ji,jj,jk) >= 0.0043 * qsr(ji,jj) )  THEN
                 neln(ji,jj) = jk+1                    ! Euphotic level : 1rst T-level strictly below Euphotic layer
                 !                                     ! nb: ensure the compatibility with nmld_trc definition in trd_mld_trc_zint
                 heup(ji,jj) = fsdepw(ji,jj,jk+1)      ! Euphotic layer depth
              ENDIF
           END DO
        END DO
      END DO
 
      heup(:,:) = MIN( 300., heup(:,:) )
           !IF( knt == nrdttrc ) THEN
           !   CALL iom_put( "Heup", heup(:,:  ) * tmask(:,:,1) )  ! euphotic layer deptht
           !ENDIF
      !
      CALL wrk_dealloc( jpkm1,         zzabovefactor, zzbelowfactor  )
      CALL wrk_dealloc( jpi, jpj, jpk, zzI_par_i, zzI_total )
    
      IF( nn_timing == 1 )  CALL timing_stop('p4z_opt')
      !
   END SUBROUTINE p4z_opt

   SUBROUTINE p4z_opt_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of tabulated attenuation coef
      !!                and of the percentage of PAR in Shortwave
      !!
      !! ** Input   :   external ascii and netcdf files
      !!----------------------------------------------------------------------
      !
      INTEGER :: ios                 ! Local integer output status for namelist read
      NAMELIST/nampisopt/ zzialpha, zzibeta, zzigamma, zzisigma, &
         &   zzitheta, zzidl, zzN2chl
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('p4z_opt_init')

      REWIND( numnatp_ref )              ! Namelist nampisopt in reference namelist : attenuation coef. and PAR
      READ  ( numnatp_ref, nampisopt, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisopt in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisopt in configuration namelist : attenuation coef. and PAR
      READ  ( numnatp_cfg, nampisopt, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisopt in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisopt )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist : nampisopt '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '      zzialpha   = ', zzialpha
         WRITE(numout,*) '      zzibeta    = ', zzibeta
         WRITE(numout,*) '      zzigamma   = ', zzigamma
         WRITE(numout,*) '      zzisigma   = ', zzisigma
         WRITE(numout,*) '      zzitheta   = ', zzitheta
         WRITE(numout,*) '      zzidl      = ', zzidl
         WRITE(numout,*) '      zzN2chl    = ', zzN2chl
      ENDIF
    
                         etot (:,:,:) = 0._wp
                         epar(:,:,:) = 0._wp
      IF( ln_qsr_bio )   etot3(:,:,:) = 0._wp
      
      IF( nn_timing == 1 )  CALL timing_stop('p4z_opt_init')
      !
   END SUBROUTINE p4z_opt_init


   INTEGER FUNCTION p4z_opt_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( epar(jpi,jpj,jpk), &
        &       etot(jpi,jpj,jpk), STAT=p4z_opt_alloc ) 
         !
      IF( p4z_opt_alloc /= 0 ) CALL ctl_warn('p4z_opt_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_opt_alloc

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE p4z_opt                   ! Empty routine
   END SUBROUTINE p4z_opt
#endif 
   !!======================================================================
END MODULE p4zopt
