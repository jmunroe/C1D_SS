MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   SMELT Compute the sources/sinks for M. rubrum due to 
   !!          heterotrophy
   !!======================================================================
   !! History :   2015 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_meso       : Compute the sources/sinks for M. rubrum grazing
   !!   p4z_meso_init  : Initialization of M. rubrum grazing the parameters 
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zint          !  interpolation and computation of various fields
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE p4zprod, ONLY: zz_rate_Si_ratio_flag
   USE trdtrc
   USE trd_oce
#if defined key_oxy
  USE p4zoxy
#endif
   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp) ::   zz_rate_mesorub_R              !:
   REAL(wp) ::   zz_rate_mesorub_flagPredSlope  !:
   REAL(wp) ::   zz_rate_mesorub_flagHalfSat    !:
   REAL(wp) ::   zz_rate_mesorub_eff            !:
   REAL(wp) ::   zz_frac_waste_FEC_NH           !:
   REAL(wp) ::   zz_frac_waste_FEC_DON          !waste fraction from mesorub grazing small phyto to DON
   REAL(wp) ::   zz_frac_waste_FEC_PON          !waste fraction from mesorub grazing small phyto to PON
   REAL(wp) ::   zz_frac_waste_FEC_BSi          !waste fraction from mesorub grazing small phyto to Bsi
      
   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmeso.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_meso( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for M. rubrum grazing
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  :: ji, jj, jk, jn
      CHARACTER (len=25) :: charout
      REAL(wp) :: zz_Mesorub_mort_flag, zz_Pflag, zz_Pmyri, zz_was_NH, zz_Mesorub_eat, &
          zz_was_DON, zz_was_PON, zz_was_BSi     
      REAL(wp), POINTER, DIMENSION(:,:,:) :: het_growth_MRub
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  ztrtra
      CHARACTER (len=20) :: strsms
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_meso')
      
      CALL wrk_alloc( jpi, jpj, jpk, het_growth_MRub )
      IF( l_trdtrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtra )  ! store fields before
         ztrtra(:,:,:,:)  = tra(:,:,:,:)
      ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                zz_Pflag = trn(ji,jj,jk,jpphy)
                zz_Pmyri = trn(ji,jj,jk,jpmyr)
                zz_Mesorub_mort_flag = zz_rate_mesorub_R * (zz_Pflag - zz_rate_mesorub_flagPredSlope) &
                    / (zz_rate_mesorub_flagHalfSat + zz_Pflag - zz_rate_mesorub_flagPredSlope &
                    + epsilon(zz_rate_mesorub_flagHalfSat)) &
                    * zz_Pmyri * tgfunc(ji,jj,jk)
                zz_Mesorub_mort_flag = max(zz_Mesorub_mort_flag, 0.d0) 

                zz_was_NH = zz_frac_waste_FEC_NH * zz_Mesorub_mort_flag * &
                     (1-zz_rate_mesorub_eff)
                zz_was_DON = zz_frac_waste_FEC_DON * zz_Mesorub_mort_flag * &
                     (1-zz_rate_mesorub_eff)
                zz_was_PON = zz_frac_waste_FEC_PON * zz_Mesorub_mort_flag * &
                     (1-zz_rate_mesorub_eff)
                zz_was_BSi = zz_frac_waste_FEC_BSi * zz_Mesorub_mort_flag * &
                     zz_rate_Si_ratio_flag ! all Si in grazed phyto directly to BSi but zero for flag
                zz_Mesorub_eat = zz_Mesorub_mort_flag * zz_rate_mesorub_eff
                
               !   Update the arrays TRA which contain the biological sources and sinks
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zz_was_NH
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zz_was_DON
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zz_was_PON
               tra(ji,jj,jk,jpmyr) = tra(ji,jj,jk,jpmyr) + zz_Mesorub_eat
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zz_Mesorub_mort_flag
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zz_was_BSi
               het_growth_MRub(ji,jj,jk) = zz_Mesorub_eat
#if defined key_oxy
               tra(ji,jj,jk,jpo2) = tra(ji,jj,jk,jpo2) - zz_o2_proreg * zz_was_NH
#endif
            END DO
         END DO
      END DO

      IF( l_trdtrc ) THEN
         strsms="MRU_"
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn)-ztrtra(:,:,:,jn), jn, strsms, kt )   ! save trends
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtra ) 
      END IF

      CALL iom_put( "HetMRub", het_growth_MRub)
      CALL wrk_dealloc( jpi, jpj, jpk, het_growth_MRub )
      
      IF( nn_timing == 1 )  CALL timing_stop('p4z_meso')
      !
   END SUBROUTINE p4z_meso

   SUBROUTINE p4z_meso_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of M. rubrum parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismes/ zz_rate_mesorub_R, zz_rate_mesorub_flagPredSlope, &
         & zz_rate_mesorub_flagHalfSat, zz_rate_mesorub_eff, zz_frac_waste_FEC_NH, &
         & zz_frac_waste_FEC_DON, zz_frac_waste_FEC_PON, zz_frac_waste_FEC_BSi 
      INTEGER :: ios           ! Local integer output status for namelist read

      REWIND( numnatp_ref )    ! Namelist nampismes in namelist_ref :  M. rubrum grazing
      READ  ( numnatp_ref, nampismes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in reference namelist', lwp )

      REWIND( numnatp_cfg )    ! Namelist nampismes in namelist_cfg : M. rubrum grazing
      READ  ( numnatp_cfg, nampismes, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismes )


      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for M. rubrum grazing, nampismes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '                zz_rate_mesorub_R =', zz_rate_mesorub_R
         WRITE(numout,*) '    zz_rate_mesorub_flagPredSlope =', zz_rate_mesorub_flagPredSlope
         WRITE(numout,*) '      zz_rate_mesorub_flagHalfSat =', zz_rate_mesorub_flagHalfSat
         WRITE(numout,*) '              zz_rate_mesorub_eff =', zz_rate_mesorub_eff
         WRITE(numout,*) '             zz_frac_waste_FEC_NH =', zz_frac_waste_FEC_NH
         WRITE(numout,*) '            zz_frac_waste_FEC_DON =', zz_frac_waste_FEC_DON
         WRITE(numout,*) '            zz_frac_waste_FEC_PON =', zz_frac_waste_FEC_PON
         WRITE(numout,*) '                FEC  refractory N =', 1.0_wp - &
                 & zz_frac_waste_FEC_NH-zz_frac_waste_FEC_DON-zz_frac_waste_FEC_PON
         WRITE(numout,*) '            zz_frac_waste_FEC_BSi =', zz_frac_waste_FEC_BSi
      ENDIF


   END SUBROUTINE p4z_meso_init


#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_meso                    ! Empty routine
   END SUBROUTINE p4z_meso
#endif 
   !!======================================================================
END MODULE p4zmeso
