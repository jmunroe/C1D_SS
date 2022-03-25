MODULE p4zmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP :   SMELT Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   2015 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_mort       :   Compute the mortality terms for phytoplankton
   !!   p4z_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  SMELT Source Minus Sink variables
   !USE p4zsink         !  vertical flux of particulate matter due to sinking
   !USE p4zprod         !  Primary productivity 
   USE prtctl_trc      !  print control for debugging
   USE p4zprod, ONLY: zz_rate_Si_ratio_diat, zz_rate_Si_ratio_myri, zz_rate_Si_ratio_flag
   USE iom
#if defined key_trdtrc
   USE trdtrc
   USE trd_oce
#endif
#if defined key_oxy
   USE p4zoxy
#endif
   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_mort    
   PUBLIC   p4z_mort_init    

   !! * Shared module variables
   REAL(wp) ::  zz_rate_flag_Rm        ! small phyto natural mort
   REAL(wp) ::  zz_frac_waste_FNM_NH   ! waste frac from small phyto natural mort to NH
   REAL(wp) ::  zz_frac_waste_FNM_DON  ! waste frac from small phyto natural mort to DON
   REAL(wp) ::  zz_frac_waste_FNM_PON  ! waste frac from small phyto natural mort to PON
   REAL(wp) ::  zz_frac_waste_FNM_Bsi  ! waste frac from small phyto natural mort to bSi
   REAL(wp) ::  zz_rate_myri_Rm        ! M rubrum natural mort
   REAL(wp) ::  zz_frac_waste_CNM_NH   ! waste frac from M rubrum natural mort to NH
   REAL(wp) ::  zz_frac_waste_CNM_DON  ! waste frac from M rubrum natural mort to DON
   REAL(wp) ::  zz_frac_waste_CNM_PON  ! waste frac from M rubrum natural mort to PON
   REAL(wp) ::  zz_frac_waste_CNM_Bsi  ! waste frac from M rubrum natural mort to bSi
   REAL(wp) ::  zz_rate_diat_Rm       ! diatom natural mort
   REAL(wp) ::  zz_frac_waste_DNM_NH   ! waste frac from diatom natural mort to NH
   REAL(wp) ::  zz_frac_waste_DNM_DON  ! waste frac from diatom natural mort to DON
   REAL(wp) ::  zz_frac_waste_DNM_PON  ! waste frac from diatom natural mort to PON
   REAL(wp) ::  zz_frac_waste_DNM_Bsi  ! waste frac from diatom natural mort to bSi
      

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER :: jn
#if defined key_trdtrc
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  ztrtra
      CHARACTER (len=20) :: strsms
      !!---------------------------------------------------------------------
      IF( l_trdtrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtra )  ! store fields before
         ztrtra(:,:,:,:)  = tra(:,:,:,:)
      ENDIF
#endif

      ! small phyto
      CALL mort( zz_rate_flag_Rm, zz_frac_waste_FNM_NH, zz_frac_waste_FNM_DON, zz_frac_waste_FNM_PON, &
                  &   zz_frac_waste_FNM_Bsi, zz_rate_Si_ratio_flag, jpphy)
      
      ! M. rubrum
      CALL mort( zz_rate_myri_Rm, zz_frac_waste_CNM_NH, zz_frac_waste_CNM_DON, zz_frac_waste_CNM_PON, &
                  &   zz_frac_waste_CNM_Bsi, zz_rate_Si_ratio_myri, jpmyr)

      ! diatoms
      CALL mort( zz_rate_diat_Rm, zz_frac_waste_DNM_NH, zz_frac_waste_DNM_DON, zz_frac_waste_DNM_PON, &
                  &   zz_frac_waste_DNM_Bsi, zz_rate_Si_ratio_diat, jpdia)

#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         strsms="MRT_"
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn)-ztrtra(:,:,:,jn), jn, strsms, kt )   ! save trends
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtra ) 
      END IF
#endif

   END SUBROUTINE p4z_mort


   SUBROUTINE mort( zz_rate_Rm, zz_frac_waste_iNM_NH, zz_frac_waste_iNM_DON, zz_frac_waste_iNM_PON, &
                  &   zz_frac_waste_iNM_Bsi, zz_rate_Si_ratio, jn)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flag  ***
      !!
      !! ** Purpose :   Compute the mortality terms for small phyto
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      REAL(wp),                               INTENT(in   ) :: zz_rate_Rm 
      REAL(wp),                               INTENT(in   ) :: zz_frac_waste_iNM_NH 
      REAL(wp),                               INTENT(in   ) :: zz_frac_waste_iNM_DON
      REAL(wp),                               INTENT(in   ) :: zz_frac_waste_iNM_PON
      REAL(wp),                               INTENT(in   ) :: zz_frac_waste_iNM_Bsi 
      REAL(wp),                               INTENT(in   ) :: zz_rate_Si_ratio
      INTEGER,                                INTENT(in   ) :: jn
      !REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT(in   ) :: ptn
      !REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT(inout) :: pta

      INTEGER  :: ji, jj, jk
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zz_NatMort(:,:,:)
      !!---------------------------------------------------------------------
      !
      CALL wrk_alloc(jpi, jpj, jpk, zz_NatMort)
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
            
               zz_NatMort(ji,jj,jk) = zz_rate_Rm * tgfunc(ji,jj,jk) * trn(ji,jj,jk,jn) 
               tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) - zz_NatMort(ji,jj,jk)
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zz_frac_waste_iNM_NH * zz_NatMort(ji,jj,jk)
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zz_frac_waste_iNM_DON * zz_NatMort(ji,jj,jk)
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zz_frac_waste_iNM_PON * zz_NatMort(ji,jj,jk)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zz_frac_waste_iNM_Bsi * zz_NatMort(ji,jj,jk) * zz_rate_Si_ratio
#if defined key_oxy
               tra(ji,jj,jk,jpo2) = tra(ji,jj,jk,jpo2) - zz_o2_proreg * zz_frac_waste_iNM_NH * zz_NatMort(ji,jj,jk)
#endif
            END DO
         END DO
      END DO

      IF( jn == jpphy .AND. iom_use("MORTPHY"))  CALL iom_put( "MORTPHY",  zz_NatMort)
      IF( jn == jpmyr .AND. iom_use("MORTMRUB")) CALL iom_put( "MORTMRUB", zz_NatMort)
      IF( jn == jpdia .AND. iom_use("MORTDIAT")) CALL iom_put( "MORTDIAT", zz_NatMort)
      CALL wrk_dealloc(jpi, jpj, jpk, zz_NatMort)
      !
   END SUBROUTINE mort

   SUBROUTINE p4z_mort_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismort/ zz_rate_flag_Rm, zz_frac_waste_FNM_NH, zz_frac_waste_FNM_DON, &
      &  zz_frac_waste_FNM_PON, zz_frac_waste_FNM_Bsi, &
      &                    zz_rate_myri_Rm, zz_frac_waste_CNM_NH, zz_frac_waste_CNM_DON, &
      &  zz_frac_waste_CNM_PON, zz_frac_waste_CNM_Bsi,  &
      &                    zz_rate_diat_Rm, zz_frac_waste_DNM_NH, zz_frac_waste_DNM_DON, &
      &  zz_frac_waste_DNM_PON, zz_frac_waste_DNM_Bsi
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampismort in reference namelist : phytoplankton mortality
      READ  ( numnatp_ref, nampismort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismort in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismort in configuration namelist : phytoplankton mortality
      READ  ( numnatp_cfg, nampismort, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismort in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mort, nampismort'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '  small phyto natural mort                               zz_rate_flag_Rm  =', zz_rate_flag_Rm
         WRITE(numout,*) '  waste frac from small phyto natural mort to NH    zz_frac_waste_FNM_NH  =', zz_frac_waste_FNM_NH
         WRITE(numout,*) '  waste frac from small phyto natural mort to DON   zz_frac_waste_FNM_DON =', zz_frac_waste_FNM_DON
         WRITE(numout,*) '  waste frac from small phyto natural mort to PON   zz_frac_waste_FNM_PON =', zz_frac_waste_FNM_PON
         WRITE(numout,*) '  waste frac small phyto natural mort to refractory zz_frac_waste_FNM_Ref =', &
                 & 1.0_wp-zz_frac_waste_FNM_NH-zz_frac_waste_FNM_DON-zz_frac_waste_FNM_PON
         WRITE(numout,*) '  waste frac from small phyto natural mort to bSi   zz_frac_waste_FNM_Bsi =', zz_frac_waste_FNM_Bsi
         WRITE(numout,*) '  small phyto si/n ratio                            zz_rate_Si_ratio_flag =', zz_rate_Si_ratio_flag
         WRITE(numout,*) '  M rubrum natural mort                               zz_rate_myri_Rm  =', zz_rate_myri_Rm
         WRITE(numout,*) '  waste frac from M rubrum natural mort to NH    zz_frac_waste_CNM_NH  =', zz_frac_waste_CNM_NH
         WRITE(numout,*) '  waste frac from M rubrum natural mort to DON   zz_frac_waste_CNM_DON =', zz_frac_waste_CNM_DON
         WRITE(numout,*) '  waste frac from M rubrum natural mort to PON   zz_frac_waste_CNM_PON =', zz_frac_waste_CNM_PON
         WRITE(numout,*) '  waste frac M rubrum natural mort to refractory zz_frac_waste_CNM_Ref =', &
                 & 1.0_wp-zz_frac_waste_CNM_NH-zz_frac_waste_CNM_DON-zz_frac_waste_CNM_PON
         WRITE(numout,*) '  waste frac from M rubrum natural mort to bSi   zz_frac_waste_CNM_Bsi =', zz_frac_waste_CNM_Bsi
         WRITE(numout,*) '  M rubrum si/n ratio                            zz_rate_Si_ratio_myri =', zz_rate_Si_ratio_myri
         WRITE(numout,*) '  diatom natural mort                           zz_rate_diat_Rm  =', zz_rate_diat_Rm
         WRITE(numout,*) '  waste frac from diatom natural mort to NH   zz_frac_waste_DNM_NH  =', zz_frac_waste_DNM_NH
         WRITE(numout,*) '  waste frac from diatom natural mort to DON  zz_frac_waste_DNM_DON =', zz_frac_waste_DNM_DON
         WRITE(numout,*) '  waste frac from diatom natural mort to PON  zz_frac_waste_DNM_PON =', zz_frac_waste_DNM_PON
         WRITE(numout,*) '  waste frac diatom natural mort to refractory zz_frac_waste_DNM_Ref=', &
                 & 1.0_wp-zz_frac_waste_DNM_NH-zz_frac_waste_DNM_DON-zz_frac_waste_DNM_PON
         WRITE(numout,*) '  waste frac from diatom natural mort to bSi  zz_frac_waste_DNM_Bsi =', zz_frac_waste_DNM_Bsi
         WRITE(numout,*) '  diatom si/n ratio                        zz_rate_Si_ratio_diat =', zz_rate_Si_ratio_diat
      ENDIF

   END SUBROUTINE p4z_mort_init

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_mort                    ! Empty routine
   END SUBROUTINE p4z_mort
#endif 
   !!======================================================================
END MODULE p4zmort
