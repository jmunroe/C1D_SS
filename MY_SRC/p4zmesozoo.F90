MODULE p4zmesozoo
   !!======================================================================
   !!                         ***  MODULE p4zmesozoo  ***
   !! TOP :   SMELT Compute the sources/sinks for mesozoozooplankton
   !!======================================================================
   !! History :   2015 (E. Olson) adapted from SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                           SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_mesozoo       :   Compute the sources/sinks for mesozoozooplankton
   !!   p4z_mesozoo_init  :   Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zlim          !  Co-limitations
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zint          !  interpolation and computation of various fields
   USE p4zprod         !  production
   USE prtctl_trc      !  print control for debugging
   USE dom_oce        ! ocean space and time domain
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

   PUBLIC   p4z_mesozoo         ! called in p4zbio.F90
   PUBLIC   p4z_mesozoo_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
      REAL(wp) ::   zz_rate_mesozoo_winterconc       !uM N mesozooplankton background concentration
      REAL(wp), dimension (1:3) :: zz_rate_mesozoo_sumpeakval !uM N magnitude of mesozooplankton summer concentration peaks
      REAL(wp), dimension (1:3) :: zz_rate_mesozoo_sumpeakpos ! year-day times of mesozooplankton summer concentration peaks
      REAL(wp), dimension (1:3) :: zz_rate_mesozoo_sumpeakwid ! year-days widths of mesozooplankton summer concentration peaks
      REAL(wp) ::   zz_rate_mesozoo_alpha            ! 0 to 1. alpha = 0: spatially uniform mesozoo; alpha = 1: mesozoo proportional to prey
      REAL(wp) ::   zz_rate_mesozoo_R                ! uM N / s mesozooplankton maximum ingestion rate @ 20 deg C
      REAL(wp) ::   zz_rate_mesozoo_Rm               ! uM N / s mesozooplankton natural mortality rate @ 20 deg C
      REAL(wp) ::   zz_rate_mesozoo_excr             ! uM N / s mesozooplankton excretion rate @ 20 deg C
      REAL(wp) ::   zz_rate_mesozoo_PredSlope        ! uM N mesozooplankton total grazing limit
      REAL(wp) ::   zz_rate_mesozoo_HalfSat          ! uM N mesozooplankton total grazing half saturation constant
      REAL(wp) ::   zz_rate_mesozoo_diatPref        !      mesozooplankton preference for diatoms
      REAL(wp) ::   zz_rate_mesozoo_diatPredslope   ! uM N mesozooplankton diatom grazing limit
      REAL(wp) ::   zz_rate_mesozoo_diatHalfSat     ! uM N mesozooplankton diatom grazing half saturation constant
      REAL(wp) ::   zz_rate_mesozoo_myriPref         !      mesozooplankton preference for M rubrum
      REAL(wp) ::   zz_rate_mesozoo_myriPredslope    ! uM N mesozooplankton M rubrum grazing limit
      REAL(wp) ::   zz_rate_mesozoo_myriHalfSat      ! uM N mesozooplankton M rubrum grazing half saturation const
      REAL(wp) ::   zz_rate_mesozoo_flagPref         !      mesozooplankton preference for small phyto
      REAL(wp) ::   zz_rate_mesozoo_flagPredslope    ! uM N mesozooplankton small phyto grazing limit
      REAL(wp) ::   zz_rate_mesozoo_flagHalfSat      ! uM N mesozooplankton small phyto grazing half saturation const
      REAL(wp) ::   zz_rate_mesozoo_PON_Pref         !      mesozooplankton preference for PON
      REAL(wp) ::   zz_rate_mesozoo_PON_Predslope    ! uM N mesozooplankton PON feeding limit
      REAL(wp) ::   zz_rate_mesozoo_PON_HalfSat      ! uM N mesozooplankton PON feeding half saturation const
      REAL(wp) ::   zz_rate_mesozoo_Z_Pref           !      mesozooplankton preference for microzooplankton
      REAL(wp) ::   zz_rate_mesozoo_Z_Predslope      ! uM N mesozooplankton microzooplankton grazing limit
      REAL(wp) ::   zz_rate_mesozoo_Z_HalfSat        ! uM N mesozooplankton microzoo grazing half saturation const
      REAL(wp) ::   zz_rate_mesozoo_eff              ! mesozooplankton grazed mass assimilation efficiency
      REAL(wp) ::   zz_frac_waste_MNM_NH             ! waste fraction from meso-zoo natural mortality to NH
      REAL(wp) ::   zz_frac_waste_MNM_DON            ! waste fraction from meso-zoo natural mortality to DON
      REAL(wp) ::   zz_frac_waste_MNM_PON            ! waste fraction from meso-zoo natural mortality to PON
      REAL(wp) ::   zz_frac_waste_MNM_BSi            ! waste fraction from meso-zoo natural mortality to bSi
      REAL(wp) ::   zz_frac_waste_MEX_NH             ! waste fraction from meso-zoo excretion to NH
      REAL(wp) ::   zz_frac_waste_MEX_DON            ! waste fraction from meso-zoo excretion to DON
      REAL(wp) ::   zz_frac_waste_MEX_PON            ! waste fraction from meso-zoo excretion to PON
      REAL(wp) ::   zz_frac_waste_MEX_Bsi            ! waste fraction from meso-zoo excretion to bSi
      REAL(wp) ::   zz_frac_waste_PEM_NH             ! waste fraction from mesozoo grazing PON to NH
      REAL(wp) ::   zz_frac_waste_PEM_DON            ! waste fraction from mesozoo grazing PON to DON
      REAL(wp) ::   zz_frac_waste_PEM_PON            ! waste fraction from mesozoo grazing PON to PON
      REAL(wp) ::   zz_frac_waste_PEM_BSi            ! waste fraction from mesozoo grazing PON to Bsi
      REAL(wp) ::   zz_frac_waste_DEM_NH             ! waste fraction from mesozoo grazing diatom to NH #  Alain, excretion rate to growth rate
      REAL(wp) ::   zz_frac_waste_DEM_DON            ! waste fraction from mesozoo grazing diatom to DON# Denman and Pena about 1/2 to dissolved
      REAL(wp) ::   zz_frac_waste_DEM_PON            ! waste fraction from mesozoo grazing diatom to PON
      REAL(wp) ::   zz_frac_waste_DEM_BSi            ! waste fraction from mesozoo grazing diatom to Bsi
      REAL(wp) ::   zz_frac_waste_FEM_NH             !  waste fraction from mesozoo grazing small phyto to NH
      REAL(wp) ::   zz_frac_waste_FEM_DON            ! waste fraction from mesozoo grazing small phyto to DON
      REAL(wp) ::   zz_frac_waste_FEM_PON            ! waste fraction from mesozoo grazing small phyto to PON
      REAL(wp) ::   zz_frac_waste_FEM_BSi            ! waste fraction from mesozoo grazing small phyto to Bsi
      REAL(wp) ::   zz_frac_waste_ZEM_NH             ! waste fraction from mesozoo grazing microzoo to NH
      REAL(wp) ::   zz_frac_waste_ZEM_DON            ! waste fraction from mesozoo grazing microzoo to DON
      REAL(wp) ::   zz_frac_waste_ZEM_PON            ! waste fraction from mesozoo grazing microzoo to PON
      REAL(wp) ::   zz_frac_waste_ZEM_BSi            ! waste fraction from mesozoo grazing microzoo to Bsi
      REAL(wp) ::   zz_frac_waste_CEM_NH             ! waste fraction from mesozoo grazing M rubrum to NH
      REAL(wp) ::   zz_frac_waste_CEM_DON            ! waste fraction from mesozoo grazing M rubrum to DON
      REAL(wp) ::   zz_frac_waste_CEM_PON            ! waste fraction from mesozoo grazing M rubrum to PON
      REAL(wp) ::   zz_frac_waste_CEM_BSi            ! waste fraction from mesozoo grazing M rubrum to bSi
      INTEGER  ::   jpk40
      REAL(wp) ::   zz_average_prey

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmesozoo.F90 3830 2013-03-06 11:00:26Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_mesozoo( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mesozoo  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozoozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step, substep
      INTEGER  :: ji, jj, jk, jn

      REAL(wp) ::   zz_day, zz_NatMort_mesozoo, zz_Excr_mesozoo, zz_food_limitation, zz_denominator, zz_MesZoBar
      
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zz_Pdiat, zz_D_PON, zz_Pmyri, zz_Pflag, zz_Z, zz_Mesozoo, &
                  zz_Meso_mort_diat, zz_Meso_mort_myri, zz_Meso_mort_flag, zz_Meso_graz_PON, zz_Meso_mort_Z, zz_masksum
#if defined key_trdtrc
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  ztrtra
      CHARACTER (len=20) :: strsms
#endif
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_mesozoo')

      CALL wrk_alloc( jpi, jpj, jpk, zz_Pdiat, zz_D_PON, zz_Pmyri, zz_Pflag, zz_Z, zz_Mesozoo, zz_masksum )
      CALL wrk_alloc( jpi, jpj, jpk, zz_Meso_mort_diat, zz_Meso_mort_myri, zz_Meso_mort_flag, zz_Meso_graz_PON, zz_Meso_mort_Z )

#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtra )  ! store fields before
         ztrtra(:,:,:,:)  = tra(:,:,:,:)
      ENDIF
#endif
      
      zz_Pdiat(:,:,:)  = trn(:,:,:,jpdia)
      zz_Pmyri(:,:,:)  = trn(:,:,:,jpmyr)
      zz_Pflag(:,:,:) = trn(:,:,:,jpphy)
      zz_Z(:,:,:) = trn(:,:,:,jpzoo)
      zz_D_PON(:,:,:) = trn(:,:,:,jppon)
      zz_day=nday_year

#if defined key_agrif
      IF( Agrif_Root() ) THEN
#endif

      ! assuming average prey is global average over upper 40m of ocean (not topography)
      ! thus also assuming that mesozoo summer and winter concs are averages over upper 40m of ocean
      zz_masksum(:,:,:) = 0.0_wp
      zz_masksum(:,:,1:jpk40) = cvol(:,:,1:jpk40)

      zz_average_prey = glob_sum((zz_Pdiat + zz_D_PON + zz_Pmyri + zz_Pflag + zz_Z)*zz_masksum) &
                        / glob_sum(zz_masksum)

#if defined key_agrif
      ELSE
        zz_average_prey = Agrif_Parent(zz_average_prey)
      ENDIF
#endif

      zz_MesZoBar = zz_rate_mesozoo_winterconc + &
                    ( sum ( zz_rate_mesozoo_sumpeakval * &
                    exp( -(zz_day-zz_rate_mesozoo_sumpeakpos)**2 / zz_rate_mesozoo_sumpeakwid**2 ) ) &
                    + sum ( zz_rate_mesozoo_sumpeakval * &
                    exp( -(zz_day-zz_rate_mesozoo_sumpeakpos-365.25_wp)**2 / zz_rate_mesozoo_sumpeakwid**2 ) ) &
                    + sum ( zz_rate_mesozoo_sumpeakval * &
                    exp( -(zz_day-zz_rate_mesozoo_sumpeakpos+365.25_wp)**2 / zz_rate_mesozoo_sumpeakwid**2 ) ) )
      zz_Mesozoo(:,:,:) = zz_MesZoBar &
                    * zz_rate_mesozoo_alpha * (zz_Pdiat(:,:,:) + zz_D_PON(:,:,:) + zz_Pmyri(:,:,:) + zz_Pflag(:,:,:) + zz_Z(:,:,:)) &
                     / ( zz_average_prey + epsilon(zz_average_prey)) + zz_MesZoBar * (1.00_wp - zz_rate_mesozoo_alpha)    
      IF( lk_iomput ) THEN
           IF( knt == nrdttrc ) THEN
              CALL iom_put( "MESZ", zz_Mesozoo)  ! euphotic layer deptht
           ENDIF
      ENDIF
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

                ! mesozoo natural mortality:
                zz_NatMort_mesozoo = zz_rate_mesozoo_Rm * tgfunc(ji,jj,jk) * zz_Mesozoo(ji,jj,jk)
                zz_Excr_mesozoo = zz_rate_mesozoo_excr * tgfunc(ji,jj,jk) * zz_Mesozoo(ji,jj,jk)

                tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + (zz_frac_waste_MNM_NH * zz_NatMort_mesozoo &
                     + zz_frac_waste_MEX_NH * zz_Excr_mesozoo) 
                tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + (zz_frac_waste_MNM_DON * zz_NatMort_mesozoo &
                     + zz_frac_waste_MEX_DON * zz_Excr_mesozoo) 
                tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + (zz_frac_waste_MNM_PON * zz_NatMort_mesozoo &
                     + zz_frac_waste_MEX_PON * zz_Excr_mesozoo) 
                ! No transfer to BSi
#if defined key_oxy
                tra(ji,jj,jk,jpo2) = tra(ji,jj,jk,jpo2) - zz_o2_proreg * (zz_frac_waste_MNM_NH * zz_NatMort_mesozoo &
                     + zz_frac_waste_MEX_NH * zz_Excr_mesozoo) 
#endif

                ! global food limitation
                zz_food_limitation = (zz_Pdiat(ji,jj,jk) + zz_D_PON(ji,jj,jk) + zz_Pmyri(ji,jj,jk) + zz_Pflag(ji,jj,jk) &
                    + zz_Z(ji,jj,jk) - zz_rate_mesozoo_PredSlope) / &
                    (zz_rate_mesozoo_HalfSat + zz_Pdiat(ji,jj,jk) + zz_D_PON(ji,jj,jk) + zz_Pmyri(ji,jj,jk) &
                    + zz_Pflag(ji,jj,jk) + zz_Z(ji,jj,jk) - zz_rate_mesozoo_PredSlope &
                    + epsilon(zz_rate_mesozoo_HalfSat))

                zz_denominator = (zz_rate_mesozoo_diatPref * zz_Pdiat(ji,jj,jk) + &
                    zz_rate_mesozoo_myriPref * zz_Pmyri(ji,jj,jk) + &
                    zz_rate_mesozoo_flagPref * zz_Pflag(ji,jj,jk) + &
                    zz_rate_mesozoo_PON_Pref * zz_D_PON(ji,jj,jk) + &
                    zz_rate_mesozoo_Z_Pref * zz_Z(ji,jj,jk) + epsilon(zz_Pdiat(ji,jj,jk)) )

                ! limitation based on diatoms
                zz_Meso_mort_diat(ji,jj,jk) = min(zz_rate_mesozoo_diatPref * zz_food_limitation &
                    * zz_Pdiat(ji,jj,jk) / zz_denominator, &
                    (zz_Pdiat(ji,jj,jk) - zz_rate_mesozoo_diatPredslope) / &
                    (zz_rate_mesozoo_diatHalfSat + zz_Pdiat(ji,jj,jk) &
                    - zz_rate_mesozoo_diatPredSlope + &
                    epsilon(zz_rate_mesozoo_diatHalfSat)))

                ! limitation based on M rubrum
                zz_Meso_mort_myri(ji,jj,jk) = min(zz_rate_mesozoo_myriPref * zz_food_limitation &
                    * zz_Pmyri(ji,jj,jk) / zz_denominator, &
                    (zz_Pmyri(ji,jj,jk) - zz_rate_mesozoo_myriPredslope) / &
                    (zz_rate_mesozoo_myriHalfSat + zz_Pmyri(ji,jj,jk) &
                    - zz_rate_mesozoo_myriPredSlope + &
                    epsilon(zz_rate_mesozoo_myriHalfSat)))

                ! limitation based on flagplankton
                zz_Meso_mort_flag(ji,jj,jk) = min(zz_rate_mesozoo_flagPref * zz_food_limitation &
                    * zz_Pflag(ji,jj,jk) / zz_denominator, &
                    (zz_Pflag(ji,jj,jk) - zz_rate_mesozoo_flagPredslope) / &
                    (zz_rate_mesozoo_flagHalfSat + zz_Pflag(ji,jj,jk) &
                    - zz_rate_mesozoo_flagPredSlope + &
                    epsilon(zz_rate_mesozoo_flagHalfSat)))

                ! limitation based on PON
                zz_Meso_graz_PON(ji,jj,jk) = min(zz_rate_mesozoo_PON_Pref * zz_food_limitation &
                    * zz_D_PON(ji,jj,jk) / zz_denominator, &
                    (zz_D_PON(ji,jj,jk) - zz_rate_mesozoo_PON_Predslope) / &
                    (zz_rate_mesozoo_PON_HalfSat + zz_D_PON(ji,jj,jk) - zz_rate_mesozoo_PON_PredSlope &
                    + epsilon(zz_rate_mesozoo_PON_HalfSat)))

                ! limitation based on Z
                zz_Meso_mort_Z(ji,jj,jk) = min(zz_rate_mesozoo_Z_Pref * zz_food_limitation &
                    * zz_Z(ji,jj,jk) / zz_denominator, &
                    (zz_Z(ji,jj,jk) - zz_rate_mesozoo_Z_Predslope) / &
                    (zz_rate_mesozoo_Z_HalfSat + zz_Z(ji,jj,jk) - zz_rate_mesozoo_Z_PredSlope &
                    + epsilon(zz_rate_mesozoo_Z_HalfSat)))

                ! global corrected by individual
                zz_food_limitation = zz_Meso_mort_diat(ji,jj,jk) + zz_Meso_mort_myri(ji,jj,jk) &
                    + zz_Meso_mort_flag(ji,jj,jk) + zz_Meso_graz_PON(ji,jj,jk) + zz_Meso_mort_Z(ji,jj,jk)

                zz_Meso_mort_diat(ji,jj,jk) = zz_rate_mesozoo_R * tgfunc(ji,jj,jk) * zz_Mesozoo(ji,jj,jk) * &
                    max(0.d0, zz_Meso_mort_diat(ji,jj,jk))

                zz_Meso_mort_myri(ji,jj,jk) = zz_rate_mesozoo_R * tgfunc(ji,jj,jk) * zz_Mesozoo(ji,jj,jk) * &
                    max(0.d0, zz_Meso_mort_myri(ji,jj,jk))

                zz_Meso_mort_flag(ji,jj,jk) = zz_rate_mesozoo_R * tgfunc(ji,jj,jk) * zz_Mesozoo(ji,jj,jk) * &
                    max(0.d0, zz_Meso_mort_flag(ji,jj,jk))

                zz_Meso_graz_PON(ji,jj,jk) = zz_rate_mesozoo_R * tgfunc(ji,jj,jk) * zz_Mesozoo(ji,jj,jk) * &
                    max(0.d0, zz_Meso_graz_PON(ji,jj,jk))

                zz_Meso_mort_Z(ji,jj,jk) = zz_rate_mesozoo_R * tgfunc(ji,jj,jk) * zz_Mesozoo(ji,jj,jk) * &
                    max(0.d0, zz_Meso_mort_Z(ji,jj,jk))

                tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zz_Meso_mort_diat(ji,jj,jk) 
                tra(ji,jj,jk,jpmyr) = tra(ji,jj,jk,jpmyr) - zz_Meso_mort_myri(ji,jj,jk) 
                tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zz_Meso_mort_flag(ji,jj,jk) 
                tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zz_Meso_mort_Z(ji,jj,jk) 
                
                tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + (zz_frac_waste_PEM_NH * zz_Meso_graz_PON(ji,jj,jk) +     &
                     zz_frac_waste_DEM_NH  * zz_Meso_mort_diat(ji,jj,jk) +               &
                     zz_frac_waste_CEM_NH  * zz_Meso_mort_myri(ji,jj,jk) +                &
                     zz_frac_waste_FEM_NH  * zz_Meso_mort_flag(ji,jj,jk) +                &
                     zz_frac_waste_ZEM_NH  * zz_Meso_mort_Z(ji,jj,jk)) * (1.0_wp-zz_rate_mesozoo_eff) 
                tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + (zz_frac_waste_PEM_DON * zz_Meso_graz_PON(ji,jj,jk) +  &
                     zz_frac_waste_DEM_DON * zz_Meso_mort_diat(ji,jj,jk) +               &
                     zz_frac_waste_CEM_DON * zz_Meso_mort_myri(ji,jj,jk) +                &
                     zz_frac_waste_FEM_DON * zz_Meso_mort_flag(ji,jj,jk) +                &
                     zz_frac_waste_ZEM_DON * zz_Meso_mort_Z(ji,jj,jk)) * (1.0_wp-zz_rate_mesozoo_eff) 
                tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + ((zz_frac_waste_PEM_PON * zz_Meso_graz_PON(ji,jj,jk) +  &
                     zz_frac_waste_DEM_PON * zz_Meso_mort_diat(ji,jj,jk) +               &
                     zz_frac_waste_CEM_PON * zz_Meso_mort_myri(ji,jj,jk) +                &
                     zz_frac_waste_FEM_PON * zz_Meso_mort_flag(ji,jj,jk) +                &
                     zz_frac_waste_ZEM_PON * zz_Meso_mort_Z(ji,jj,jk)) * (1.0_wp-zz_rate_mesozoo_eff) -             &
                     zz_Meso_graz_PON(ji,jj,jk)) 
                tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + ( &
                     zz_frac_waste_DEM_BSi * zz_Meso_mort_diat(ji,jj,jk) * zz_rate_Si_ratio_diat + &
                     zz_frac_waste_CEM_BSi * zz_Meso_mort_myri(ji,jj,jk) * zz_rate_Si_ratio_myri +   &
                     zz_frac_waste_FEM_BSi * zz_Meso_mort_flag(ji,jj,jk) * zz_rate_Si_ratio_flag   &
                     )   ! assuming Si not taken up by mesozoo;
#if defined key_oxy
                tra(ji,jj,jk,jpo2) = tra(ji,jj,jk,jpo2) - zz_o2_proreg * ( (zz_frac_waste_PEM_NH * zz_Meso_graz_PON(ji,jj,jk) +     &
                     zz_frac_waste_DEM_NH  * zz_Meso_mort_diat(ji,jj,jk) +               &
                     zz_frac_waste_CEM_NH  * zz_Meso_mort_myri(ji,jj,jk) +                &
                     zz_frac_waste_FEM_NH  * zz_Meso_mort_flag(ji,jj,jk) +                &
                     zz_frac_waste_ZEM_NH  * zz_Meso_mort_Z(ji,jj,jk)) * (1.0_wp-zz_rate_mesozoo_eff) )
#endif

            END DO
         END DO
      END DO
     
      CALL iom_put( "GRMESZDIAT", zz_Meso_mort_diat)
      CALL iom_put( "GRMESZMRUB", zz_Meso_mort_myri)
      CALL iom_put( "GRMESZPHY", zz_Meso_mort_flag)
      CALL iom_put( "GRMESZPON", zz_Meso_graz_PON)
      CALL iom_put( "GRMESZMICZ", zz_Meso_mort_Z)

#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         strsms="MSZ_"
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn)-ztrtra(:,:,:,jn), jn, strsms, kt )   ! save trends
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtra ) 
      END IF
#endif
 
      CALL wrk_dealloc( jpi, jpj, jpk, zz_Pdiat, zz_D_PON, zz_Pmyri, zz_Pflag, zz_Z, zz_Mesozoo, zz_masksum )
      CALL wrk_dealloc( jpi, jpj, jpk, zz_Meso_mort_diat, zz_Meso_mort_myri, zz_Meso_mort_flag, zz_Meso_graz_PON, zz_Meso_mort_Z )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_mesozoo')
      !
   END SUBROUTINE p4z_mesozoo


   SUBROUTINE p4z_mesozoo_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mesozoo_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismezo
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismezo/ zz_rate_mesozoo_winterconc, zz_rate_mesozoo_sumpeakval, &
        &  zz_rate_mesozoo_sumpeakpos, zz_rate_mesozoo_sumpeakwid, zz_rate_mesozoo_alpha, zz_rate_mesozoo_R, zz_rate_mesozoo_Rm, &
        &  zz_rate_mesozoo_excr, zz_rate_mesozoo_PredSlope, zz_rate_mesozoo_HalfSat, zz_rate_mesozoo_diatPref, &
        &  zz_rate_mesozoo_diatPredslope, zz_rate_mesozoo_diatHalfSat, zz_rate_mesozoo_myriPref, &
        &  zz_rate_mesozoo_myriPredslope, zz_rate_mesozoo_myriHalfSat, zz_rate_mesozoo_flagPref, &
        &  zz_rate_mesozoo_flagPredslope, zz_rate_mesozoo_flagHalfSat, zz_rate_mesozoo_PON_Pref, &
        &  zz_rate_mesozoo_PON_Predslope, zz_rate_mesozoo_PON_HalfSat, zz_rate_mesozoo_Z_Pref, &
        &  zz_rate_mesozoo_Z_Predslope, zz_rate_mesozoo_Z_HalfSat, zz_rate_mesozoo_eff, zz_frac_waste_MNM_NH, zz_frac_waste_MNM_DON, &
        &  zz_frac_waste_MNM_PON, zz_frac_waste_MNM_BSi, zz_frac_waste_MEX_NH, zz_frac_waste_MEX_DON, &
        &  zz_frac_waste_MEX_PON, zz_frac_waste_MEX_Bsi, zz_frac_waste_PEM_NH, zz_frac_waste_PEM_DON, &
        &  zz_frac_waste_PEM_PON, zz_frac_waste_PEM_BSi, zz_frac_waste_DEM_NH, zz_frac_waste_DEM_DON, &
        &  zz_frac_waste_DEM_PON, zz_frac_waste_DEM_BSi, zz_frac_waste_FEM_NH, zz_frac_waste_FEM_DON, &
        &  zz_frac_waste_FEM_PON, zz_frac_waste_FEM_BSi, zz_frac_waste_ZEM_NH, zz_frac_waste_ZEM_DON, &
        &  zz_frac_waste_ZEM_PON, zz_frac_waste_ZEM_BSi, zz_frac_waste_CEM_NH, zz_frac_waste_CEM_DON, &
        &  zz_frac_waste_CEM_PON, zz_frac_waste_CEM_BSi 
      INTEGER :: ios, jk                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampismezo in reference namelist : mesozooplankton parameters
      READ  ( numnatp_ref, nampismezo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismezo in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismezo in configuration namelist : mesozooplankton parameters
      READ  ( numnatp_cfg, nampismezo, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismezo in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismezo )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, nampismezo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) 'uM N mesozoo background conc                zz_rate_mesozoo_winterconc=', zz_rate_mesozoo_winterconc
         WRITE(numout,*) 'uM N magnitude of mesozoo summer conc peaks zz_rate_mesozoo_sumpeakval=', zz_rate_mesozoo_sumpeakval
         WRITE(numout,*) 'year-day times of mesozoo summer conc peaks zz_rate_mesozoo_sumpeakpos=', zz_rate_mesozoo_sumpeakpos
         WRITE(numout,*) 'year-day width of mesozoo summer conc peaks zz_rate_mesozoo_sumpeakwid=', zz_rate_mesozoo_sumpeakwid
         WRITE(numout,*) 'uM N/s mesozoo max ingestion @ 20 deg C     zz_rate_mesozoo_R=', zz_rate_mesozoo_R
         WRITE(numout,*) 'uM N/s mesozoo natural mortality @ 20 deg C zz_rate_mesozoo_Rm=', zz_rate_mesozoo_Rm
         WRITE(numout,*) 'uM N/s mesozoo excretion rate @ 20 deg C    zz_rate_mesozoo_excr=', zz_rate_mesozoo_excr
         WRITE(numout,*) 'uM N mesozoo total grazing limit            zz_rate_mesozoo_PredSlope=', zz_rate_mesozoo_PredSlope
         WRITE(numout,*) 'uM N mesozoo total grazing half saturation  zz_rate_mesozoo_HalfSat=', zz_rate_mesozoo_HalfSat
         WRITE(numout,*) 'mesozoo preference for diatoms              zz_rate_mesozoo_diatPref=', zz_rate_mesozoo_diatPref
         WRITE(numout,*) 'uM N mesozoo diatom grazing limit           zz_rate_mesozoo_diatPredslope=', zz_rate_mesozoo_diatPredslope
         WRITE(numout,*) 'uM N mesozoo diatom grazing half saturation zz_rate_mesozoo_diatHalfSat=', zz_rate_mesozoo_diatHalfSat
         WRITE(numout,*) 'mesozoo preference for M rubrum           zz_rate_mesozoo_myriPref=', zz_rate_mesozoo_myriPref
         WRITE(numout,*) 'uM N mesozoo M rubrum grazing limit       zz_rate_mesozoo_myriPredslope=', zz_rate_mesozoo_myriPredslope
         WRITE(numout,*) 'uM N mesozoo M rubrum grazing half saturation const zz_rate_mesozoo_myriHalfSat=', zz_rate_mesozoo_myriHalfSat
         WRITE(numout,*) 'mesozoo preference for small phyto           zz_rate_mesozoo_flagPref=', zz_rate_mesozoo_flagPref
         WRITE(numout,*) 'uM N mesozoo small phyto grazing limit       zz_rate_mesozoo_flagPredslope=', zz_rate_mesozoo_flagPredslope
         WRITE(numout,*) 'uM N mesozoo small phyto grazing half saturation const zz_rate_mesozoo_flagHalfSat=', zz_rate_mesozoo_flagHalfSat
         WRITE(numout,*) 'mesozoo preference for PON                  zz_rate_mesozoo_PON_Pref=', zz_rate_mesozoo_PON_Pref
         WRITE(numout,*) 'uM N mesozoo PON feeding limit              zz_rate_mesozoo_PON_Predslope=', zz_rate_mesozoo_PON_Predslope
         WRITE(numout,*) 'uM N mesozoo PON feeding half saturation const zz_rate_mesozoo_PON_HalfSat=', zz_rate_mesozoo_PON_HalfSat
         WRITE(numout,*) 'mesozoo preference for microzoo             zz_rate_mesozoo_Z_Pref=', zz_rate_mesozoo_Z_Pref
         WRITE(numout,*) 'uM N mesozoo microzoo grazing limit         zz_rate_mesozoo_Z_Predslope=', zz_rate_mesozoo_Z_Predslope
         WRITE(numout,*) 'uM N mesozoo microzoo grazing half saturation const zz_rate_mesozoo_Z_HalfSat=', zz_rate_mesozoo_Z_HalfSat
         WRITE(numout,*) 'mesozo grazed mass assimilation efficiency  zz_rate_mesozoo_eff=', zz_rate_mesozoo_eff
         WRITE(numout,*) 'waste fraction from meso-zoo natural mortality to NH   zz_frac_waste_MNM_NH =', zz_frac_waste_MNM_NH
         WRITE(numout,*) 'waste fraction from meso-zoo natural mortality to DON  zz_frac_waste_MNM_DON=', zz_frac_waste_MNM_DON
         WRITE(numout,*) 'waste fraction from meso-zoo natural mortality to PON  zz_frac_waste_MNM_PON=', zz_frac_waste_MNM_PON
         WRITE(numout,*) 'waste fraction from meso-zoo natural mortality to bSi  zz_frac_waste_MNM_BSi=', zz_frac_waste_MNM_BSi
         WRITE(numout,*) 'waste fraction from meso-zoo excretion to NH     zz_frac_waste_MEX_NH =', zz_frac_waste_MEX_NH
         WRITE(numout,*) 'waste fraction from meso-zoo excretion to DON    zz_frac_waste_MEX_DON=', zz_frac_waste_MEX_DON
         WRITE(numout,*) 'waste fraction from meso-zoo excretion to PON    zz_frac_waste_MEX_PON=', zz_frac_waste_MEX_PON
         WRITE(numout,*) 'waste fraction from meso-zoo excretion to bSi    zz_frac_waste_MEX_Bsi=', zz_frac_waste_MEX_Bsi
         WRITE(numout,*) 'waste fraction from mesozoo grazing PON to NH    zz_frac_waste_PEM_NH =', zz_frac_waste_PEM_NH
         WRITE(numout,*) 'waste fraction from mesozoo grazing PON to DON   zz_frac_waste_PEM_DON=', zz_frac_waste_PEM_DON
         WRITE(numout,*) 'waste fraction from mesozoo grazing PON to PON   zz_frac_waste_PEM_PON=', zz_frac_waste_PEM_PON
         WRITE(numout,*) 'waste fraction from mesozoo grazing PON to Bsi   zz_frac_waste_PEM_BSi=', zz_frac_waste_PEM_BSi
         WRITE(numout,*) 'waste fraction from mesozoo grazing diatom to NH   zz_frac_waste_DEM_NH =',  zz_frac_waste_DEM_NH
         WRITE(numout,*) 'waste fraction from mesozoo grazing diatom to DON  zz_frac_waste_DEM_DON=',  zz_frac_waste_DEM_DON
         WRITE(numout,*) 'waste fraction from mesozoo grazing diatom to PON  zz_frac_waste_DEM_PON=',  zz_frac_waste_DEM_PON
         WRITE(numout,*) 'waste fraction from mesozoo grazing diatom to Bsi  zz_frac_waste_DEM_BSi=',  zz_frac_waste_DEM_BSi
         WRITE(numout,*) 'waste fraction from mesozoo grazing small phyto to NH    zz_frac_waste_FEM_NH =',  zz_frac_waste_FEM_NH
         WRITE(numout,*) 'waste fraction from mesozoo grazing small phyto to DON   zz_frac_waste_FEM_DON=', zz_frac_waste_FEM_DON
         WRITE(numout,*) 'waste fraction from mesozoo grazing small phyto to PON   zz_frac_waste_FEM_PON=', zz_frac_waste_FEM_PON
         WRITE(numout,*) 'waste fraction from mesozoo grazing small phyto to Bsi   zz_frac_waste_FEM_BSi=', zz_frac_waste_FEM_BSi
         WRITE(numout,*) 'waste fraction from mesozoo grazing microzoo to NH     zz_frac_waste_ZEM_NH =', zz_frac_waste_ZEM_NH
         WRITE(numout,*) 'waste fraction from mesozoo grazing microzoo to DON    zz_frac_waste_ZEM_DON=', zz_frac_waste_ZEM_DON
         WRITE(numout,*) 'waste fraction from mesozoo grazing microzoo to PON    zz_frac_waste_ZEM_PON=', zz_frac_waste_ZEM_PON
         WRITE(numout,*) 'waste fraction from mesozoo grazing microzoo to Bsi    zz_frac_waste_ZEM_BSi=', zz_frac_waste_ZEM_BSi
         WRITE(numout,*) 'waste fraction from mesozoo grazing M rubrum to NH    zz_frac_waste_CEM_NH =',  zz_frac_waste_CEM_NH
         WRITE(numout,*) 'waste fraction from mesozoo grazing M rubrum to DON   zz_frac_waste_CEM_DON=', zz_frac_waste_CEM_DON
         WRITE(numout,*) 'waste fraction from mesozoo grazing M rubrum to PON   zz_frac_waste_CEM_PON=', zz_frac_waste_CEM_PON
         WRITE(numout,*) 'waste fraction from mesozoo grazing M rubrum to bSi   zz_frac_waste_CEM_BSi=', zz_frac_waste_CEM_BSi
         WRITE(numout,*) 'diatom silicon/nitrogen ratio          zz_rate_Si_ratio_diat=', zz_rate_Si_ratio_diat
         WRITE(numout,*) 'M rubrum silicon/nitrogen ratio      zz_rate_Si_ratio_myri=', zz_rate_Si_ratio_myri
         WRITE(numout,*) 'small phyto silicon/nitrogen ratio      zz_rate_Si_ratio_flag=', zz_rate_Si_ratio_flag
         WRITE(numout,*) '                       FEM refractory N =', &
                 & 1.0_wp - zz_frac_waste_FEM_NH - zz_frac_waste_FEM_DON - zz_frac_waste_FEM_PON
         WRITE(numout,*) '                       CEM refractory N =', &
                 & 1.0_wp - zz_frac_waste_CEM_NH - zz_frac_waste_CEM_DON - zz_frac_waste_CEM_PON
         WRITE(numout,*) '                       DEM refractory N =', &
                 & 1.0_wp - zz_frac_waste_DEM_NH - zz_frac_waste_DEM_DON - zz_frac_waste_DEM_PON
         WRITE(numout,*) '                       PEM refractory N =', &
                 & 1.0_wp - zz_frac_waste_PEM_NH - zz_frac_waste_PEM_DON - zz_frac_waste_PEM_PON
         WRITE(numout,*) '                       ZEM refractory N =', &
                 & 1.0_wp - zz_frac_waste_ZEM_NH - zz_frac_waste_ZEM_DON - zz_frac_waste_ZEM_PON
      ENDIF
      jpk40=0
      DO jk=1, jpkm1
         IF (gdept_1d(jk) .LE. 40.0_wp) THEN
            jpk40=jk
         ENDIF
      ENDDO
      zz_average_prey = 0._wp
   END SUBROUTINE p4z_mesozoo_init
   

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_mesozoo                    ! Empty routine
   END SUBROUTINE p4z_mesozoo
#endif 
   !!======================================================================
END MODULE p4zmesozoo
