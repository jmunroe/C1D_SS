MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   SMELT Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   2015 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_micro       :   Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init  :   Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  SMELT Source Minus Sink variables
   USE p4zint          !  interpolation and computation of various fields
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging
   USE p4zprod, ONLY:  zz_rate_Si_ratio_diat, zz_rate_Si_ratio_myri, zz_rate_Si_ratio_flag
#if defined key_trdtrc
   USE trdtrc
   USE trd_oce
#endif
#if defined key_oxy
   USE p4zoxy
#endif
   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp) ::   zz_rate_uzoo_Rm                 ! uM N / s microzo natural mortality rate @ 20 deg C
   REAL(wp) ::   zz_rate_uzoo_excr               ! uM N / s microzo excretion rate @ 20 deg C
   REAL(wp) ::   zz_frac_waste_ZNM_NH            
   REAL(wp) ::   zz_frac_waste_ZEX_NH            
   REAL(wp) ::   zz_rate_uzoo_PredSlope          ! uM N microzo total grazing limit
   REAL(wp) ::   zz_rate_uzoo_HalfSat            ! uM N microzo total grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_diatPref          ! microzo preference for diatom
   REAL(wp) ::   zz_rate_uzoo_diatPredslope     ! uM N microzo diatom grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_diatHalfSat       ! uM N microzo diatom grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_myriPref           ! microzo preference for M rubrum
   REAL(wp) ::   zz_rate_uzoo_myriPredslope      ! uM N microzo M rubrum grazing limit
   REAL(wp) ::   zz_rate_uzoo_myriHalfSat        ! uM N microzo M rubrum grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_flagPref           ! microzo preference for small phyto
   REAL(wp) ::   zz_rate_uzoo_flagPredslope      ! uM N microzo small phyto grazing limit
   REAL(wp) ::   zz_rate_uzoo_flagHalfSat        ! uM N microzo small phyto grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_PON_Pref           ! microzo preference for PON
   REAL(wp) ::   zz_rate_uzoo_PON_Predslope      ! uM N microzo microzoo grazing limit
   REAL(wp) ::   zz_rate_uzoo_PON_HalfSat        ! uM N microzo PON grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_Z_Pref             ! microzo preference for microzoo cannibalism
   REAL(wp) ::   zz_rate_uzoo_Z_Predslope        ! uM N microzo microzoo grazing limit
   REAL(wp) ::   zz_rate_uzoo_Z_HalfSat          ! uM N microzo microzoo grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_R                  ! match SOG 0.600E-04_wp       ! 1/s diatom maximum growth rate # Hitchcock 1980, 1.4 d-1 for T. nordenskelii at 10degrees
   REAL(wp) ::   zz_frac_waste_PEZ_NH            ! waste frac from microzoo grazing PON to NH
   REAL(wp) ::   zz_frac_waste_DEZ_NH            ! waste frac from microzoo grazing diatom to NH
   REAL(wp) ::   zz_frac_waste_CEZ_NH            ! waste frac from microzoo grazing M rubrum to NH
   REAL(wp) ::   zz_frac_waste_FEZ_NH            ! waste frac from microzoo grazing small phyto to NH
   REAL(wp) ::   zz_frac_waste_ZEZ_NH            ! match SOG 1._wp               ! waste frac from micro-zoo excretion to NH
   REAL(wp) ::   zz_rate_uZoo_eff                ! match SOG
   REAL(wp) ::   zz_frac_waste_ZNM_DON           ! waste frac from micro-zoo natural mortality to DON
   REAL(wp) ::   zz_frac_waste_ZEX_DON           ! waste frac from micro-zoo excretion to DON
   REAL(wp) ::   zz_frac_waste_DEZ_DON           ! waste frac from microzoo grazing diatom to DON
   REAL(wp) ::   zz_frac_waste_FEZ_DON           ! waste frac from microzoo grazing small phyto to DON
   REAL(wp) ::   zz_frac_waste_CEZ_DON           ! waste frac from microzoo grazing M rubrum to DON
   REAL(wp) ::   zz_frac_waste_PEZ_DON           ! waste frac from microzoo grazing PON to DON
   REAL(wp) ::   zz_frac_waste_ZEZ_DON           ! waste frac from microzoo grazing microzoo to DON
   REAL(wp) ::   zz_frac_waste_ZNM_PON           ! waste frac from micro-zoo natural mortality to PON
   REAL(wp) ::   zz_frac_waste_ZEX_PON           ! waste frac from micro-zoo excretion to PON
   REAL(wp) ::   zz_frac_waste_DEZ_PON           ! waste frac from microzoo grazing diatom to PON
   REAL(wp) ::   zz_frac_waste_FEZ_PON           ! waste frac from microzoo grazing small phyto to PON
   REAL(wp) ::   zz_frac_waste_CEZ_PON           ! waste frac from microzoo grazing M rubrum to PON
   REAL(wp) ::   zz_frac_waste_PEZ_PON           ! waste frac from microzoo grazing PON to PON
   REAL(wp) ::   zz_frac_waste_ZEZ_PON           ! waste frac from microzoo grazing microzoo to PON
   REAL(wp) ::   zz_frac_waste_ZNM_BSi           ! waste frac from micro-zoo natural mortality to bSi
   REAL(wp) ::   zz_frac_waste_ZEX_Bsi           ! waste frac from micro-zoo excretion to bSi
   REAL(wp) ::   zz_frac_waste_DEZ_Bsi           ! waste frac from microzoo grazing diatom to Bsi
   REAL(wp) ::   zz_frac_waste_FEZ_BSi           ! waste frac from microzoo grazing small phyto to Bsi
   REAL(wp) ::   zz_frac_waste_CEZ_BSi           ! waste frac from microzoo grazing M rubrum to Bsi
   REAL(wp) ::   zz_frac_waste_PEZ_BSi           ! waste frac from microzoo grazing PON to Bsi
   REAL(wp) ::   zz_frac_waste_ZEZ_BSi           ! waste frac from microzoo grazing microzoo to Bsi
   
   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmicro.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_micro( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::  kt  ! ocean time step
      INTEGER, INTENT(in) ::  knt 
      INTEGER  :: ji, jj, jk, jn
      REAL(wp) ::   zz_NatMort_uzoo, zz_Excr_uzoo, zz_food_limitation, zz_denominator, &
        zz_uzoo_eat
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zz_P_diat, zz_P_myri, zz_P_flag, zz_Z, zz_D_PON, zz_uZoo_mort_diat, &
        zz_uZoo_mort_myri, zz_uZoo_mort_flag, zz_uZoo_graz_PON, zz_uZoo_graz_Z
#if defined key_trdtrc
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  ztrtra
      CHARACTER (len=20) :: strsms      
#endif
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_micro')
      CALL wrk_alloc( jpi, jpj, jpk, zz_P_diat, zz_P_myri, zz_P_flag, zz_Z, zz_D_PON, zz_uZoo_mort_diat )
      CALL wrk_alloc( jpi, jpj, jpk, zz_uZoo_mort_myri, zz_uZoo_mort_flag, zz_uZoo_graz_PON, zz_uZoo_graz_Z)

#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtra )  ! store fields before
         ztrtra(:,:,:,:)  = tra(:,:,:,:)
      ENDIF
#endif
      
      zz_P_diat(:,:,:) = trn(:,:,:,jpdia)
      zz_P_myri(:,:,:) = trn(:,:,:,jpmyr)
      zz_P_flag(:,:,:) = trn(:,:,:,jpphy)
      zz_Z(:,:,:)      = trn(:,:,:,jpzoo)
      zz_D_PON(:,:,:)  = trn(:,:,:,jppon)

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

                ! microzoo natural mortality:
                zz_NatMort_uzoo = zz_rate_uzoo_Rm * tgfunc(ji,jj,jk) * zz_Z(ji,jj,jk)
                zz_Excr_uzoo = zz_rate_uzoo_excr * tgfunc(ji,jj,jk) * zz_Z(ji,jj,jk)

                tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - (zz_NatMort_uzoo + zz_Excr_uzoo) 
                tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + (zz_frac_waste_ZNM_NH * zz_NatMort_uzoo &
                    + zz_frac_waste_ZEX_NH * zz_Excr_uzoo) 
                    
                tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + (zz_frac_waste_ZNM_DON * zz_NatMort_uzoo &
                    + zz_frac_waste_ZEX_DON * zz_Excr_uzoo) 
                tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + (zz_frac_waste_ZNM_PON * zz_NatMort_uzoo &
                    + zz_frac_waste_ZEX_PON * zz_Excr_uzoo) 
                !tra(ji,jj,jk,(Ref)) = tra(ji,jj,jk,(Ref)) + zz_frac_waste_ZNM_Ref * zz_NatMort_uzoo &
                !    + zz_frac_waste_ZEX_Ref * zz_Excr_uzoo
                ! no transfer to BSi
#if defined key_oxy
                tra(ji,jj,jk,jpo2) = tra(ji,jj,jk,jpo2) - zz_o2_proreg * (zz_frac_waste_ZNM_NH * zz_NatMort_uzoo &
                    + zz_frac_waste_ZEX_NH * zz_Excr_uzoo) 
#endif
                ! microzoo grazing:
                zz_food_limitation = (zz_P_diat(ji,jj,jk) + zz_P_myri(ji,jj,jk) + zz_P_flag(ji,jj,jk) + &
                    zz_D_PON(ji,jj,jk) + zz_Z(ji,jj,jk) - zz_rate_uzoo_PredSlope) / &
                    (zz_rate_uzoo_HalfSat + zz_P_diat(ji,jj,jk) + zz_P_myri(ji,jj,jk) + zz_P_flag(ji,jj,jk) &
                    + zz_D_PON(ji,jj,jk) + zz_Z(ji,jj,jk) - zz_rate_uzoo_Predslope &
                    + epsilon(zz_rate_uzoo_HalfSat))

                zz_denominator = (zz_rate_uzoo_diatPref * zz_P_diat(ji,jj,jk) + &
                    zz_rate_uzoo_myriPref * zz_P_myri(ji,jj,jk) + &
                    zz_rate_uzoo_flagPref * zz_P_flag(ji,jj,jk) + &
                    zz_rate_uzoo_PON_Pref * zz_D_PON(ji,jj,jk) + &
                    zz_rate_uzoo_Z_Pref * zz_Z(ji,jj,jk) + &
                    epsilon(zz_P_diat(ji,jj,jk)) )

                zz_uZoo_mort_diat(ji,jj,jk) = min(zz_rate_uzoo_diatPref * zz_food_limitation &
                    * zz_P_diat(ji,jj,jk) / zz_denominator, &
                    (zz_P_diat(ji,jj,jk) - zz_rate_uzoo_diatPredslope) / &
                    (zz_rate_uzoo_diatHalfSat + zz_P_diat(ji,jj,jk) &
                    - zz_rate_uzoo_diatPredslope + &
                    epsilon(zz_rate_uzoo_diatHalfSat)) )

                zz_uZoo_mort_myri(ji,jj,jk) = min(zz_rate_uzoo_myriPref * zz_food_limitation &
                    * zz_P_myri(ji,jj,jk) / zz_denominator, &
                    (zz_P_myri(ji,jj,jk) - zz_rate_uzoo_myriPredslope) / &
                    (zz_rate_uzoo_myriHalfSat + zz_P_myri(ji,jj,jk) &
                    - zz_rate_uzoo_myriPredslope + &
                    epsilon(zz_rate_uzoo_myriHalfSat)) )

                 zz_uZoo_mort_flag(ji,jj,jk) = min(zz_rate_uzoo_flagPref * zz_food_limitation &
                    * zz_P_flag(ji,jj,jk) / zz_denominator, &
                    (zz_P_flag(ji,jj,jk) - zz_rate_uzoo_flagPredslope) / &
                    (zz_rate_uzoo_flagHalfSat + zz_P_flag(ji,jj,jk) &
                    - zz_rate_uzoo_flagPredslope + &
                    epsilon(zz_rate_uzoo_flagHalfSat)) )

                zz_uZoo_graz_PON(ji,jj,jk) = min(zz_rate_uzoo_PON_Pref * zz_food_limitation &
                    * zz_D_PON(ji,jj,jk) / zz_denominator, &
                    (zz_D_PON(ji,jj,jk) - zz_rate_uzoo_PON_Predslope) / &
                    (zz_rate_uzoo_PON_HalfSat + zz_D_PON(ji,jj,jk) &
                    - zz_rate_uzoo_PON_Predslope + &
                    epsilon(zz_rate_uzoo_PON_HalfSat)) )

                zz_uZoo_graz_Z(ji,jj,jk) = min(zz_rate_uzoo_Z_Pref * zz_food_limitation &
                    * zz_Z(ji,jj,jk) / zz_denominator, &
                    (zz_Z(ji,jj,jk) - zz_rate_uzoo_Z_Predslope) / &
                    (zz_rate_uzoo_Z_HalfSat + zz_Z(ji,jj,jk) &
                    - zz_rate_uzoo_Z_Predslope + &
                    epsilon(zz_rate_uzoo_Z_HalfSat)) )

                zz_uZoo_mort_diat(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_mort_diat(ji,jj,jk))
                zz_uZoo_mort_myri(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_mort_myri(ji,jj,jk))
                zz_uZoo_mort_flag(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_mort_flag(ji,jj,jk))
                zz_uZoo_graz_PON(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_graz_PON(ji,jj,jk))
                zz_uZoo_graz_Z(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_graz_Z(ji,jj,jk))
                
                zz_uzoo_eat = (zz_uZoo_mort_diat(ji,jj,jk) + zz_uZoo_mort_myri(ji,jj,jk) + zz_uZoo_mort_flag(ji,jj,jk) &
                     + zz_uZoo_graz_PON(ji,jj,jk) + zz_uZoo_graz_Z(ji,jj,jk)) * zz_rate_uZoo_eff
                tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + (zz_uzoo_eat - zz_uZoo_graz_Z(ji,jj,jk)) 
                tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zz_uZoo_mort_diat(ji,jj,jk) 
                tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zz_uZoo_mort_flag(ji,jj,jk) 
                tra(ji,jj,jk,jpmyr) = tra(ji,jj,jk,jpmyr) - zz_uZoo_mort_myri(ji,jj,jk) 

                tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + ( &
                     zz_frac_waste_PEZ_NH * zz_UZoo_graz_PON(ji,jj,jk) + &
                     zz_frac_waste_DEZ_NH * zz_uZoo_mort_diat(ji,jj,jk) + &
                     zz_frac_waste_CEZ_NH * zz_uZoo_mort_myri(ji,jj,jk) + &
                     zz_frac_waste_FEZ_NH * zz_uZoo_mort_flag(ji,jj,jk) + &
                     zz_frac_waste_ZEZ_NH * zz_uZoo_graz_Z(ji,jj,jk)) * (1-zz_rate_uZoo_eff) 
                tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + ( &
                     zz_frac_waste_PEZ_DON * zz_UZoo_graz_PON(ji,jj,jk) + &
                     zz_frac_waste_DEZ_DON * zz_uZoo_mort_diat(ji,jj,jk) + &
                     zz_frac_waste_CEZ_DON * zz_uZoo_mort_myri(ji,jj,jk) + &
                     zz_frac_waste_FEZ_DON * zz_uZoo_mort_flag(ji,jj,jk) + &
                     zz_frac_waste_ZEZ_DON * zz_uZoo_graz_Z(ji,jj,jk)) * (1-zz_rate_uZoo_eff) 
                tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + ( &
                     zz_frac_waste_PEZ_PON * zz_UZoo_graz_PON(ji,jj,jk) + &
                     zz_frac_waste_DEZ_PON * zz_uZoo_mort_diat(ji,jj,jk) + &
                     zz_frac_waste_CEZ_PON * zz_uZoo_mort_myri(ji,jj,jk) + &
                     zz_frac_waste_FEZ_PON * zz_uZoo_mort_flag(ji,jj,jk) + &
                     zz_frac_waste_ZEZ_PON * zz_uZoo_graz_Z(ji,jj,jk))  * &
                     (1-zz_rate_uZoo_eff) - zz_uZoo_graz_PON(ji,jj,jk) 
                tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + ( &
                     zz_frac_waste_DEZ_BSi * zz_uZoo_mort_diat(ji,jj,jk) * zz_rate_Si_ratio_diat + &
                     zz_frac_waste_CEZ_BSi * zz_uZoo_mort_myri(ji,jj,jk) * zz_rate_Si_ratio_myri + &
                     zz_frac_waste_FEZ_BSi * zz_uZoo_mort_flag(ji,jj,jk) * zz_rate_Si_ratio_flag &
                     ) 
#if defined key_oxy
                tra(ji,jj,jk,jpo2) = tra(ji,jj,jk,jpo2) - zz_o2_proreg * ( &
                     zz_frac_waste_PEZ_NH * zz_UZoo_graz_PON(ji,jj,jk) + &
                     zz_frac_waste_DEZ_NH * zz_uZoo_mort_diat(ji,jj,jk) + &
                     zz_frac_waste_CEZ_NH * zz_uZoo_mort_myri(ji,jj,jk) + &
                     zz_frac_waste_FEZ_NH * zz_uZoo_mort_flag(ji,jj,jk) + &
                     zz_frac_waste_ZEZ_NH * zz_uZoo_graz_Z(ji,jj,jk)) * (1-zz_rate_uZoo_eff) 
#endif
            END DO
         END DO
      END DO
      
      CALL iom_put( "GRMICZDIAT", zz_uZoo_mort_diat)
      CALL iom_put( "GRMICZMRUB", zz_uZoo_mort_myri)
      CALL iom_put( "GRMICZPHY", zz_uZoo_mort_flag)
      CALL iom_put( "GRMICZPON", zz_uZoo_graz_PON)
      CALL iom_put( "GRMICZMICZ", zz_uZoo_graz_Z)

#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         strsms="MIZ_"
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn)-ztrtra(:,:,:,jn), jn, strsms, kt )   ! save trends
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtra ) 
      END IF
#endif
      CALL wrk_dealloc( jpi, jpj, jpk, zz_P_diat, zz_P_myri, zz_P_flag, zz_Z, zz_D_PON, zz_uZoo_mort_diat )
      CALL wrk_dealloc( jpi, jpj, jpk, zz_uZoo_mort_myri, zz_uZoo_mort_flag, zz_uZoo_graz_PON, zz_uZoo_graz_Z )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_micro')
      !
   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampiszoo/ zz_rate_uzoo_Rm, zz_rate_uzoo_excr, zz_frac_waste_ZNM_NH, zz_frac_waste_ZEX_NH, &
         &   zz_rate_uzoo_PredSlope, zz_rate_uzoo_HalfSat, zz_rate_uzoo_diatPref, zz_rate_uzoo_diatPredslope, &
         &   zz_rate_uzoo_diatHalfSat, zz_rate_uzoo_myriPref, zz_rate_uzoo_myriPredslope, zz_rate_uzoo_myriHalfSat, &
         &   zz_rate_uzoo_flagPref, zz_rate_uzoo_flagPredslope, zz_rate_uzoo_flagHalfSat, zz_rate_uzoo_PON_Pref, &
         &   zz_rate_uzoo_PON_Predslope, zz_rate_uzoo_PON_HalfSat, zz_rate_uzoo_Z_Pref, zz_rate_uzoo_Z_Predslope, &
         &   zz_rate_uzoo_Z_HalfSat, zz_rate_uzoo_R, zz_frac_waste_PEZ_NH, zz_frac_waste_DEZ_NH, zz_frac_waste_CEZ_NH, &
         &   zz_frac_waste_FEZ_NH, zz_frac_waste_ZEZ_NH, zz_rate_uZoo_eff, zz_frac_waste_ZNM_DON, zz_frac_waste_ZEX_DON, &
         &   zz_frac_waste_DEZ_DON, zz_frac_waste_FEZ_DON, zz_frac_waste_CEZ_DON, zz_frac_waste_PEZ_DON, &
         &   zz_frac_waste_ZEZ_DON, zz_frac_waste_ZNM_PON, zz_frac_waste_ZEX_PON, zz_frac_waste_DEZ_PON, &
         &   zz_frac_waste_FEZ_PON, zz_frac_waste_CEZ_PON, zz_frac_waste_PEZ_PON, zz_frac_waste_ZEZ_PON, &
         &   zz_frac_waste_ZNM_BSi, zz_frac_waste_ZEX_Bsi, zz_frac_waste_DEZ_Bsi, zz_frac_waste_FEZ_BSi, &
         &   zz_frac_waste_CEZ_BSi, zz_frac_waste_PEZ_BSi, zz_frac_waste_ZEZ_BSi 
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampiszoo in reference namelist : microzooplankton
      READ  ( numnatp_ref, nampiszoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiszoo in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampiszoo in configuration namelist : microzooplankton
      READ  ( numnatp_cfg, nampiszoo, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiszoo in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampiszoo )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzo, nampiszoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '  uM N / s microzo natural mort @ 20 deg C            zz_rate_uzoo_Rm=', zz_rate_uzoo_Rm
         WRITE(numout,*) '  uM N / s microzo excretion rate @ 20 deg C        zz_rate_uzoo_excr=', zz_rate_uzoo_excr
         WRITE(numout,*) '                                                 zz_frac_waste_ZNM_NH=', zz_frac_waste_ZNM_NH
         WRITE(numout,*) '                                                 zz_frac_waste_ZEX_NH=', zz_frac_waste_ZEX_NH
         WRITE(numout,*) '  uM N microzo total grazing limit             zz_rate_uzoo_PredSlope=', zz_rate_uzoo_PredSlope
         WRITE(numout,*) '  uM N microzo total grazing half sat const      zz_rate_uzoo_HalfSat=', zz_rate_uzoo_HalfSat
         WRITE(numout,*) '  microzo preference for diatom           zz_rate_uzoo_diatPref=', zz_rate_uzoo_diatPref
         WRITE(numout,*) '  uM N miczo micphyt graz half sat const  zz_rate_uzoo_diatPredslope=', zz_rate_uzoo_diatPredslope
         WRITE(numout,*) '  uM N miczo micphyt grazing half sat const zz_rate_uzoo_diatHalfSat=', zz_rate_uzoo_diatHalfSat
         WRITE(numout,*) '  microzo preference for M rubrum              zz_rate_uzoo_myriPref=', zz_rate_uzoo_myriPref
         WRITE(numout,*) '  uM N microzo M rubrum grazing limit      zz_rate_uzoo_myriPredslope=', zz_rate_uzoo_myriPredslope
         WRITE(numout,*) '  uM N microzo M rubrum graz half sat const  zz_rate_uzoo_myriHalfSat=', zz_rate_uzoo_myriHalfSat
         WRITE(numout,*) '  microzo preference for small phyto              zz_rate_uzoo_flagPref=', zz_rate_uzoo_flagPref
         WRITE(numout,*) '  uM N microzo small phyto grazing limit      zz_rate_uzoo_flagPredslope=', zz_rate_uzoo_flagPredslope
         WRITE(numout,*) '  uM N microzo small phyto graz half sat const  zz_rate_uzoo_flagHalfSat=', zz_rate_uzoo_flagHalfSat
         WRITE(numout,*) '  microzo preference for PON                    zz_rate_uzoo_PON_Pref=', zz_rate_uzoo_PON_Pref
         WRITE(numout,*) '  uM N microzo microzo grazing limit       zz_rate_uzoo_PON_Predslope=', zz_rate_uzoo_PON_Predslope
         WRITE(numout,*) '  uM N microzo PON grazing half sat const    zz_rate_uzoo_PON_HalfSat=', zz_rate_uzoo_PON_HalfSat
         WRITE(numout,*) '  microzo preference for microzo cannibalism      zz_rate_uzoo_Z_Pref=', zz_rate_uzoo_Z_Pref
         WRITE(numout,*) '  uM N microzo microzo grazing limit         zz_rate_uzoo_Z_Predslope=', zz_rate_uzoo_Z_Predslope
         WRITE(numout,*) '  uM N microzo microzo grazing half sat const  zz_rate_uzoo_Z_HalfSat=', zz_rate_uzoo_Z_HalfSat
         WRITE(numout,*) '  match SOG 0.600E-04  1/s diatom max growth rate zz_rate_uzoo_R=', zz_rate_uzoo_R
         WRITE(numout,*) '  waste frac from microzo graz PON to NH         zz_frac_waste_PEZ_NH=', zz_frac_waste_PEZ_NH
         WRITE(numout,*) '  waste frac from microzo graz diatom to NH   zz_frac_waste_DEZ_NH=', zz_frac_waste_DEZ_NH
         WRITE(numout,*) '  waste frac from microzo graz M rubrum to NH    zz_frac_waste_CEZ_NH=', zz_frac_waste_CEZ_NH
         WRITE(numout,*) '  waste frac from microzo graz small phyto to NH    zz_frac_waste_FEZ_NH=', zz_frac_waste_FEZ_NH
         WRITE(numout,*) '  waste frac from microzo excretion to NH        zz_frac_waste_ZEZ_NH=', zz_frac_waste_ZEZ_NH
         WRITE(numout,*) '                                                     zz_rate_uZoo_eff=', zz_rate_uZoo_eff
         WRITE(numout,*) '  waste frac from microzo natural mort to DON   zz_frac_waste_ZNM_DON=', zz_frac_waste_ZNM_DON
         WRITE(numout,*) '  waste frac from microzo excretion to DON      zz_frac_waste_ZEX_DON=', zz_frac_waste_ZEX_DON
         WRITE(numout,*) '  waste frac from microzo graz diatom to DON zz_frac_waste_DEZ_DON=', zz_frac_waste_DEZ_DON
         WRITE(numout,*) '  waste frac from microzo graz small phyto to DON  zz_frac_waste_FEZ_DON=', zz_frac_waste_FEZ_DON
         WRITE(numout,*) '  waste frac from microzo graz M rubrum to DON  zz_frac_waste_CEZ_DON=', zz_frac_waste_CEZ_DON
         WRITE(numout,*) '  waste frac from microzo graz PON to DON       zz_frac_waste_PEZ_DON=', zz_frac_waste_PEZ_DON
         WRITE(numout,*) '  waste frac from microzo graz microzoo to DON  zz_frac_waste_ZEZ_DON=', zz_frac_waste_ZEZ_DON
         WRITE(numout,*) '  waste frac from microzo natural mort to PON   zz_frac_waste_ZNM_PON=', zz_frac_waste_ZNM_PON
         WRITE(numout,*) '  waste frac from microzo excretion to PON      zz_frac_waste_ZEX_PON=', zz_frac_waste_ZEX_PON
         WRITE(numout,*) '  waste frac from microzo graz diatom to PON zz_frac_waste_DEZ_PON=', zz_frac_waste_DEZ_PON
         WRITE(numout,*) '  waste frac from microzo graz small phyto to PON  zz_frac_waste_FEZ_PON=', zz_frac_waste_FEZ_PON
         WRITE(numout,*) '  waste frac from microzo graz M rubrum to PON  zz_frac_waste_CEZ_PON=', zz_frac_waste_CEZ_PON
         WRITE(numout,*) '  waste frac from microzo graz PON to PON       zz_frac_waste_PEZ_PON=', zz_frac_waste_PEZ_PON
         WRITE(numout,*) '  waste frac from microzo graz microzoo to PON  zz_frac_waste_ZEZ_PON=', zz_frac_waste_ZEZ_PON
         WRITE(numout,*) '  waste frac from microzo natural mort to bSi   zz_frac_waste_ZNM_BSi=', zz_frac_waste_ZNM_BSi
         WRITE(numout,*) '  waste frac from microzo excretion to bSi      zz_frac_waste_ZEX_Bsi=', zz_frac_waste_ZEX_Bsi
         WRITE(numout,*) '  waste frac from microzo graz diatom to Bsi zz_frac_waste_DEZ_Bsi=', zz_frac_waste_DEZ_Bsi
         WRITE(numout,*) '  waste frac from microzo graz small phyto to Bsi  zz_frac_waste_FEZ_BSi=', zz_frac_waste_FEZ_BSi
         WRITE(numout,*) '  waste frac from microzo graz M rubrum to Bsi  zz_frac_waste_CEZ_BSi=', zz_frac_waste_CEZ_BSi
         WRITE(numout,*) '  waste frac from microzo graz PON to Bsi       zz_frac_waste_PEZ_BSi=', zz_frac_waste_PEZ_BSi
         WRITE(numout,*) '  waste frac from microzo graz microzoo to Bsi  zz_frac_waste_ZEZ_BSi=', zz_frac_waste_ZEZ_BSi
         WRITE(numout,*) '  diatom silicon/nitrogen ratio                 zz_rate_Si_ratio_diat=', zz_rate_Si_ratio_diat
         WRITE(numout,*) '  M rubrum silicon/nitrogen ratio             zz_rate_Si_ratio_myri=', zz_rate_Si_ratio_myri
         WRITE(numout,*) '  small phyto silicon/nitrogen ratio             zz_rate_Si_ratio_flag=', zz_rate_Si_ratio_flag
         WRITE(numout,*) '                       FEZ refractory N =', &
                 & 1.0_wp - zz_frac_waste_FEZ_NH - zz_frac_waste_FEZ_DON - zz_frac_waste_FEZ_PON
         WRITE(numout,*) '                       CEZ refractory N =', &
                 & 1.0_wp - zz_frac_waste_CEZ_NH - zz_frac_waste_CEZ_DON - zz_frac_waste_CEZ_PON
         WRITE(numout,*) '                       DEZ refractory N =', &
                 & 1.0_wp - zz_frac_waste_DEZ_NH - zz_frac_waste_DEZ_DON - zz_frac_waste_DEZ_PON
         WRITE(numout,*) '                       PEZ refractory N =', &
                 & 1.0_wp - zz_frac_waste_PEZ_NH - zz_frac_waste_PEZ_DON - zz_frac_waste_PEZ_PON
         WRITE(numout,*) '                       ZEZ refractory N =', &
                 & 1.0_wp - zz_frac_waste_ZEZ_NH - zz_frac_waste_ZEZ_DON - zz_frac_waste_ZEZ_PON
      ENDIF

   END SUBROUTINE p4z_micro_init

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_micro                    ! Empty routine
   END SUBROUTINE p4z_micro
#endif 
   !!======================================================================
END MODULE p4zmicro
