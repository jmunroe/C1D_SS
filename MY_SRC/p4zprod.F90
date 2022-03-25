MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP : SMELT  Growth Rates - autotrophic growth 
   !!======================================================================
   !! History :  2015 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                              SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  :   Initialization of the parameters for growth
   !!   p4z_prod_alloc :   Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
#if defined key_trdtrc
   USE trdtrc
   USE trd_oce
#endif
#if defined key_skog
   USE p4zcar 
#endif
#if defined key_oxy  
   USE p4zoxy
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_prod_alloc

   !! * Shared module variables
         ! SOG params:
      REAL(wp), PUBLIC  ::   zz_rate_R_diat              !: 1/s
      REAL(wp), PUBLIC  ::   zz_rate_R_myri              !: 1/s
      REAL(wp), PUBLIC  ::   zz_rate_R_flag              !: 1/s
      REAL(wp), PUBLIC  ::   zz_rate_maxtemp_diat        !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_maxtemp_myri        !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_maxtemp_flag        !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_temprange_diat      !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_temprange_myri      !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_temprange_flag      !: deg C
      REAL(wp), PUBLIC  ::   zz_PE_Iopt_diat           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_Iopt_myri           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_Iopt_flag           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_a_diat           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_a_myri           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_a_flag           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_b_diat           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_b_myri           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_b_flag           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_f_diat           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_f_myri           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_PE_f_flag           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_rate_K_Si_diat           !:
      REAL(wp), PUBLIC  ::   zz_rate_K_Si_myri           !:
      REAL(wp), PUBLIC  ::   zz_rate_K_Si_flag           !:
      REAL(wp), PUBLIC  ::   zz_rate_kapa_diat           !:
      REAL(wp), PUBLIC  ::   zz_rate_kapa_myri           !:
      REAL(wp), PUBLIC  ::   zz_rate_kapa_flag           !:
      REAL(wp), PUBLIC  ::   zz_rate_k_diat              !:
      REAL(wp), PUBLIC  ::   zz_rate_k_myri              !:
      REAL(wp), PUBLIC  ::   zz_rate_k_flag              !:
      REAL(wp), PUBLIC  ::   zz_rate_Si_ratio_diat       !:
      REAL(wp), PUBLIC  ::   zz_rate_Si_ratio_myri       !:
      REAL(wp), PUBLIC  ::   zz_rate_Si_ratio_flag       !:
      REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zz_diat_Nlimit   !: diatom N limitation for p4zsink 

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zprod.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Compute PP depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER  ::   ji, jj, jk, jn
      
      REAL(wp) ::   zproreg, zproreg2, zproregm
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zprorca, zprorcad, zprorcam, zpronew, zpronewd, zpronewm

      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zz_Si, zz_NH, zz_NO, zz_temp, &
         zz_plank_growth_diat, zz_plank_growth_myri, zz_plank_growth_flag, &
         zz_P_diat, zz_P_myri, zz_P_flag, &
         zz_uptake_NO_diat, zz_uptake_NO_myri, zz_uptake_NO_flag, &
         zz_uptake_NH_diat, zz_uptake_NH_myri, zz_uptake_NH_flag, &
         zz_uptake_PC_diat, zz_uptake_PC_myri, zz_uptake_PC_flag, &
         zz_dummy
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  ztrtra
      CHARACTER (len=20) :: strsms

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_prod')
      !
      !  Allocate temporary workspace
#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtra )  ! store fields before
         ztrtra(:,:,:,:)  = tra(:,:,:,:)
      ENDIF
#endif
      CALL wrk_alloc( jpi, jpj, jpk, zprorca, zprorcad, zprorcam, zpronew, zpronewd, zpronewm )
      CALL wrk_alloc( jpi, jpj, jpk, zz_Si, zz_NH, zz_NO, zz_temp, zz_plank_growth_diat, zz_plank_growth_myri, zz_plank_growth_flag)
      CALL wrk_alloc( jpi, jpj, jpk, zz_P_diat, zz_P_myri, zz_P_flag, zz_uptake_NO_diat, zz_uptake_NO_myri )
      CALL wrk_alloc( jpi, jpj, jpk, zz_uptake_NO_flag, zz_dummy )
      CALL wrk_alloc( jpi, jpj, jpk, zz_uptake_NH_diat, zz_uptake_PC_diat, zz_uptake_NH_myri)
      CALL wrk_alloc( jpi, jpj, jpk, zz_uptake_PC_myri, zz_uptake_NH_flag, zz_uptake_PC_flag)     
      !
      zprorca (:,:,:) = 0._wp
      zprorcad(:,:,:) = 0._wp
      zprorcam(:,:,:) = 0._wp
      zpronew (:,:,:) = 0._wp
      zpronewd(:,:,:) = 0._wp
      zpronewm(:,:,:) = 0._wp
      zz_plank_growth_diat(:,:,:) = 0._wp
      zz_plank_growth_myri(:,:,:) = 0._wp
      zz_plank_growth_flag(:,:,:) = 0._wp
      zz_uptake_NO_diat(:,:,:) = 0._wp
      zz_uptake_NO_myri(:,:,:) = 0._wp
      zz_uptake_NO_flag(:,:,:) = 0._wp

      ! set SMELT vars:
      zz_P_diat(:,:,:)  = trn(:,:,:,jpdia)
      zz_P_myri(:,:,:)  = trn(:,:,:,jpmyr)
      zz_P_flag(:,:,:) = trn(:,:,:,jpphy)
      zz_Si(:,:,:) = trn(:,:,:,jpsil)
      zz_NH(:,:,:) = trn(:,:,:,jpnh4)
      zz_NO(:,:,:) = trn(:,:,:,jpno3)
      zz_temp(:,:,:) = tsn(:,:,:,jp_tem)

      !epar is PAR

      !diatoms:
      CALL p4z_growthSOG(zz_NO, zz_NH, zz_Si, zz_P_diat, epar, zz_temp, zz_rate_R_diat, zz_rate_maxtemp_diat, &
        zz_rate_temprange_diat, zz_PE_a_diat, zz_PE_b_diat, zz_PE_f_diat, zz_rate_K_Si_diat, &
        zz_rate_kapa_diat, zz_rate_k_diat, zz_rate_Si_ratio_diat, zz_plank_growth_diat, zz_uptake_NO_diat, zz_uptake_NH_diat, zz_uptake_PC_diat, zz_diat_Nlimit)
      !M. rubrum:
      CALL p4z_growthSOG(zz_NO, zz_NH, zz_Si, zz_P_myri, epar, zz_temp, zz_rate_R_myri, zz_rate_maxtemp_myri, &
        zz_rate_temprange_myri, zz_PE_a_myri, zz_PE_b_myri, zz_PE_f_myri, zz_rate_K_Si_myri, &
        zz_rate_kapa_myri, zz_rate_k_myri, zz_rate_Si_ratio_myri, zz_plank_growth_myri, zz_uptake_NO_myri, zz_uptake_NH_myri, zz_uptake_PC_myri, zz_dummy)
      !small phyto:
      CALL p4z_growthSOG(zz_NO, zz_NH, zz_Si, zz_P_flag, epar, zz_temp, zz_rate_R_flag, zz_rate_maxtemp_flag, &
        zz_rate_temprange_flag, zz_PE_a_flag, zz_PE_b_flag, zz_PE_f_flag, zz_rate_K_Si_flag, &
        zz_rate_kapa_flag, zz_rate_k_flag, zz_rate_Si_ratio_flag, zz_plank_growth_flag, zz_uptake_NO_flag, zz_uptake_NH_flag, zz_uptake_PC_flag, zz_dummy)

#if defined key_skog      
      CALL p4z_car( kt, knt, zz_uptake_NO_diat, zz_uptake_NO_flag, zz_uptake_NO_myri, &
        zz_uptake_NH_diat, zz_uptake_NH_flag, zz_uptake_NH_myri, & 
        zz_uptake_PC_diat, zz_uptake_PC_flag, zz_uptake_PC_myri)
  
      CALL p4z_ta( kt, knt, zz_uptake_NO_diat, zz_uptake_NO_flag, zz_uptake_NO_myri, &
        zz_uptake_NH_diat, zz_uptake_NH_flag, zz_uptake_NH_myri )    
#endif

      ! Computation of the various production terms
!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               !  production terms for small phyto zprbio=zz_plank_growth_flag
               zprorca(ji,jj,jk) = zz_plank_growth_flag(ji,jj,jk)  * zz_P_flag(ji,jj,jk) 
               zpronew(ji,jj,jk) = zz_uptake_NO_flag(ji,jj,jk) 

               !  production terms for Mesodinium rubrum
               zprorcam(ji,jj,jk) = zz_plank_growth_myri(ji,jj,jk) * zz_P_myri(ji,jj,jk) 
               zpronewm(ji,jj,jk) = zz_uptake_NO_myri(ji,jj,jk) 
                 
               !  production terms for diatoms
               zprorcad(ji,jj,jk) = zz_plank_growth_diat(ji,jj,jk) * zz_P_diat(ji,jj,jk) 
               zpronewd(ji,jj,jk) = zz_uptake_NO_diat(ji,jj,jk) 
            END DO
         END DO
      END DO

      CALL iom_put( "PAR",epar) ! PAR (as used in bio model)
      CALL iom_put( "PPDIAT",zz_plank_growth_diat * zz_P_diat) ! Diatom primary productivity (umol N/s)
      CALL iom_put( "PPPHY",zz_plank_growth_flag * zz_P_flag) ! Small phyto primary productivity (umol N/s)
      CALL iom_put( "PPMRUB",zz_plank_growth_myri * zz_P_myri) ! M. Rubrum primary productivity (umol N/s)
      CALL iom_put( "PPDIATNO3",zz_uptake_NO_diat) ! Diatom primary productivity (umol N/s)
      CALL iom_put( "PPPHYNO3",zz_uptake_NO_flag) ! Small phyto primary productivity (umol N/s)
      CALL iom_put( "PPMRUBNO3",zz_uptake_NO_myri) ! M. Rubrum primary productivity (umol N/s)
      CALL iom_put( "PPDIATNO3V",zz_uptake_NO_diat*cvol) ! Diatom primary productivity (umol N/s)
      CALL iom_put( "PPPHYNO3V",zz_uptake_NO_flag*cvol) ! Small phyto primary productivity (umol N/s)
      CALL iom_put( "PPMRUBNO3V",zz_uptake_NO_myri*cvol) ! M. Rubrum primary productivity (umol N/s)
#if defined key_trdtrc && defined key_skog
      IF ( iom_use("UNC_DIC" )) THEN 
      strsms = "UNC_"        
      CALL trd_trc( zz_uptake_PC_diat(:,:,:)+zz_uptake_PC_flag(:,:,:)+zz_uptake_PC_myri(:,:,:), jpdic, strsms, kt )
      ENDIF
#endif
      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              zproreg  = zprorca(ji,jj,jk) - zpronew(ji,jj,jk)
              zproreg2 = zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)
              zproregm = zprorcam(ji,jj,jk) - zpronewm(ji,jj,jk)
              tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronew(ji,jj,jk) - zpronewd(ji,jj,jk) - zpronewm(ji,jj,jk)
              tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproreg - zproreg2 - zproregm
              tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorca(ji,jj,jk)
              tra(ji,jj,jk,jpmyr) = tra(ji,jj,jk,jpmyr) + zprorcam(ji,jj,jk)
              tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk)
              tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - (zprorcad(ji,jj,jk)*zz_rate_Si_ratio_diat + zprorca(ji,jj,jk)*zz_rate_Si_ratio_flag + zprorcam(ji,jj,jk)*zz_rate_Si_ratio_myri)
#if defined key_oxy
              tra(ji,jj,jk,jpo2)  = tra(ji,jj,jk,jpo2) + zz_o2_proreg * (zproreg + zproreg2 + zproregm) &
                               & + (zz_o2_proreg + zz_o2_nitr) * (zpronew(ji,jj,jk) + zpronewd(ji,jj,jk) + zpronewm(ji,jj,jk))
#endif
          END DO
        END DO
      END DO

#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         strsms="PRD_"
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn)-ztrtra(:,:,:,jn), jn, strsms, kt )   ! save trends
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtra ) 
      END IF
#endif
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zprorca, zprorcad, zprorcam, zpronew, zpronewd, zpronewm )
      CALL wrk_dealloc( jpi, jpj, jpk, zz_Si, zz_NH, zz_NO, zz_temp, zz_plank_growth_diat, zz_plank_growth_myri, zz_plank_growth_flag)
      CALL wrk_dealloc( jpi, jpj, jpk,  zz_P_diat, zz_P_myri, zz_P_flag, zz_uptake_NO_diat, zz_uptake_NO_myri )
      CALL wrk_dealloc( jpi, jpj, jpk,  zz_uptake_NO_flag, zz_dummy )
      CALL wrk_dealloc(jpi, jpj, jpk, zz_uptake_NH_diat, zz_uptake_PC_diat, zz_uptake_NH_myri, zz_uptake_PC_myri, zz_uptake_NH_flag, zz_uptake_PC_flag)

     IF( nn_timing == 1 )  CALL timing_stop('p4z_prod')
     !
   END SUBROUTINE p4z_prod

   SUBROUTINE p4z_growthSOG(zz_NO, zz_NH, zz_Si, zz_P, zz_I_par, zz_temp, zz_rate_R, zz_rate_maxtemp, &
        zz_rate_temprange, zz_PE_a, zz_PE_b, zz_PE_f, zz_rate_K_Si, &
        zz_rate_kapa, zz_rate_k, zz_rate_Si_ratio, zz_plank_growth, zz_uptake_NO, &
        zz_uptake_NH, zz_uptake_PC, zz_plank_Nlimit)
     !! calculate growth rates and nutrient utilization
        !! light and nutrient limitation
     INTEGER  ::   ji, jj, jk
     REAL(wp), INTENT(in)   :: zz_NO(:,:,:), zz_NH(:,:,:), zz_Si(:,:,:), zz_P(:,:,:), &
        zz_I_par(:,:,:), zz_temp(:,:,:)        ! 3d input arrays, DIMENSION(jpi,jpj,jpk)
     REAL(wp), INTENT(in)  :: zz_rate_R, zz_rate_maxtemp, zz_rate_temprange, & 
        zz_PE_a, zz_PE_b, zz_PE_f, zz_rate_K_Si, zz_rate_kapa, zz_rate_k, &
        zz_rate_Si_ratio  ! parameters
     REAL(wp), INTENT(out), DIMENSION(jpi,jpj,jpk)   :: zz_plank_growth, zz_uptake_NO, &
        zz_uptake_NH, zz_uptake_PC, zz_plank_Nlimit ! 3d output arrays
   
     REAL(wp), POINTER, DIMENSION(:,:,:) ::   zz_Rmax(:,:,:), zz_plank_growth_light(:,:,:), &
        zz_Sc(:,:,:), zz_Oup_cell(:,:,:), zz_Hup_cell(:,:,:)
   
         ! Allocate temporary workspace vars
     CALL wrk_alloc( jpi, jpj, jpk, zz_Rmax, zz_plank_growth_light, zz_Sc, zz_Oup_cell, zz_Hup_cell)
        
   ! Q10 effect SOG; replaced temp_Q10 with tgfunc and
      ! KtoC(temp(j)) with tsn(ji,jj,jk,jp_tem)
      ! diatoms:
      zz_Rmax(:,:,:)= zz_rate_R * tgfunc(:,:,:) &
            * min(max(zz_rate_maxtemp - zz_temp, 0.0_wp), &
                  zz_rate_temprange) / &
            (zz_rate_temprange + epsilon(zz_rate_temprange))

      
   ! Computation of SOG limitation terms: light, nutrients
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1, jpi
                  ! LIGHT (now based on Platt 1981 but scaled to nearly match previous version
                  ! and preserve maximum growth rate)
                  zz_plank_growth_light(ji,jj,jk) = &
                       (1.0_wp - exp(-zz_I_par(ji,jj,jk) *zz_PE_a) ) * &
                       (exp(-zz_I_par(ji,jj,jk) *zz_PE_b )) * zz_PE_f

                  ! Si
                  zz_Sc(ji,jj,jk) = zz_Si(ji,jj,jk) / (zz_rate_K_Si + zz_Si(ji,jj,jk)+epsilon(zz_rate_K_Si))

                  ! Nitrate and Ammonium
                  IF (zz_NO(ji,jj,jk) > epsilon(zz_NO(ji,jj,jk))) THEN
                     zz_Oup_cell(ji,jj,jk) = zz_NO(ji,jj,jk) * zz_rate_kapa / &
                          (zz_rate_k + zz_NO(ji,jj,jk) * zz_rate_kapa + &
                          zz_NH(ji,jj,jk)+epsilon(zz_NH))
                  ELSE
                     zz_Oup_cell(ji,jj,jk) = 0._wp
                  ENDIF
                  IF (zz_NH(ji,jj,jk) > epsilon(zz_NH(ji,jj,jk))) THEN
                     zz_Hup_cell(ji,jj,jk) = zz_NH(ji,jj,jk) / &
                          (zz_rate_k + zz_NO(ji,jj,jk) * zz_rate_kapa + &
                          zz_NH(ji,jj,jk))
                  ELSE
                     zz_Hup_cell(ji,jj,jk) = 0._wp
                  ENDIF

                  IF (zz_Oup_cell(ji,jj,jk) < 0._wp) THEN
                     WRITE(numout,*) "Oup_cell(ji,jj,jk) < 0. in NPZD.f90"
                     WRITE(numout,*) zz_Oup_cell(ji,jj,jk)
                     CALL EXIT(1)
                  ENDIF
                  IF (zz_Hup_cell(ji,jj,jk) < 0._wp) THEN
                     WRITE(numout,*) "Hup_cell(ji,jj,jk) < 0. in NPZD.f90"
                     WRITE(numout,*) zz_Hup_cell(ji,jj,jk)
                     CALL EXIT(1)
                  ENDIF

                  ! exponent of 1/5 follows Alain
                  zz_plank_Nlimit(ji,jj,jk) = min((zz_Oup_cell(ji,jj,jk) + &
                     zz_Hup_cell(ji,jj,jk)), zz_Sc(ji,jj,jk))**0.2_wp

                  ! Choose light limitation or nutrient limitation
                  IF (zz_plank_growth_light(ji,jj,jk) < 0._wp ) THEN
                     zz_plank_growth(ji,jj,jk) = 0._wp
                  ELSE

                     IF (min(zz_plank_growth_light(ji,jj,jk),zz_Sc(ji,jj,jk)) >= &
                          zz_Oup_cell(ji,jj,jk) + zz_Hup_cell(ji,jj,jk)) THEN
                        !N LIMITING
                        zz_plank_growth(ji,jj,jk) = zz_Rmax(ji,jj,jk) * &
                             (zz_Oup_cell(ji,jj,jk) + zz_Hup_cell(ji,jj,jk))

                        IF (zz_plank_growth(ji,jj,jk) < 0._wp) THEN
                           zz_plank_growth(ji,jj,jk) = 0._wp
                        ENDIF
                        zz_uptake_NO(ji,jj,jk) = zz_Rmax(ji,jj,jk) * zz_Oup_cell(ji,jj,jk) * &
                             zz_P(ji,jj,jk)
                        zz_uptake_NH(ji,jj,jk) = zz_Rmax(ji,jj,jk) * zz_Hup_cell(ji,jj,jk) * &
                             zz_P(ji,jj,jk)

                     ELSE
                        IF (zz_plank_growth_light(ji,jj,jk) < zz_Sc(ji,jj,jk)) THEN
                        !LIGHT LIMITING
                           zz_plank_growth(ji,jj,jk) = zz_Rmax(ji,jj,jk) * zz_plank_growth_light(ji,jj,jk)
                        ELSE
                        ! Si limitation
                           zz_plank_growth(ji,jj,jk) = zz_Rmax(ji,jj,jk) * zz_Sc(ji,jj,jk)
                        ENDIF
                        ! split the nitrogen uptake between NH and NO
                        IF (zz_plank_growth(ji,jj,jk) <= &
                             zz_Rmax(ji,jj,jk) * zz_Hup_cell(ji,jj,jk)) THEN
                           zz_uptake_NH(ji,jj,jk) = zz_plank_growth(ji,jj,jk) * &
                                zz_P(ji,jj,jk) 
                           zz_uptake_NO(ji,jj,jk) = 0.0_wp 
                        ELSE
                           zz_uptake_NH(ji,jj,jk) = zz_Rmax(ji,jj,jk) * &
                                zz_Hup_cell(ji,jj,jk) * zz_P(ji,jj,jk)
                           zz_uptake_NO(ji,jj,jk) = (zz_plank_growth(ji,jj,jk) - &
                                zz_Rmax(ji,jj,jk) * zz_Hup_cell(ji,jj,jk)) &
                                * zz_P(ji,jj,jk)

                        ENDIF
                     ENDIF

                     ! PC Differential Carbon Uptake
                     ! Chlr Redux Factor = 0.2 (Ianson and Allen 2002)
                     zz_uptake_PC(ji,jj,jk) = (zz_Rmax(ji,jj,jk) * zz_plank_growth_light(ji,jj,jk) - &
                          zz_plank_growth(ji,jj,jk)) &
                          * 0.2_wp * zz_P(ji,jj,jk)
                  ENDIF
           END DO
        END DO
     END DO

     
      ! Deallocate temporary workspace SOG vars 
     CALL wrk_dealloc( jpi, jpj, jpk, zz_Rmax, zz_plank_growth_light, zz_Sc, zz_Oup_cell, zz_Hup_cell)
   END SUBROUTINE P4z_growthSOG

   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of PP parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      !
      NAMELIST/nampisprod/ zz_rate_R_diat, zz_rate_R_myri, zz_rate_R_flag, zz_rate_maxtemp_diat, &
      &        zz_rate_maxtemp_myri, zz_rate_maxtemp_flag, zz_rate_temprange_diat, zz_rate_temprange_myri, &
      &        zz_rate_temprange_flag, zz_PE_a_diat, zz_PE_a_myri, zz_PE_a_flag, zz_PE_Iopt_diat, zz_PE_Iopt_myri, &
      &        zz_PE_Iopt_flag, zz_rate_K_Si_diat, zz_rate_K_Si_myri, zz_rate_K_Si_flag, &
      &        zz_rate_kapa_diat, zz_rate_kapa_myri, zz_rate_kapa_flag, zz_rate_k_diat, zz_rate_k_myri, &
      &        zz_rate_k_flag, &
      &        zz_rate_Si_ratio_diat, zz_rate_Si_ratio_myri, zz_rate_Si_ratio_flag
      INTEGER :: ios                 ! Local integer output status for namelist read
      REAL(wp) :: zz_PE_alpha_diat, zz_PE_beta_diat, zz_PE_alpha_myri, zz_PE_beta_myri, zz_PE_alpha_flag, zz_PE_beta_flag
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampisprod in reference namelist : SMELT phytoplankton production
      READ  ( numnatp_ref, nampisprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisprod in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisprod in configuration namelist : SMELT phytoplankton production
      READ  ( numnatp_cfg, nampisprod, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisprod in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisprod )
      
      IF(zz_PE_a_diat>0 .AND. zz_PE_Iopt_diat>0 .AND. zz_rate_R_diat>0) THEN
         zz_PE_alpha_diat=zz_PE_a_diat*zz_rate_R_diat
         zz_PE_beta_diat=zz_PE_alpha_diat/(exp(zz_PE_Iopt_diat*zz_PE_alpha_diat/zz_rate_R_diat)-1._wp)
         zz_PE_b_diat=zz_PE_beta_diat/zz_rate_R_diat
         zz_PE_f_diat=(zz_PE_alpha_diat+zz_PE_beta_diat)/zz_PE_alpha_diat*((zz_PE_alpha_diat+zz_PE_beta_diat)/zz_PE_beta_diat)**(zz_PE_beta_diat/zz_PE_alpha_diat)
      ELSE
         zz_PE_b_diat=0.0_wp
         zz_PE_f_diat=0.0_wp
         zz_PE_beta_diat=0.0_wp
      ENDIF
      IF(zz_PE_a_myri>0 .AND. zz_PE_Iopt_myri>0 .AND. zz_rate_R_myri>0) THEN
         zz_PE_alpha_myri=zz_PE_a_myri*zz_rate_R_myri
         zz_PE_beta_myri=zz_PE_alpha_myri/(exp(zz_PE_Iopt_myri*zz_PE_alpha_myri/zz_rate_R_myri)-1._wp)
         zz_PE_b_myri=zz_PE_beta_myri/zz_rate_R_myri
         zz_PE_f_myri=(zz_PE_alpha_myri+zz_PE_beta_myri)/zz_PE_alpha_myri*((zz_PE_alpha_myri+zz_PE_beta_myri)/zz_PE_beta_myri)**(zz_PE_beta_myri/zz_PE_alpha_myri)
      ELSE
         zz_PE_b_myri=0.0_wp
         zz_PE_f_myri=0.0_wp
         zz_PE_beta_myri=0.0_wp
      ENDIF
      IF(zz_PE_a_flag>0 .AND. zz_PE_Iopt_flag>0 .AND. zz_rate_R_flag>0) THEN
         zz_PE_alpha_flag=zz_PE_a_flag*zz_rate_R_flag
         zz_PE_beta_flag=zz_PE_alpha_flag/(exp(zz_PE_Iopt_flag*zz_PE_alpha_flag/zz_rate_R_flag)-1._wp)
         zz_PE_b_flag=zz_PE_beta_flag/zz_rate_R_flag
         zz_PE_f_flag=(zz_PE_alpha_flag+zz_PE_beta_flag)/zz_PE_alpha_flag*((zz_PE_alpha_flag+zz_PE_beta_flag)/zz_PE_beta_flag)**(zz_PE_beta_flag/zz_PE_alpha_flag) 
      ELSE
         zz_PE_b_flag=0.0_wp
         zz_PE_f_flag=0.0_wp
         zz_PE_beta_flag=0.0_wp
      ENDIF

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, nampisprod'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) ' zz_rate_R_diat =', zz_rate_R_diat
         WRITE(numout,*) ' zz_rate_R_myri =',  zz_rate_R_myri
         WRITE(numout,*) ' zz_rate_R_flag =', zz_rate_R_flag
         WRITE(numout,*) ' zz_rate_maxtemp_diat =', zz_rate_maxtemp_diat
         WRITE(numout,*) ' zz_rate_maxtemp_myri =', zz_rate_maxtemp_myri
         WRITE(numout,*) ' zz_rate_maxtemp_flag =', zz_rate_maxtemp_flag
         WRITE(numout,*) ' zz_rate_temprange_diat =', zz_rate_temprange_diat
         WRITE(numout,*) ' zz_rate_temprange_myri =', zz_rate_temprange_myri
         WRITE(numout,*) ' zz_rate_temprange_flag =', zz_rate_temprange_flag
         WRITE(numout,*) ' zz_PE_Iopt_diat =', zz_PE_Iopt_diat
         WRITE(numout,*) ' zz_PE_Iopt_myri =', zz_PE_Iopt_myri
         WRITE(numout,*) ' zz_PE_Iopt_flag =', zz_PE_Iopt_flag
         WRITE(numout,*) ' zz_PE_a_diat =', zz_PE_a_diat
         WRITE(numout,*) ' zz_PE_a_myri =', zz_PE_a_myri
         WRITE(numout,*) ' zz_PE_a_flag =', zz_PE_a_flag
         WRITE(numout,*) ' zz_rate_K_Si_diat =', zz_rate_K_Si_diat
         WRITE(numout,*) ' zz_rate_K_Si_myri =', zz_rate_K_Si_myri
         WRITE(numout,*) ' zz_rate_K_Si_flag =', zz_rate_K_Si_flag
         WRITE(numout,*) ' zz_rate_kapa_diat =', zz_rate_kapa_diat
         WRITE(numout,*) ' zz_rate_kapa_myri =', zz_rate_kapa_myri
         WRITE(numout,*) ' zz_rate_kapa_flag =', zz_rate_kapa_flag
         WRITE(numout,*) ' zz_rate_k_diat =', zz_rate_k_diat
         WRITE(numout,*) ' zz_rate_k_myri =', zz_rate_k_myri
         WRITE(numout,*) ' zz_rate_k_flag =', zz_rate_k_flag
         WRITE(numout,*) ' zz_rate_Si_ratio_diat =', zz_rate_Si_ratio_diat
         WRITE(numout,*) ' zz_rate_Si_ratio_myri =', zz_rate_Si_ratio_myri
         WRITE(numout,*) ' zz_rate_Si_ratio_flag =', zz_rate_Si_ratio_flag
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) ' Calculated light dependence parameters: '
         WRITE(numout,*) ' Diatom      Initial PE slope (mmol N/s/(W/m2)) =', zz_PE_alpha_diat
         WRITE(numout,*) ' M. rubrum   Initial PE slope (mmol N/s/(W/m2)) =', zz_PE_alpha_myri
         WRITE(numout,*) ' Small Phyto Initial PE slope (mmol N/s/(W/m2)) =', zz_PE_alpha_flag
         WRITE(numout,*) ' zz_PE_b_diat =', zz_PE_b_diat
         WRITE(numout,*) ' zz_PE_b_myri =', zz_PE_b_myri
         WRITE(numout,*) ' zz_PE_b_flag =', zz_PE_b_flag
         WRITE(numout,*) ' zz_PE_f_diat =', zz_PE_f_diat
         WRITE(numout,*) ' zz_PE_f_myri =', zz_PE_f_myri
         WRITE(numout,*) ' zz_PE_f_flag =', zz_PE_f_flag
       ENDIF
      !zz_PE_a_diat = 0.07215007215007214_wp
      !zz_PE_a_myri  = 0.00_wp
      !zz_PE_a_flag  = 0.303030303030303_wp
      !zz_PE_b_diat = 0.0007936507936507937_wp
      !zz_PE_b_myri  = 0.00_wp
      !zz_PE_b_flag  = 0.0033333333333333335_wp
      !zz_PE_f_diat = 1.06_wp
      !zz_PE_f_myri  = 1.06_wp
      !zz_PE_f_flag  = 1.06_wp

   END SUBROUTINE p4z_prod_init

   INTEGER FUNCTION p4z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( zz_diat_Nlimit(jpi,jpj,jpk), STAT = p4z_prod_alloc )
      !
      IF( p4z_prod_alloc /= 0 ) CALL ctl_warn('p4z_prod_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_prod_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_prod                    ! Empty routine
   END SUBROUTINE p4z_prod
#endif 
   !!======================================================================
END MODULE p4zprod
