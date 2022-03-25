MODULE p4zflx
   !!======================================================================
   !!                         ***  MODULE p4zflx  ***
   !! TOP :   PISCES CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr) Include total atm P correction
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
   !! dummy module to remove this code
   !!======================================================================
#if defined key_skog || defined key_oxy
   USE oce_trc                      !  shared variables between ocean and passive tracers
   USE trc                          !  passive tracers common variables
   USE sms_pisces                   !  PISCES Source Minus Sink variables
   USE prtctl_trc                   !  print control for debugging
   USE iom                          !  I/O manager
   USE fldread                      !  read input fields
#if defined key_skog
   USE mocsywrapper                !  passes variables
   USE gastransfer
   USE mocsy_gasflux
#endif
   USE sbc_oce, ONLY :  sst_m    !  in situ temperature (surface)
#if defined key_trdtrc
   USE trdtrc
   USE trd_oce
#endif
   USE sbc_oce, ONLY :  wndm        !  wind
   USE sbcapr, ONLY  :  apr         ! atmospheric pressure
   USE par_oce, ONLY :  jp_tem, jp_sal ! indices temp sal
   USE oce, ONLY: tsb, tsn, tsa, rhop
   IMPLICIT NONE
   PRIVATE

#if defined key_skog || defined key_oxy
   PUBLIC   p4z_flx
#endif
#if defined key_skog
   PUBLIC   p4z_flx_init
   !                                 !!** Namelist  nampisext  **
   REAL(wp)          ::  zz_carbon_year
   REAL(wp)          ::  zz_LR_slope !: slope of co2 trend
   REAL(wp)          ::  zz_LR_int   !: int of co2 trend
   REAL(wp)          ::  zz_ctr      !: gaussian parameters for yearly cycle
   REAL(wp)          ::  zz_amp      !: see namelist_smelt_skog for details
   REAL(wp)          ::  zz_wid
   REAL(wp)          ::  zz_ctr2
   REAL(wp)          ::  zz_amp2
   REAL(wp)          ::  zz_wid2
   REAL(wp)          ::  zz_ctr3
   REAL(wp)          ::  zz_amp3
   REAL(wp)          ::  zz_wid3

   LOGICAL           ::  ln_co2int  !: flag to read in a file and interpolate atmospheric pco2 or not
   CHARACTER(len=34) ::  clname     !: filename of pco2 values
#endif
#if defined key_oxy
   !oxygen calculations, from p4zche from PISCES
   REAL(wp) ::  xconv  = 0.01_wp / 3600._wp
   REAL(wp) ::   atcox  = 0.20946         ! units atm
   REAL(wp) ::   o2atm  = 1. / ( 1000. * 0.20946 )
   !                                    ! volumetric solubility constants for o2 in ml/L
   REAL(wp) ::   oxyco  = 1. / 22.4144   ! converts from liters of an ideal gas
   REAL(wp) ::   ox0    =  2.00856      ! from Table 1 for Eq 8 of Garcia and Gordon, 1992.
   REAL(wp) ::   ox1    =  3.22400      ! corrects for moisture and fugacity, but not total atmospheric pressure
   REAL(wp) ::   ox2    =  3.99063      !      Original PISCES code noted this was a solubility, but
   REAL(wp) ::   ox3    =  4.80299      ! was in fact a bunsen coefficient with units L-O2/(Lsw atm-O2)
   REAL(wp) ::   ox4    =  9.78188e-1   ! Hence, need to divide EXP( zoxy ) by 1000, ml-O2 => L-O2
   REAL(wp) ::   ox5    =  1.71069      ! and atcox = 0.20946 to add the 1/atm dimension.
   REAL(wp) ::   ox6    = -6.24097e-3
   REAL(wp) ::   ox7    = -6.93498e-3
   REAL(wp) ::   ox8    = -6.90358e-3
   REAL(wp) ::   ox9    = -4.29155e-3
   REAL(wp) ::   ox10   = -3.11680e-7
#endif
   !!* Substitution
#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE p4z_flx ( kt)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx  ***
      !!
      !! ** Purpose :   CALL GAS EXCHANGE AT SEA SURFACE
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt

      INTEGER :: jn
#if defined key_trdtrc
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  ztrtra
      CHARACTER (len=20) :: strsms
#endif
      IF( nn_timing == 1 )  CALL timing_start('p4z_flx')
#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtra )  ! store fields before
         ztrtra(:,:,:,:)  = tra(:,:,:,:)
      ENDIF
#endif
#if defined key_skog
      CALL p4z_flx_c ( kt)
#endif
#if defined key_oxy
      CALL p4z_flx_o2( kt)
#endif
#if defined key_trdtrc
      IF( l_trdtrc ) THEN
         strsms="FLX_"
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_trc( tra(:,:,:,jn)-ztrtra(:,:,:,jn), jn, strsms, kt )   ! save trends
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtra ) 
      END IF
#endif
      IF( nn_timing == 1 )  CALL timing_stop('p4z_flx')
   END SUBROUTINE p4z_flx

#if defined key_skog
   SUBROUTINE p4z_flx_c ( kt)
!!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx  ***
      !!
      !! ** Purpose :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
      !!
      !! ** Method  :
      !!              - Include total atm P correction via Esbensen & Kushnir (1981)
      !!              - Pressure correction NOT done for key_cpl_carbon_cycle
      !!              - Remove Wanninkhof chemical enhancement;
      !!              - Add option for time-interpolation of atcco2.txt
      !!---------------------------------------------------------------------
      !

      INTEGER, INTENT( in ) ::   kt
      INTEGER  ::   ji, jj, jm, iind, iindm1
      REAL(wp) ::   ztc, ztc2, ztc3, ztc4, zws, zkgwan
      REAL(wp) ::   zfld16, zflu16, zfact
      REAL(wp) ::   zsch_o2
      REAL(wp), DIMENSION(1) :: dv_atmco2, dv_patm
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:) :: chemo2, zkgo2
      REAL(wp) ::   zoflx
      REAL(wp) ::   zz_year, zz_day
      REAL(wp) ::   zz_yearcyc, zz_atcco2, ztkel
      REAL(wp) ::   zsal
      REAL(wp) ::   zsal2
      REAL(wp) ::   ztgg, ztgg2, ztgg3 , ztgg4 , ztgg5
      REAL(wp) ::   zoxy
      REAL(wp) ::   zt, zt2, zlogt, zcek1
      REAL(wp) ::   sal_psu

      !
      CALL wrk_alloc( jpi, jpj, zkgo2, chemo2 )

      IF( zz_carbon_year .lt. 0) THEN
         zz_year = nyear
      ELSE
         zz_year = zz_carbon_year
      ENDIF
      zz_day = nday_year

      zz_yearcyc =  zz_amp * exp( -((zz_day - zz_ctr)/zz_wid)**2)&
              + zz_amp2 * exp( -((zz_day - zz_ctr2)/zz_wid2)**2)&
              + zz_amp3 * exp( -((zz_day - zz_ctr3)/zz_wid3)**2)

      zz_atcco2 = (zz_year+(zz_day/365))*zz_LR_slope+zz_LR_int + zz_yearcyc

      DO jj = 1, jpj
        DO ji = 1, jpi
          CALL gas_transfer( wndm(ji,jj), 1, 7,  f_kw660(ji,jj) )
        END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi

         dv_atmco2(1) = zz_atcco2
         dv_patm(1) =  apr(ji,jj)/101325

         !use mocsy_interface with in_situ T (sst_m)
         !use mocsy_interface with salinity in psu
         sal_psu = tsn(ji,jj,1,jp_sal)*35/35.16504

         IF (sst_m(ji, jj) > -3) THEN

            CALL mocsy_interface( sst_m(ji,jj), sal_psu, &
                 trn(ji,jj,1,jpta),trn(ji,jj,1,jpdic),trn(ji,jj,1,jpsil),&
                 trn(ji,jj,1,jpno3)/16.,dv_patm(1),0.,50.,f_kw660(ji,jj),&
                 dv_atmco2(1),1,&
                 f_ph(ji,jj),f_pco2w(ji,jj),f_fco2w(ji,jj),f_co2w(ji,jj),&
                 f_hco3(ji,jj),f_co3(ji,jj),f_omarg(ji,jj),f_omcal(ji,jj),&
                 f_BetaD(ji,jj),f_rhosw(ji,jj),f_opres(ji,jj),       &
                 f_insitut(ji,jj),f_pco2atm(ji,jj),f_fco2atm(ji,jj),f_schmidtco2(ji,jj), &
                 f_kwco2(ji,jj), f_K0(ji,jj), f_co2starair(ji,jj), f_co2flux(ji,jj), &
                 f_dpco2(ji,jj))
         ELSE
            f_co2flux(ji, jj) = 0.
         ENDIF

         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            !contribution of air-sea CO2 flux to dic conc
            tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic) + f_co2flux(ji,jj)/fse3t(ji,jj,1)
         END DO
      END DO

   CALL iom_put("f_co2flux",f_co2flux)

   CALL wrk_dealloc( jpi, jpj, zkgo2, chemo2)


   END SUBROUTINE p4z_flx_c
#endif

#if defined key_oxy
   SUBROUTINE p4z_flx_o2 ( kt)
!!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx  ***
      !!
      !! ** Purpose :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
      !!
      !! ** Method  :
      !!              - Include total atm P correction via Esbensen & Kushnir (1981)
      !!              - Pressure correction NOT done for key_cpl_carbon_cycle
      !!              - Remove Wanninkhof chemical enhancement;
      !!---------------------------------------------------------------------
      !

      INTEGER, INTENT( in ) ::   kt
      INTEGER  ::   ji, jj, jm, iind, iindm1
      REAL(wp) ::   ztc, ztc2, ztc3, ztc4, zws, zkgwan
      REAL(wp) ::   zfld16, zflu16, zfact
      REAL(wp) ::   zsch_o2
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:) :: chemo2, zkgo2
      REAL(wp) ::   zoflx
      REAL(wp) ::   zz_year, zz_day
      REAL(wp) ::   ztkel
      REAL(wp) ::   zsal
      REAL(wp) ::   zsal2
      REAL(wp) ::   ztgg, ztgg2, ztgg3 , ztgg4 , ztgg5
      REAL(wp) ::   zoxy
      REAL(wp) ::   zt, zt2, zlogt, zcek1
      REAL(wp) ::   sal_psu

      IF( nn_timing == 1 )  CALL timing_start('p4z_flx_o2')
      !
      CALL wrk_alloc( jpi, jpj, zkgo2, chemo2 )

      !DO jj = 1, jpj
      !  DO ji = 1, jpi
      !    CALL gas_transfer( wndm(ji,jj), 1, 7,  f_kw660(ji,jj) )
      !  END DO
      !END DO
!~~ separate section just for O2, code modified from PISCES/P4Z~~~~~~


      ! FIRST COMPUTE GAS EXCHANGE COEFFICIENTS
      ! -------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi
            ztc  = MIN( 35., tsn(ji,jj,1,jp_tem) )
            ztc2 = ztc * ztc
            ztc3 = ztc * ztc2
            ztc4 = ztc2 * ztc2
!~~~~~~~~~~~~~ F_B compute schmidt number~~~~~~~~~~~~~~~~~~
           ! Compute the schmidt Number for  O2
            zsch_o2  = 1920.4 - 135.6  * ztc + 5.2122 * ztc2 - 0.109390 * ztc3 + 0.0009377 * ztc4
!~~~~~~~~~~~~~ F_C compute piston velocity~~~~~~~~~~~~~~~~~
!            !  wind speed
           zws = wndm(ji,jj) * wndm(ji, jj)

           ! Compute the piston velocity for O2
           zkgwan = 0.251 * zws
           zkgwan = zkgwan * xconv * tmask(ji,jj,1)
# if defined key_degrad
            zkgwan = zkgwan * facvol(ji,jj,1)
#endif
           ! compute gas exchange for O2
           zkgo2(ji,jj) = zkgwan * SQRT( 660./ zsch_o2 )

         END DO
      END DO

! compute chemo2, solubility of o2 in mol/(L atm)

         DO jj= 1, jpj
            DO ji = 1, jpi
              ztkel = tsn(ji,jj,1,jp_tem) + 273.15
              zsal  = tsn(ji,jj,1,jp_sal) + ( 1.- tmask(ji,jj,1) ) * 35.
              zsal2 = zsal * zsal
              ztgg  = LOG( ( 298.15 - tsn(ji,jj,1,jp_tem) ) / ztkel )  ! Set the GORDON & GARCIA scaled temperature
              ztgg2 = ztgg  * ztgg
              ztgg3 = ztgg2 * ztgg
              ztgg4 = ztgg3 * ztgg
              ztgg5 = ztgg4 * ztgg
              zoxy  = ox0 + ox1 * ztgg + ox2 * ztgg2 + ox3 * ztgg3 + ox4 * ztgg4 + ox5 * ztgg5   &
                     + zsal * ( ox6 + ox7 * ztgg + ox8 * ztgg2 + ox9 * ztgg3 ) +  ox10 * zsal2
              chemo2(ji,jj) = ( EXP( zoxy ) * o2atm ) * oxyco * atcox     ! mol/(L atm)
            END DO
          END DO
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Compute O2 flux
            ! convert N/m2 to atm
            zfld16 = apr(ji,jj) /101325. * chemo2(ji,jj) * tmask(ji,jj,1) * zkgo2(ji,jj)*1000000    ! (mol/L) * (m/s)
            zflu16 = trn(ji,jj,1,jpo2) * tmask(ji,jj,1) * zkgo2(ji,jj)
            zoflx = zfld16 - zflu16
            !contribution of air-sea O2 flux to oxygen conc
            tra(ji,jj,1,jpo2) = tra(ji,jj,1,jpo2) + zoflx / fse3t(ji,jj,1)
         END DO
      END DO

   CALL wrk_dealloc( jpi, jpj, zkgo2, chemo2)

   IF( nn_timing == 1 )  CALL timing_stop('p4z_flx_o2')

   END SUBROUTINE p4z_flx_o2
#endif
#if defined key_skog
   SUBROUTINE p4z_flx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_flx_init  ***
      !!
      !! ** Purpose :   Initialization of atmospheric conditions
      !!
      !! ** Method  :   Read the nampisext namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !! ** input   :   Namelist nampisext
      NAMELIST/nampisext/ ln_co2int, clname,&
           zz_carbon_year, zz_LR_slope, zz_LR_int, zz_ctr, zz_amp,&
           zz_wid, zz_ctr2, zz_amp2,&
           zz_wid2, zz_ctr3, zz_amp3, zz_wid3

      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampisext in reference namelist : Pisces atm. conditions
      READ  ( numnatp_ref, nampisext, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisext in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisext in configuration namelist : Pisces atm. conditions
      READ  ( numnatp_cfg, nampisext, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisext in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisext )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for air-sea exchange, nampisext'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) 'ln_co2int = ', ln_co2int
         WRITE(numout,*) 'clname = ', clname
         WRITE(numout,*) 'zz_carbon_year = ', zz_carbon_year
         WRITE(numout,*) 'zz_LR_slope = ', zz_LR_slope
         WRITE(numout,*) 'zz_LR_int = ', zz_LR_int
         WRITE(numout,*) 'zz_ctr, zz_ctr2, zz_ctr3 = ', zz_ctr, zz_ctr2, zz_ctr3
         WRITE(numout,*) 'zz_amp, zz_amp2, zz_amp3 = ', zz_amp, zz_amp2, zz_amp3
         WRITE(numout,*) 'zz_wid, zz_wid2, zz_wid3 = ', zz_wid, zz_wid2, zz_wid3
      ENDIF
   END SUBROUTINE p4z_flx_init
#endif

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================

CONTAINS
  SUBROUTINE p4z_flx
  END SUBROUTINE p4z_flx

  SUBROUTINE p4z_flx_init
  END SUBROUTINE p4z_flx_init

  SUBROUTINE p4z_flx_o2
  END SUBROUTINE p4z_flx_o2

#endif


END MODULE p4zflx
