MODULE p4zoxy
   !!======================================================================
   !!                         ***  MODULE p4zoxy  ***
   !! TOP :   SMELT oxygen parameter definition 
   !!=========================================================================
   !! History :   2020 (E. Olson and T. Jarnikova) adapted TJ's (and PISCES) 
   !!                   oxygen code to include additional processes and run  
   !!                   independently from SKOG
   !!----------------------------------------------------------------------
#if defined key_pisces  && defined key_oxy
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                           SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_oxy_init  :  Initialisation of parameters for oxygen calculations
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  SMELT Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_oxy_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::   zz_o2_nitr       ! ratio of O2:N uptake during nitrification (NH4->NO3)
   REAL(wp), PUBLIC ::   zz_o2_proreg     ! -O2:N oxygen evolution during conversion of ammonium to organic matter
   REAL(wp), PUBLIC ::   zz_alpha_SOD     ! fraction of N flux to sediment leading to oxygen consumption

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zoxy.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_oxy_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_oxy_init  ***
      !!
      !! ** Purpose :   Initialization of oxygen model parameters
      !!
      !! ** Method  :   Read the nampisoxy namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisoxy
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisoxy/ zz_o2_nitr, zz_o2_proreg, zz_alpha_SOD
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )     
      READ  ( numnatp_ref, nampisoxy, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisoxy in reference namelist', lwp )

      REWIND( numnatp_cfg )    
      READ  ( numnatp_cfg, nampisoxy, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisoxy in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisoxy )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for oxygen model, nampisoxy'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    O2:N ratio, nitrification                    zz_o2_nitr = ', zz_o2_nitr
         WRITE(numout,*) '    -O2:N ratio, conv of ammonium to org matter  zz_o2_proreg = ', zz_o2_proreg
         WRITE(numout,*) '    Fraction of sinking N contributing to SOD    zz_alpha_sod = ', zz_alpha_sod
      ENDIF

   END SUBROUTINE p4z_oxy_init

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_oxy                    ! Empty routine
   END SUBROUTINE p4z_oxy
#endif 

   !!======================================================================
END MODULE p4zoxy
