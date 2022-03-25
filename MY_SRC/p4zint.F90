MODULE p4zint
   !!======================================================================
   !!                         ***  MODULE p4zint  ***
   !! TOP :   interpolation and computation of various accessory fields
   !!            - temperature dependence of biological rates
   !!======================================================================
   !! History :   2015 (E. Olson) adapted to include SOG T-dependence (Q10)
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_int        :  interpolation and computation of various accessory fields
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_int  
   
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zint.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_int( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_int  ***
      !!
      !! ** Purpose :   interpolation and computation of various accessory fields
      !!                     (temperature dependence of biological rates)
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      INTEGER  :: ji, jj                 ! dummy loop indices
      REAL(wp) :: zvar                   ! local variable
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_int')
      !
      ! Computation of phyto and zoo metabolic rate
      ! -------------------------------------------
      tgfunc (:,:,:) = EXP( 0.07 * (tsn(:,:,:,jp_tem)-20.0_wp) ) ! altered to match SOG temp_Q10
      IF (iom_use("TQ10")) THEN
        CALL iom_put( "TQ10",tgfunc) ! temperature Q10 dependence
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_int')
      !
   END SUBROUTINE p4z_int

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_int                   ! Empty routine
      WRITE(*,*) 'p4z_int: You should not have seen this print! error?'
   END SUBROUTINE p4z_int
#endif 
   !!======================================================================
END MODULE p4zint
