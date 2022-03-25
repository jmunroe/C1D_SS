MODULE sms_pisces   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   SMELT Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!             2014 (E. Olson) adapted for SMELT
   !!----------------------------------------------------------------------
#if defined key_pisces || defined key_pisces_reduced 
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                  SMELT / PISCES model
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc

   IMPLICIT NONE
   PUBLIC

   INTEGER ::   numnatp_ref = -1           !! Logical units for namelist pisces
   INTEGER ::   numnatp_cfg = -1           !! Logical units for namelist pisces
   INTEGER ::   numonp      = -1           !! Logical unit for namelist pisces output

   !!*  Biological fluxes for light : variables shared by pisces & lobster
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  neln  !: number of T-levels + 1 in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup  !: euphotic layer depth
   
#if defined key_pisces 
   !!*  Time variables
   INTEGER  ::   nrdttrc           !: ???
   INTEGER  ::   ndayflxtr         !: ???
   REAL(wp) ::   rfact , rfactr    !: ???
   REAL(wp) ::   rfact2, rfact2r   !: ???
   REAL(wp) ::   xstep             !: Time step duration for biology
   REAL(wp) ::   ryyss             !: number of seconds per year 
   REAL(wp) ::   r1_ryyss          !: inverse number of seconds per year 

   !!* restoring
   LOGICAL  ::  ln_pisdmp          !: restoring or not of nutrients to a mean value
   INTEGER  ::  nn_pisdmp          !: frequency of relaxation or not of nutrients to a mean value

   !!* Mass conservation
   LOGICAL  ::  ln_check_mass      !: Flag to check mass conservation

   !!* Temperature dependancy of SMS terms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc    !: Temp. dependancy of various biological rates

#endif
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: sms_pisces.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sms_pisces_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_pisces_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_warn
      INTEGER ::   ierr(5)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light : shared variables for pisces & lobster
      ALLOCATE( neln(jpi,jpj), heup(jpi,jpj), STAT=ierr(1) )
      !
#if defined key_pisces
      !* Temperature dependancy of SMS terms
      ALLOCATE( tgfunc(jpi,jpj,jpk)  , STAT=ierr(5) )
         !
#endif
      !
      sms_pisces_alloc = MAXVAL( ierr )
      !
      IF( sms_pisces_alloc /= 0 )   CALL ctl_warn('sms_pisces_alloc: failed to allocate arrays') 
      !
   END FUNCTION sms_pisces_alloc

#else
   !!----------------------------------------------------------------------   
   !!  Empty module :                                     NO bio model
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================   
END MODULE sms_pisces    
