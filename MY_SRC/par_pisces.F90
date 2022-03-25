MODULE par_pisces
   !!======================================================================
   !!                        ***  par_pisces  ***
   !! TOP :   set the SMELT parameters
   !!======================================================================
   !! History :   2014 (E. Olson) adapted for SMELT
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_pisces.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

#if defined key_pisces
   !!---------------------------------------------------------------------
   !!   'key_pisces'   :                         SMELT / PISCES bio-model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .TRUE.  !: SMELT / PISCES flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_p4z        = .TRUE.  !: p4z flag true
   LOGICAL, PUBLIC, PARAMETER ::   lk_kriest     = .FALSE. !: Kriest flag false
#if defined key_skog && defined key_oxy
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 14      !: number of SMELT passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 14      !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 12      !: additional 3d output
#elif defined key_skog
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 13      !: number of SMELT passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 13      !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 11      !: additional 3d output
#elif defined key_oxy
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 12      !: number of SMELT passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 14      !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 12      !: additional 3d output
#else 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 11      !: number of SMELT passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 13      !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 11      !: additional 3d output
#endif
#if defined key_oxy
   LOGICAL, PUBLIC, PARAMETER ::   lk_oxy        = .TRUE.  !: oxygen flag 
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_oxy        = .FALSE. !: oxygen flag 
#endif
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_trd =  1      !: number of sms trends for SMELT

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   !    WARNING: be carefull about the order when reading the restart
        !   !!gm  this warning should be obsolet with IOM
   INTEGER, PUBLIC, PARAMETER ::   jpno3 =  1    !: Nitrates Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnh4 =  2    !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsil =  3    !: silicate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdia =  4    !: Diatoms Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  5    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpmyr =  6    !: M. rubra Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpzoo =  7    !: Microzooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdon =  8    !: dissolved organic nitrogen concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppon =  9    !: particulate organic n concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdsi = 10    !: Biogenic Silica Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpriv = 11    !: River tracer
#if defined key_skog
   INTEGER, PUBLIC, PARAMETER ::   jpdic = 12    !: dissolved inorganic carbon
   INTEGER, PUBLIC, PARAMETER ::   jpta = 13    !: total alkalinity
#endif
#if defined key_oxy && defined key_skog
   INTEGER, PUBLIC, PARAMETER ::   jpo2 = 14    !: dissolved oxygen
#elif defined key_oxy 
   INTEGER, PUBLIC, PARAMETER ::   jpo2 = 12    !: dissolved oxygen
#endif 
#else
   !!---------------------------------------------------------------------
   !!   Default                                   No CFC geochemical model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .FALSE.  !: SMELT / PISCES flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_p4z        = .FALSE.  !: p4z flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     =  0       !: No CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  =  0       !: No CFC additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  =  0       !: No CFC additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_trd =  0       !: number of sms trends for SMELT / PISCES
#endif

   ! Starting/ending bio model do-loop indices (N.B. no bio : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0     = 1                  !: First index of SMELT / PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1     = jp_pisces          !: Last  index of SMELT / PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_2d  = 1               !: First index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_2d  = jp_pisces_2d    !: Last  index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_3d  = 1               !: First index of 3D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_3d  = jp_pisces_3d    !: Last  index of 3d diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_trd = 1              !: First index of bio diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_trd = jp_pisces_trd  !: Last  index of bio diag


   !!======================================================================
END MODULE par_pisces
