MODULE trcnam_pisces
   !!======================================================================
   !!                      ***  MODULE trcnam_pisces  ***
   !! TOP :   initialisation of some run parameters for SMELT bio-model
   !!======================================================================
   !! History :    -   !  1999-10 (M.A. Foujols, M. Levy) original code
   !!              -   !  2000-01 (L. Bopp) hamocc3, p3zd
   !!             1.0  !  2003-08 (C. Ethe)  module F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.pisces.h90i
   !!              2014 (E. Olson) adapted for SMELT
   !!----------------------------------------------------------------------
#if defined key_pisces 
   !!----------------------------------------------------------------------
   !!   'key_pisces'   :                           SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_nam_pisces       : SMELT model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_pisces      ! sms trends
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_pisces   ! called by trcnam.F90 module


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_pisces.F90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_pisces
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_pisces  ***  
      !!
      !! ** Purpose :   read SMELT namelist
      !!----------------------------------------------------------------------
      !!
      INTEGER :: jl, jn
      INTEGER :: ios                 ! Local integer output status for namelist read
      CHARACTER(LEN=20)   ::   clname
      !!
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      clname = 'namelist_smelt'
      IF(lwp) WRITE(numout,*) ' trc_nam_pisces : read SMELT namelist'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      CALL ctl_opn( numnatp_ref, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfg, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonp     , 'output.namelist.sme' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

   END SUBROUTINE trc_nam_pisces

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_pisces                      ! Empty routine
   END  SUBROUTINE  trc_nam_pisces
#endif  

   !!======================================================================
END MODULE trcnam_pisces
