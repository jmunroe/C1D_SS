MODULE p4zsbc
   !!======================================================================
   !!                         ***  MODULE p4sbc  ***
   !! TOP :   SMELT surface boundary conditions on nutrients
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, C. Ethe) Original code
   !!                     2015 (E. Olson) adapted handling of river tracers 
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sbc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_sbc_init   :  Initialization of p4z_sbc
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  SMELT Source Minus Sink variables
   USE iom             !  I/O manager
   USE fldread         !  time interpolation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sbc
   PUBLIC   p4z_sbc_init   

   !! * Shared module variables
   LOGICAL , PUBLIC  :: ln_river    !: boolean for river input of nutrients
   LOGICAL , PUBLIC  :: ln_turb     !: boolean for setting river turbidity
   TYPE(FLD_N)       :: sn_turb     !! information about river turbidity file
 
   LOGICAL , PUBLIC  :: ll_sbc

   INTEGER , PUBLIC, PARAMETER  :: jp_nriv  = jp_pisces  !: Maximum number of river input fields: treat all pisces vars

   !! * Module variables
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_river  ! structure of input river nutrients
   !$AGRIF_DO_NOT_TREAT
   TYPE(FLD), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   sf_turb ! structure: river runoff (file information, fields read)
   !$AGRIF_END_DO_NOT_TREAT
   INTEGER , PARAMETER :: nbtimes = 365  !: maximum number of times record in a file

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: rnf_pis, rnf_pis_b    !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: rnf_turb ! river turbidity field

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsbc.F90 5507 2015-06-29 15:19:49Z aumont $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sbc  ***
      !!
      !! ** purpose :   read and interpolate the external sources of nutrients
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * local declarations
      INTEGER  :: jn
      REAL(wp) :: znfact
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sbc')

      ! Handling of rivers adapted to SMELT
      ! Compute river at nit000 or only if there is more than 1 time record in river file
      ! -----------------------------------------
      IF( ln_river ) THEN
         IF( kt /= nit000) THEN
            rnf_pis_b(:,:,:) = rnf_pis(:,:,:) ! swap 
         ENDIF
         CALL fld_read( kt, nn_fsbc, sf_river )
         IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN
           DO jn = jp_pcs0, jp_pcs1
                znfact=1.0_wp
                rnf_pis(:,:,jn) = ( sf_river(jn)%fnow(:,:,1) ) * rnf(:,:) * r1_rau0 * znfact
              WHERE (sf_river(jn)%fnow(:,:,1) < -9._wp) 
                 rnf_pis(:,:,jn) = trn(:,:,1,jn) * rnf(:,:) * r1_rau0
              END WHERE
           END DO
         ENDIF
         IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
            IF( ln_rsttr .AND.    &                               !* Restart: read in restart file
                             & iom_varid( numrtr, 'rnf_pis_'//TRIM(ctrcnm(1))//'_b', ldstop = .FALSE. ) > 0 ) THEN
                IF(lwp) WRITE(numout,*) ' p4zsbc      nit000-1 runoff forcing fields read in the restart file'
                DO jn = 1, jp_nriv                          ! before runoff !
                    CALL iom_get( numrtr, jpdom_autoglo, 'rnf_pis_'//TRIM(ctrcnm(jn))//'_b', rnf_pis_b(:,:,jn) )
                END DO
            ELSE                                                   !* no restart: set from nit000 values
                IF(lwp) WRITE(numout,*) ' p4zsbc      nit000-1 runoff forcing fields set to nit000'
                rnf_pis_b(:,:,:) = rnf_pis(:,:,:)
            ENDIF
         ENDIF
      ENDIF
      ! set turbidity in Fraser based on input file
      IF( ln_turb ) THEN
         CALL fld_read( kt, nn_fsbc, sf_turb )
      ENDIF
      IF( lrst_trc ) THEN
      ! EO: save rnf_pis_b in restart file (save current rnf_pis since this will be _b at restart)
         IF (ln_river) THEN
             DO jn = 1, jp_nriv
                 CALL iom_rstput( kt, nitrst, numrtw, 'rnf_pis_'//TRIM(ctrcnm(jn))//'_b', rnf_pis(:,:,jn))
             END DO
         ENDIF
      ENDIF

      IF( nn_timing == 1 )  CALL timing_stop('p4z_sbc')
      !
   END SUBROUTINE p4z_sbc

   SUBROUTINE p4z_sbc_init

      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sbc_init  ***
      !!
      !! ** purpose :   initialization of the external sources of nutrients
      !!
      !! ** method  :   read the files and compute the budget
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER  :: ifpr
      INTEGER  :: ierr, ierr1
      INTEGER  :: ios                 ! Local integer output status for namelist read
      !
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION(jp_nriv) ::  slf_river    ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_riverno3, sn_rivernh4, sn_riversil   ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_riverdia, sn_riverphy, sn_rivermyr, sn_riverzoo
      TYPE(FLD_N) ::   sn_riverdoc, sn_riverpoc, sn_riverdsi, sn_rivertra
#if defined key_oxy
      TYPE(FLD_N) ::   sn_rivero2
#endif
#if defined key_skog
      TYPE(FLD_N) ::   sn_riverdic, sn_riverta
#endif
      !
#if defined key_skog && defined key_oxy
      NAMELIST/nampissbc/ cn_dir, sn_riverno3, sn_rivernh4, sn_riversil, sn_riverdia,    &
           &                sn_riverphy, sn_rivermyr, sn_riverzoo, sn_riverdoc, &
        &                   sn_riverpoc, sn_riverdsi,   &
        &                sn_rivertra, sn_turb, ln_river, ln_turb, &
        &                sn_riverdic, sn_riverta, sn_rivero2
#elif defined key_oxy
      NAMELIST/nampissbc/ cn_dir, sn_riverno3, sn_rivernh4, sn_riversil, sn_riverdia,    &
           &                sn_riverphy, sn_rivermyr, sn_riverzoo, sn_riverdoc, &
        &                   sn_riverpoc, sn_riverdsi,   &
        &                sn_rivertra, sn_turb, ln_river, ln_turb, &
        &                sn_rivero2
#else
      NAMELIST/nampissbc/ cn_dir, sn_riverno3, sn_rivernh4, sn_riversil, sn_riverdia,    &
           &                sn_riverphy, sn_rivermyr, sn_riverzoo, sn_riverdoc, &
           &                sn_riverpoc, sn_riverdsi,   &
        &                sn_rivertra, sn_turb, ln_river, ln_turb
#endif
      
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sbc_init')
      !
      !                            !* set file information
      REWIND( numnatp_ref )              ! Namelist nampissbc in reference namelist : Pisces external sources of nutrients
      READ  ( numnatp_ref, nampissbc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissbc in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampissbc in configuration namelist : Pisces external sources of nutrients
      READ  ( numnatp_cfg, nampissbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampissbc )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist : nampissbc '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '    river input of nutrients                 ln_river    = ', ln_river
         WRITE(numout,*) '    river turbidity                          ln_turb     = ', ln_turb
            END IF
      IF ( ln_river .OR. ln_turb) ll_sbc=.True.
      IF ( ln_turb .AND. .NOT. ln_rnf) THEN
         CALL ctl_stop( 'STOP', 'p4z_sbc_init: cannot have river turbidity without rivers' )
      ENDIF

      ! set the number of level over which river runoffs are applied 
      ! online configuration : computed in sbcrnf
      IF( lk_offline ) THEN
        nk_rnf(:,:) = 1
        h_rnf (:,:) = fsdept(:,:,1)
      ENDIF

      ! nutrient input from rivers
      ! --------------------------
      IF( ln_river ) THEN
         !
         slf_river(jpno3) = sn_riverno3  ;  slf_river(jpnh4) = sn_rivernh4  ;  slf_river(jpsil) = sn_riversil 
         slf_river(jpdia) = sn_riverdia  ;  slf_river(jpphy) = sn_riverphy  ;  slf_river(jpmyr) = sn_rivermyr
         slf_river(jpzoo) = sn_riverzoo  ;  slf_river(jpdon) = sn_riverdoc  ;  slf_river(jppon) = sn_riverpoc  
         slf_river(jpdsi) = sn_riverdsi  ;  slf_river(jpriv) = sn_rivertra
#if defined key_skog
         slf_river(jpdic) = sn_riverdic  ;  slf_river(jpta) = sn_riverta 
#endif
#if defined key_oxy
         slf_river(jpo2) = sn_rivero2
#endif          
         !
         ALLOCATE( rnf_pis(jpi,jpj,jp_nriv), rnf_pis_b(jpi,jpj,jp_nriv) ) 
         ALLOCATE( sf_river(jp_nriv), STAT=ierr1 )           !* allocate and fill sf_river with sn_river_

         IF( ierr1 > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_irver structure' )
         !
         CALL fld_fill( sf_river, slf_river, cn_dir, 'p4z_sed_init', 'Input from river ', 'nampissed' )
         DO ifpr = 1, jp_nriv
                                          ALLOCATE( sf_river(ifpr)%fnow(jpi,jpj,1  ) )
            IF( slf_river(ifpr)%ln_tint ) ALLOCATE( sf_river(ifpr)%fdta(jpi,jpj,1,2) )
         
         END DO
      END IF 
      ! river turbidity
      ! ---------------
      IF( ln_turb) THEN
         ALLOCATE( sf_turb(1), STAT=ierr )   ! create sf_turb structure (river turbidity)
 
         IF ( ierr > 0 ) THEN
            CALL ctl_stop( 'p4zsbc: unable to allocate sf_turb structure') ;  RETURN
         ENDIF
         CALL fld_fill(sf_turb, (/ sn_turb /), cn_dir, 'p4z_sbc_init','read river turbidity', 'nampisriv') 
         ALLOCATE( sf_turb(1)%fnow(jpi,jpj,1)  )
         IF(  sn_turb%ln_tint ) ALLOCATE( sf_turb(1)%fdta(jpi,jpj,1,2) )
         ENDIF
      IF( ll_sbc ) CALL p4z_sbc( nit000 ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sbc_init')
      !
   END SUBROUTINE p4z_sbc_init

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sbc                         ! Empty routine
   END SUBROUTINE p4z_sbc
#endif 
   !!======================================================================
END MODULE p4zsbc
