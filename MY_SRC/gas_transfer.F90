MODULE gastransfer
   !!======================================================================
   !!                         ***  MODULE trcdms_medusa  ***
   !! TOP :   MEDUSA
   !!======================================================================
   !! History :
   !!  -   !  2015-06  (A. Yool)             added for UKESM1 project
   !!----------------------------------------------------------------------
#if defined key_skog
!#if defined key_medusa && defined key_roam
      USE oce_trc
      USE trc
      USE sms_pisces                   !  PISCES Source Minus Sink variables
      USE prtctl_trc                   !  print control for debugging
      USE iom                          !  I/O manager
      USE fldread                      !  read input fields
 !     USE lbclnk
 !     USE prtctl_trc      ! Print control for debugging
 !     USE in_out_manager  ! I/O manager

      IMPLICIT NONE

      PUBLIC   gas_transfer  ! called by trcbio_medusa.F90 module

   !!* Substitution
!#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

      subroutine gas_transfer(wind, N, eqn, kw660)
! --------------------------------------------------------------------
!  Gas transfer velocity
! --------------------------------------------------------------------
!
! Title  : Calculates gas transfer velocity
! Author : Andrew Yool
! Date   : 15/10/04
! 
! This subroutine uses near-surface wind speed to calculate gas
! transfer velocity for use in CO2 and O2 exchange calculations.
!
! Note that the parameterisation of Wanninkhof quoted here is a
! truncation of the original equation.  It excludes a chemical
! enhancement function (based on temperature), although such
! temperature dependence is reported negligible by Etcheto & 
! Merlivat (1988).
!
! Note also that in calculating scalar wind, the variance of the
! wind over the period of a timestep is ignored.  Some authors,
! for instance OCMIP-2, favour including some reference to the
! variability of wind.  However, their wind fields are averaged
! over relatively long time periods, and so this issue may be
! safely (!) ignored here.
! 
! AXY (12/06/2015)
! UPDATED: revised formulation from Wanninkhof (2014) update to
! original 1992 paper. Full reference is:
!
! Wanninkhof, R. (2014). Relationship between wind speed and gas 
! exchange over the ocean revisited. LIMNOLOGY AND OCEANOGRAPHY-METHODS
! 12, 351-362, doi:10.4319/lom.2014.12.351
! 
! Subroutine inputs are (in order) : 
!     wind      wind velocity at 10 m (m/s)
!     N         size of input array (value 1 at this time)
!     eqn       choice of parameterisation:
!               1 = Liss & Merlivat (1986)    [approximated]
!               2 = Wanninkhof (1992)         [sans enhancement]
!               3 = Nightingale et al. (2000) [good]
!               4 = Nightingale et al. (2000) [better]
!               5 = Nightingale et al. (2000) [best]
!               6 = OCMIP-2                   [sans variability]
!               7 = Wanninkhof (2014)         [assumes 6h avg winds]
! (*) k         gas transfer velocity (m/s)
! 
! Where (*) is the function output and (+) is a diagnostic output.
!
      implicit none

      INTEGER, INTENT(in) :: N, eqn
! Input variables
!     real(kind=wp), INTENT(in),  DIMENSION(N) :: wind
      real(wp), INTENT(in)  :: wind
!
! Output variables
!     real(kind=wp), INTENT(out), DIMENSION(N) :: kw660
      real(wp), INTENT(out) :: kw660
!
!     INTEGER :: eqn
!
! Coefficients for various parameterisations
      real(wp) :: a(7)
      real(wp) :: b(7)
!
!     real(wp), DIMENSION(N) :: tmp_k
      real(wp) :: tmp_k
!
! Values of coefficients
      data a(1) / 0.166 /  ! Liss & Merlivat (1986)    [approximated]
      data a(2) / 0.3   /  ! Wanninkhof (1992)         [sans enhancement]
      data a(3) / 0.23  /  ! Nightingale et al. (2000) [good]
      data a(4) / 0.23  /  ! Nightingale et al. (2000) [better]
      data a(5) / 0.222 /  ! Nightingale et al. (2000) [best]
      data a(6) / 0.337 /  ! OCMIP-2                   [sans variability]
      data a(7) / 0.251 /  ! Wanninkhof (2014)         [assumes 6h avg winds]
!
      data b(1) / 0.133 /
      data b(2) / 0.0   /
      data b(3) / 0.0   /
      data b(4) / 0.1   /
      data b(5) / 0.333 /
      data b(6) / 0.0   /
      data b(7) / 0.0   /
!
! Which parameterisation is to be used?
!     eqn = 7
!
! Calculate gas transfer velocity (cm/h)
      tmp_k = (a(eqn) * wind**2) + (b(eqn) * wind)
!
! Convert tmp_k from cm/h to m/s
      kw660 = tmp_k / (100. * 3600.)
!
      return

    end subroutine gas_transfer

!=======================================================================
!=======================================================================
!=======================================================================
#endif
   !!======================================================================
END MODULE gastransfer
