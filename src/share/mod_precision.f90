!-------------------------------------------------------------------------------
!> module PRECISION
!!
!! @par Description
!!          precision module
!!          Imported from SCALE library
!!
!! @author Team SCALE
!!
!<
module mod_precision
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: SP = kind(0.E0) ! Single Precision: kind(0.E0)
  integer, public, parameter :: DP = kind(0.D0) ! Double Precision: kind(0.D0)

  integer, public, parameter :: SP_PREC = precision(0.E0)
  integer, public, parameter :: DP_PREC = precision(0.D0)

#ifdef SINGLE
  integer, public, parameter :: RP      = SP
  integer, public, parameter :: RP_PREC = SP_PREC
#else
  integer, public, parameter :: RP      = DP
  integer, public, parameter :: RP_PREC = DP_PREC
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
end module mod_precision
