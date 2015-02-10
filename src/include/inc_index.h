!-------------------------------------------------------------------------------
!> Pre-defined grid index
!!
!! @par Description
!!          Include file of grid index
!!          If you set environment "ENABLE_FIXEDINDEX=T", this file is used
!!          If you set this file, execute following:
!!          "make fixedindex glevel=xxx rlevel=xxx nmpi=xxx zlayer=xxx diamond=xxx".
!<
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Total inner grid (2D) = 10242
  ! Total inner grid (3D) = 962748
  ! Total grid       (2D) = 1156
  ! Total grid       (3D) = 1109760
  !-----------------------------------------------------------------------------

  ! main parameter
  integer, public, parameter :: ADM_glevel      = 5
  integer, public, parameter :: ADM_rlevel      = 0
  integer, public, parameter :: ADM_vlayer      = 94
  integer, public, parameter :: ADM_DMD         = 10
  integer, public            :: ADM_prc_all     = 10

  ! region
  integer, public, parameter :: ADM_rgn_nmax    = 10
  integer, public, parameter :: ADM_lall        = 1
  integer, public, parameter :: ADM_rgn_nmax_pl = 2 ! number of pole    region
  integer, public, parameter :: ADM_lall_pl     = 2 ! number of pole    region per process

  ! horizontal grid
  integer, public, parameter :: ADM_gall        = 1156
  integer, public, parameter :: ADM_gall_in     = 1089
  integer, public, parameter :: ADM_gall_1d     = 34
  integer, public, parameter :: ADM_gmin        = 2
  integer, public, parameter :: ADM_gmax        = 33

  integer, public, parameter :: ADM_gall_pl     = 6
  integer, public, parameter :: ADM_gslf_pl     = 1
  integer, public, parameter :: ADM_gmin_pl     = 2
  integer, public, parameter :: ADM_gmax_pl     = 6

  ! vertical grid
  integer, public, parameter :: ADM_kall        = 96
  integer, public, parameter :: ADM_kmin        = 2
  integer, public, parameter :: ADM_kmax        = 95

  ! List vectors
  integer, public, parameter :: ADM_IooJoo_nmax = 1024
  integer, public, parameter :: ADM_IooJmo_nmax = 1056
  integer, public, parameter :: ADM_IooJop_nmax = 1056
  integer, public, parameter :: ADM_IooJmp_nmax = 1088
  integer, public, parameter :: ADM_ImoJoo_nmax = 1056
  integer, public, parameter :: ADM_ImoJmo_nmax = 1089
  integer, public, parameter :: ADM_ImoJop_nmax = 1089
  integer, public, parameter :: ADM_ImoJmp_nmax = 1122
  integer, public, parameter :: ADM_IopJoo_nmax = 1056
  integer, public, parameter :: ADM_IopJmo_nmax = 1089
  integer, public, parameter :: ADM_IopJop_nmax = 1089
  integer, public, parameter :: ADM_IopJmp_nmax = 1122
  integer, public, parameter :: ADM_ImpJoo_nmax = 1088
  integer, public, parameter :: ADM_ImpJmo_nmax = 1122
  integer, public, parameter :: ADM_ImpJop_nmax = 1122
  integer, public, parameter :: ADM_ImpJmp_nmax = 1156

!-------------------------------------------------------------------------------
