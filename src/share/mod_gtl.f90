!-------------------------------------------------------------------------------
!> Module generic tool
!!
!! @par Description
!!         This module is for the generic subroutine, e.g., global mean.
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_gtl
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: GTL_global_sum
  public :: GTL_global_sum_srf
  public :: GTL_global_sum_eachlayer
  public :: GTL_global_mean
  public :: GTL_global_mean_eachlayer
  public :: GTL_max
  public :: GTL_max_k
  public :: GTL_min
  public :: GTL_min_k

  interface GTL_global_sum
     module procedure GTL_global_sum_SP
     module procedure GTL_global_sum_DP
  end interface GTL_global_sum

  interface GTL_global_sum_srf
     module procedure GTL_global_sum_srf_SP
     module procedure GTL_global_sum_srf_DP
  end interface GTL_global_sum_srf

  interface GTL_max
     module procedure GTL_max_SP
     module procedure GTL_max_DP
  end interface GTL_max

  interface GTL_min
     module procedure GTL_min_SP
     module procedure GTL_min_DP
  end interface GTL_min

  public :: GTL_mk_rigidrotation

  public :: GTL_clip_region
  public :: GTL_clip_region2
  public :: GTL_clip_region_1layer
  public :: GTL_clip_region_1layer_k

  public :: GTL_clip_region_SP
  public :: GTL_clip_region_DP
  public :: GTL_clip_region2_SP
  public :: GTL_clip_region2_DP
  public :: GTL_clip_region_1layer_SP
  public :: GTL_clip_region_1layer_DP
  public :: GTL_clip_region_1layer_k_SP
  public :: GTL_clip_region_1layer_k_DP

  interface GTL_clip_region
     module procedure GTL_clip_region_SP
     module procedure GTL_clip_region_DP
  end interface GTL_clip_region

  interface GTL_clip_region2
     module procedure GTL_clip_region2_SP
     module procedure GTL_clip_region2_DP
  end interface GTL_clip_region2

  interface GTL_clip_region_1layer
     module procedure GTL_clip_region_1layer_SP
     module procedure GTL_clip_region_1layer_DP
  end interface GTL_clip_region_1layer

  interface GTL_clip_region_1layer_k
     module procedure GTL_clip_region_1layer_k_SP
     module procedure GTL_clip_region_1layer_k_DP
  end interface GTL_clip_region_1layer_k

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  function GTL_global_sum_SP( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl, &
       ADM_kmin,    &
       ADM_kmax
    use mod_comm, only: &
       COMM_Stat_sum
    use mod_vmtr, only: &
       VMTR_VOLUME,   &
       VMTR_VOLUME_pl
    implicit none

    real(SP), intent(in) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(SP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(SP)             :: sum_g

    real(SP) :: sum
    integer  :: i, j, k, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    sum = 0.0_SP
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do k = ADM_kmin, ADM_kmax
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       sum = sum + var(suf(i,j),k,l) * VMTR_VOLUME(suf(i,j),k,l)
    enddo
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,        ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
          sum = sum + var_pl(ADM_gslf_pl,k,l) * VMTR_VOLUME_pl(ADM_gslf_pl,k,l)
       enddo
       enddo
    endif

    call COMM_Stat_sum( sum, sum_g )

    return
  end function GTL_global_sum_SP

  !-----------------------------------------------------------------------------
  function GTL_global_sum_DP( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl, &
       ADM_kmin,    &
       ADM_kmax
    use mod_comm, only: &
       COMM_Stat_sum
    use mod_vmtr, only: &
       VMTR_VOLUME,   &
       VMTR_VOLUME_pl
    implicit none

    real(DP), intent(in) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(DP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(DP)             :: sum_g

    real(DP) :: sum
    integer  :: i, j, k, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    sum = 0.0_DP
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do k = ADM_kmin, ADM_kmax
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       sum = sum + var(suf(i,j),k,l) * VMTR_VOLUME(suf(i,j),k,l)
    enddo
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,        ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
          sum = sum + var_pl(ADM_gslf_pl,k,l) * VMTR_VOLUME_pl(ADM_gslf_pl,k,l)
       enddo
       enddo
    endif

    call COMM_Stat_sum( sum, sum_g )

    return
  end function GTL_global_sum_DP

  !-----------------------------------------------------------------------------
  function GTL_global_sum_srf_SP( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_comm, only: &
       COMM_Stat_sum
    use mod_gmtr, only: &
       GMTR_area,   &
       GMTR_area_pl
    implicit none

    real(SP), intent(in) :: var   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(SP), intent(in) :: var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(SP)             :: sum_g

    real(SP) :: sum
    integer  :: i, j, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    sum = 0.0_SP
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       sum = sum + var(suf(i,j),ADM_KNONE,l) * GMTR_area(suf(i,j),l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          sum = sum + var_pl(ADM_gslf_pl,ADM_KNONE,l) * GMTR_area_pl(ADM_gslf_pl,l)
       enddo
    endif

    call COMM_Stat_sum( sum, sum_g )

    return
  end function GTL_global_sum_srf_SP

  !-----------------------------------------------------------------------------
  function GTL_global_sum_srf_DP( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_comm, only: &
       COMM_Stat_sum
    use mod_gmtr, only: &
       GMTR_area,   &
       GMTR_area_pl
    implicit none

    real(DP), intent(in) :: var   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(DP), intent(in) :: var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(DP)             :: sum_g

    real(DP) :: sum
    integer  :: i, j, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    sum = 0.0_DP
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       sum = sum + var(suf(i,j),ADM_KNONE,l) * GMTR_area(suf(i,j),l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          sum = sum + var_pl(ADM_gslf_pl,ADM_KNONE,l) * GMTR_area_pl(ADM_gslf_pl,l)
       enddo
    endif

    call COMM_Stat_sum( sum, sum_g )

    return
  end function GTL_global_sum_srf_DP

  !-----------------------------------------------------------------------------
  subroutine GTL_global_sum_eachlayer( var, var_pl, sum_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_comm, only: &
       COMM_Stat_sum_eachlayer
    use mod_gmtr, only: &
       GMTR_area,   &
       GMTR_area_pl
    use mod_vmtr, only: &
       VMTR_RGAM,   &
       VMTR_RGAM_pl
    implicit none

    real(RP), intent(in)  :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: sum_g (ADM_kall)

    real(RP) :: sum(ADM_kall)
    integer  :: i, j, k, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    sum(:) = 0.0_RP
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do k = 1,        ADM_kall
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       sum(k) = sum(k) + var      (suf(i,j),k,l)    &
                       * GMTR_area(suf(i,j),l)      &
                       / VMTR_RGAM(suf(i,j),k,l)**2
    enddo
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          sum(k) = sum(k) + var_pl      (ADM_gslf_pl,k,l)    &
                          * GMTR_area_pl(ADM_gslf_pl,l)      &
                          / VMTR_RGAM_pl(ADM_gslf_pl,k,l)**2
       enddo
       enddo
    endif

    call COMM_Stat_sum_eachlayer( ADM_kall, sum(:), sum_g(:) )

    return
  end subroutine GTL_global_sum_eachlayer

  !-----------------------------------------------------------------------------
  function GTL_global_mean( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    implicit none

    real(RP), intent(in) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP)             :: sum_g

    real(RP)       :: one   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP)       :: one_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    logical,  save :: first = .true.
    real(RP), save :: volume_g
    !---------------------------------------------------------------------------

    if ( first ) then
       !--- calc global volume at first time
       one   (:,:,:) = 1.0_RP
       one_pl(:,:,:) = 1.0_RP

       volume_g = GTL_global_sum( one(:,:,:), one_pl(:,:,:) )

       first = .false.
    endif

    sum_g = GTL_global_sum( var(:,:,:), var_pl(:,:,:) )

    sum_g = sum_g / volume_g

    return
  end function GTL_global_mean

  !-----------------------------------------------------------------------------
  subroutine GTL_global_mean_eachlayer( var, var_pl, sum_g )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    implicit none

    real(RP), intent(in)  :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: sum_g (ADM_kall)

    real(RP)       :: one   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP)       :: one_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP)       :: area_g(ADM_kall)
    !---------------------------------------------------------------------------

    one   (:,:,:) = 1.0_RP
    one_pl(:,:,:) = 1.0_RP

    call GTL_global_sum_eachlayer( one(:,:,:), one_pl(:,:,:), area_g(:) )
    call GTL_global_sum_eachlayer( var(:,:,:), var_pl(:,:,:), sum_g (:) )

    sum_g(:) = sum_g(:) / area_g(:)

    return
  end subroutine GTL_global_mean_eachlayer

  !-----------------------------------------------------------------------------
  function GTL_max_SP( var, var_pl, kdim, kstart, kend ) result( vmax_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_const, only: &
       CONST_HUGE
    use mod_comm, only: &
       COMM_Stat_max
    implicit none

    integer,  intent(in) :: kdim
    integer,  intent(in) :: kstart
    integer,  intent(in) :: kend
    real(SP), intent(in) :: var   (ADM_gall,   kdim,ADM_lall   )
    real(SP), intent(in) :: var_pl(ADM_gall_pl,kdim,ADM_lall_pl)
    real(SP)             :: vmax_g

    real(SP) :: vmax
    integer  :: i, j, k, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    vmax = -CONST_HUGE
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do k = kstart,   kend
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       vmax = max( vmax, var(suf(i,j),k,l) )
    enddo
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,      ADM_lall_pl
       do k = kstart, kend
          vmax = max( vmax, var_pl(ADM_gslf_pl,k,l) )
       enddo
       enddo
    endif

    call COMM_Stat_max( vmax, vmax_g )

    return
  end function GTL_max_SP

  !-----------------------------------------------------------------------------
  function GTL_max_DP( var, var_pl, kdim, kstart, kend ) result( vmax_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_const, only: &
       CONST_HUGE
    use mod_comm, only: &
       COMM_Stat_max
    implicit none

    integer,  intent(in) :: kdim
    integer,  intent(in) :: kstart
    integer,  intent(in) :: kend
    real(DP), intent(in) :: var   (ADM_gall,   kdim,ADM_lall   )
    real(DP), intent(in) :: var_pl(ADM_gall_pl,kdim,ADM_lall_pl)
    real(DP)             :: vmax_g

    real(DP) :: vmax
    integer  :: i, j, k, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    vmax = -CONST_HUGE
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do k = kstart,   kend
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       vmax = max( vmax, var(suf(i,j),k,l) )
    enddo
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,      ADM_lall_pl
       do k = kstart, kend
          vmax = max( vmax, var_pl(ADM_gslf_pl,k,l) )
       enddo
       enddo
    endif

    call COMM_Stat_max( vmax, vmax_g )

    return
  end function GTL_max_DP

  !-----------------------------------------------------------------------------
  function GTL_max_k( var, var_pl, k ) result( vmax_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_const, only: &
       CONST_HUGE
    use mod_comm, only: &
       COMM_Stat_max
    implicit none

    real(RP), intent(in) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer,  intent(in) :: k
    real(RP)             :: vmax_g

    real(RP) :: vmax
    integer  :: i, j, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    vmax = -CONST_HUGE
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       vmax = max( vmax, var(suf(i,j),k,l) )
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,        ADM_lall_pl
          vmax = max( vmax, var_pl(ADM_gslf_pl,k,l) )
       enddo
    endif

    call COMM_Stat_max( vmax, vmax_g )

    return
  end function GTL_max_k

  !-----------------------------------------------------------------------------
  function GTL_min_SP( var, var_pl, kdim, kstart, kend, nonzero ) result( vmin_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_const, only: &
       CONST_HUGE
    use mod_comm, only: &
       COMM_Stat_min
    implicit none

    integer,  intent(in) :: kdim
    integer,  intent(in) :: kstart
    integer,  intent(in) :: kend
    real(SP), intent(in) :: var   (ADM_gall,   kdim,ADM_lall   )
    real(SP), intent(in) :: var_pl(ADM_gall_pl,kdim,ADM_lall_pl)
    real(SP)             :: vmin_g

    logical,  intent(in), optional :: nonzero

    real(SP) :: vmin
    integer  :: i, j, k, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    if ( present(nonzero) ) then
       if ( nonzero ) then

          vmin = CONST_HUGE
!OCL SERIAL
          do l = 1,        ADM_lall
!OCL PARALLEL
          do k = kstart, kend
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             if (       var(suf(i,j),k,l) > 0.0_SP &
                  .AND. var(suf(i,j),k,l) < vmin ) then

                vmin = var(suf(i,j),k,l)

             endif
          enddo
          enddo
          enddo
          enddo

       if ( ADM_have_pl ) then
!OCL SERIAL
          do l = 1,        ADM_lall
!OCL PARALLEL
          do k = kstart, kend
             if (       var_pl(ADM_gslf_pl,k,l) > 0.0_SP &
                  .AND. var_pl(ADM_gslf_pl,k,l) < vmin ) then

                vmin = var_pl(ADM_gslf_pl,k,l)

             endif
          enddo
          enddo
       endif

       call COMM_Stat_min( vmin, vmin_g )
       return

       endif
    endif

    vmin = CONST_HUGE
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do k = kstart,   kend
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       vmin = min( vmin, var(suf(i,j),k,l) )
    enddo
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,      ADM_lall_pl
       do k = kstart, kend
          vmin = min( vmin, var_pl(ADM_gslf_pl,k,l) )
       enddo
       enddo
    endif

    call COMM_Stat_min( vmin, vmin_g )

    return
  end function GTL_min_SP

  !-----------------------------------------------------------------------------
  function GTL_min_DP( var, var_pl, kdim, kstart, kend, nonzero ) result( vmin_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_const, only: &
       CONST_HUGE
    use mod_comm, only: &
       COMM_Stat_min
    implicit none

    integer,  intent(in) :: kdim
    integer,  intent(in) :: kstart
    integer,  intent(in) :: kend
    real(DP), intent(in) :: var   (ADM_gall,   kdim,ADM_lall   )
    real(DP), intent(in) :: var_pl(ADM_gall_pl,kdim,ADM_lall_pl)
    real(DP)             :: vmin_g

    logical,  intent(in), optional :: nonzero

    real(DP) :: vmin
    integer  :: i, j, k, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    if ( present(nonzero) ) then
       if ( nonzero ) then

          vmin = CONST_HUGE
!OCL SERIAL
          do l = 1,        ADM_lall
!OCL PARALLEL
          do k = kstart, kend
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             if (       var(suf(i,j),k,l) > 0.0_DP &
                  .AND. var(suf(i,j),k,l) < vmin ) then

                vmin = var(suf(i,j),k,l)

             endif
          enddo
          enddo
          enddo
          enddo

       if ( ADM_have_pl ) then
!OCL SERIAL
          do l = 1,        ADM_lall
!OCL PARALLEL
          do k = kstart, kend
             if (       var_pl(ADM_gslf_pl,k,l) > 0.0_DP &
                  .AND. var_pl(ADM_gslf_pl,k,l) < vmin ) then

                vmin = var_pl(ADM_gslf_pl,k,l)

             endif
          enddo
          enddo
       endif

       call COMM_Stat_min( vmin, vmin_g )
       return

       endif
    endif

    vmin = CONST_HUGE
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do k = kstart,   kend
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       vmin = min( vmin, var(suf(i,j),k,l) )
    enddo
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,      ADM_lall_pl
       do k = kstart, kend
          vmin = min( vmin, var_pl(ADM_gslf_pl,k,l) )
       enddo
       enddo
    endif

    call COMM_Stat_min( vmin, vmin_g )

    return
  end function GTL_min_DP

  !-----------------------------------------------------------------------------
  function GTL_min_k( var, var_pl, k ) result( vmin_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_1d, &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_const, only: &
       CONST_HUGE
    use mod_comm, only: &
       COMM_Stat_min
    implicit none

    real(RP), intent(in) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer,  intent(in) :: k
    real(RP)             :: vmin_g

    real(RP) :: vmin
    integer  :: i, j, l

    integer  :: suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    vmin = +CONST_HUGE
!OCL SERIAL
    do l = 1,        ADM_lall
!OCL PARALLEL
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       vmin = min( vmin, var(suf(i,j),k,l) )
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,        ADM_lall_pl
          vmin = min( vmin, var_pl(ADM_gslf_pl,k,l) )
       enddo
    endif

    call COMM_Stat_min( vmin, vmin_g )

    return
  end function GTL_min_k

  !-----------------------------------------------------------------------------
  subroutine GTL_mk_rigidrotation( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl, &
       alpha,     &
       vmax       )
    use mod_adm, only: &
       ADM_KNONE,       &
       ADM_have_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_kall,        &
       ADM_gslf_pl
    use mod_grd, only: &
       GRD_LAT,    &
       GRD_LAT_pl, &
       GRD_LON,    &
       GRD_LON_pl
    use mod_gmtr, only: &
       P_IX => GMTR_p_IX, &
       P_IY => GMTR_p_IY, &
       P_IZ => GMTR_p_IZ, &
       P_JX => GMTR_p_JX, &
       P_JY => GMTR_p_JY, &
       P_JZ => GMTR_p_JZ, &
       GMTR_p,            &
       GMTR_p_pl
    implicit none

    real(RP), intent(inout) :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: alpha
    real(RP), intent(in)    :: vmax

    real(RP) :: u, v

    integer  :: g, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       u =  vmax * ( cos(GRD_LAT(g,l)) * cos(alpha) &
                   + sin(GRD_LAT(g,l))              &
                   * cos(GRD_LON(g,l)) * sin(alpha) )
       v = -vmax * ( sin(GRD_LON(g,l)) * sin(alpha) )

       vx(g,k,l) = u * GMTR_p(g,k0,l,P_IX) &
                 + v * GMTR_p(g,k0,l,P_JX)
       vy(g,k,l) = u * GMTR_p(g,k0,l,P_IY) &
                 + v * GMTR_p(g,k0,l,P_JY)
       vz(g,k,l) = u * GMTR_p(g,k0,l,P_IZ) &
                 + v * GMTR_p(g,k0,l,P_JZ)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       g = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          u =  vmax * ( cos(GRD_LAT_pl(g,l)) * cos(alpha) &
                      + sin(GRD_LAT_pl(g,l))              &
                      * cos(GRD_LON_pl(g,l)) * sin(alpha) )
          v = -vmax * ( sin(GRD_LON_pl(g,l)) * sin(alpha) )

          vx_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,P_IX) &
                       + v * GMTR_p_pl(g,k0,l,P_JX)
          vy_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,P_IY) &
                       + v * GMTR_p_pl(g,k0,l,P_JY)
          vz_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,P_IZ) &
                       + v * GMTR_p_pl(g,k0,l,P_JZ)
       enddo
       enddo
    endif

    return
  end subroutine GTL_mk_rigidrotation

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_SP( var_orig, var_clip, kmin, kmax )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gall_1d, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(RP), intent(in)  :: var_orig(ADM_gall   ,ADM_kall       ,ADM_lall)
    real(SP), intent(out) :: var_clip(ADM_gall_in,1:(kmax-kmin+1),ADM_lall)

    integer  :: gmin, gmax, lall, gall_1d
    integer  :: i, j, k, l, n, ij
    !---------------------------------------------------------------------------

    gmin    = ADM_gmin
    gmax    = ADM_gmax
    gall_1d = ADM_gall_1d
    lall    = ADM_lall

!OCL XFILL
    !$omp parallel do default(none),private(i,j,k,l,n,ij), &
    !$omp shared(gmin,gmax,gall_1d,lall,kmin,kmax,var_clip,var_orig), &
    !$omp collapse(3)
    do l = 1,    lall
    do k = kmin, kmax
    do j = gmin, gmax+1
    do i = gmin, gmax+1
       n  = (gall_1d-1) * (j-2) + (i-1)
       ij = (gall_1d  ) * (j-1) + (i  )

       var_clip(n,k-kmin+1,l) = var_orig(ij,k,l)
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine GTL_clip_region_SP

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_DP( var_orig, var_clip, kmin, kmax )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gall_1d, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(RP), intent(in)  :: var_orig(ADM_gall   ,ADM_kall       ,ADM_lall)
    real(DP), intent(out) :: var_clip(ADM_gall_in,1:(kmax-kmin+1),ADM_lall)

    integer  :: gmin, gmax, lall, gall_1d
    integer  :: i, j, k, l, n, ij
    !---------------------------------------------------------------------------

    gmin    = ADM_gmin
    gmax    = ADM_gmax
    gall_1d = ADM_gall_1d
    lall    = ADM_lall

!OCL XFILL
    !$omp parallel do default(none),private(i,j,k,l,n,ij), &
    !$omp shared(gmin,gmax,gall_1d,lall,kmin,kmax,var_clip,var_orig), &
    !$omp collapse(3)
    do l = 1,    lall
    do k = kmin, kmax
    do j = gmin, gmax+1
    do i = gmin, gmax+1
       n  = (gall_1d-1) * (j-2) + (i-1)
       ij = (gall_1d  ) * (j-1) + (i  )

       var_clip(n,k-kmin+1,l) = var_orig(ij,k,l)
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine GTL_clip_region_DP

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region2_SP( var_orig, var_clip, kmin, kmax, vmax )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gall_1d, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    integer,  intent(in)  :: vmax
    real(RP), intent(in)  :: var_orig(ADM_gall   ,ADM_kall       ,vmax,ADM_lall)
    real(SP), intent(out) :: var_clip(ADM_gall_in,1:(kmax-kmin+1),vmax,ADM_lall)

    integer  :: gmin, gmax, lall, gall_1d
    integer  :: i, j, k, v, l, n, ij
    !---------------------------------------------------------------------------

    gmin    = ADM_gmin
    gmax    = ADM_gmax
    gall_1d = ADM_gall_1d
    lall    = ADM_lall

!OCL XFILL
    !$omp parallel do default(none),private(i,j,k,v,l,n,ij), &
    !$omp shared(gmin,gmax,gall_1d,lall,kmin,kmax,vmax,var_clip,var_orig), &
    !$omp collapse(4)
    do l = 1,    lall
    do v = 1,    vmax
    do k = kmin, kmax
    do j = gmin, gmax+1
    do i = gmin, gmax+1
       n  = (gall_1d-1) * (j-2) + (i-1)
       ij = (gall_1d  ) * (j-1) + (i  )

       var_clip(n,k-kmin+1,v,l) = var_orig(ij,k,v,l)
    enddo
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine GTL_clip_region2_SP

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region2_DP( var_orig, var_clip, kmin, kmax, vmax )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gall_1d, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    integer,  intent(in)  :: vmax
    real(RP), intent(in)  :: var_orig(ADM_gall   ,ADM_kall       ,vmax,ADM_lall)
    real(DP), intent(out) :: var_clip(ADM_gall_in,1:(kmax-kmin+1),vmax,ADM_lall)

    integer  :: gmin, gmax, lall, gall_1d
    integer  :: i, j, k, v, l, n, ij
    !---------------------------------------------------------------------------

    gmin    = ADM_gmin
    gmax    = ADM_gmax
    gall_1d = ADM_gall_1d
    lall    = ADM_lall

!OCL XFILL
    !$omp parallel do default(none),private(i,j,k,v,l,n,ij), &
    !$omp shared(gmin,gmax,gall_1d,lall,kmin,kmax,vmax,var_clip,var_orig), &
    !$omp collapse(4)
    do l = 1,    lall
    do v = 1,    vmax
    do k = kmin, kmax
    do j = gmin, gmax+1
    do i = gmin, gmax+1
       n  = (gall_1d-1) * (j-2) + (i-1)
       ij = (gall_1d  ) * (j-1) + (i  )

       var_clip(n,k-kmin+1,v,l) = var_orig(ij,k,v,l)
    enddo
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine GTL_clip_region2_DP

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer_SP( var_orig, var_clip )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    real(RP), intent(in)  :: var_orig(ADM_gall   ,ADM_lall)
    real(SP), intent(out) :: var_clip(ADM_gall_in,ADM_lall)

    integer  :: gmin, gmax, lall, gall_1d
    integer  :: i, j, l, n, ij
    !---------------------------------------------------------------------------

    gmin    = ADM_gmin
    gmax    = ADM_gmax
    gall_1d = ADM_gall_1d
    lall    = ADM_lall

!OCL XFILL
    !$omp parallel do default(none),private(i,j,l,n,ij), &
    !$omp shared(gmin,gmax,gall_1d,lall,var_clip,var_orig), &
    !$omp collapse(2)
    do l = 1,    lall
    do j = gmin, gmax+1
    do i = gmin, gmax+1
       n  = (gall_1d-1) * (j-2) + (i-1)
       ij = (gall_1d  ) * (j-1) + (i  )

       var_clip(n,l) = var_orig(ij,l)
    enddo
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine GTL_clip_region_1layer_SP

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer_DP( var_orig, var_clip )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    real(RP), intent(in)  :: var_orig(ADM_gall   ,ADM_lall)
    real(DP), intent(out) :: var_clip(ADM_gall_in,ADM_lall)

    integer  :: gmin, gmax, lall, gall_1d
    integer  :: i, j, l, n, ij
    !---------------------------------------------------------------------------

    gmin    = ADM_gmin
    gmax    = ADM_gmax
    gall_1d = ADM_gall_1d
    lall    = ADM_lall

!OCL XFILL
    !$omp parallel do default(none),private(i,j,l,n,ij), &
    !$omp shared(gmin,gmax,gall_1d,lall,var_clip,var_orig), &
    !$omp collapse(2)
    do l = 1,    lall
    do j = gmin, gmax+1
    do i = gmin, gmax+1
       n  = (gall_1d-1) * (j-2) + (i-1)
       ij = (gall_1d  ) * (j-1) + (i  )

       var_clip(n,l) = var_orig(ij,l)
    enddo
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine GTL_clip_region_1layer_DP

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer_k_SP( var_orig, var_clip, ksize, k )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    integer,  intent(in)  :: ksize
    real(RP), intent(in)  :: var_orig(ADM_gall   ,ksize,ADM_lall)
    real(SP), intent(out) :: var_clip(ADM_gall_in,      ADM_lall)
    integer,  intent(in)  :: k

    integer  :: gmin, gmax, lall, gall_1d
    integer  :: i, j, l, n, ij
    !---------------------------------------------------------------------------

    gmin    = ADM_gmin
    gmax    = ADM_gmax
    gall_1d = ADM_gall_1d
    lall    = ADM_lall

!OCL XFILL
    !$omp parallel do default(none),private(i,j,l,n,ij), &
    !$omp shared(k,gmin,gmax,gall_1d,lall,var_clip,var_orig), &
    !$omp collapse(2)
    do l = 1,    lall
    do j = gmin, gmax+1
    do i = gmin, gmax+1
       n  = (gall_1d-1) * (j-2) + (i-1)
       ij = (gall_1d  ) * (j-1) + (i  )

       var_clip(n,l) = var_orig(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine GTL_clip_region_1layer_k_SP

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer_k_DP( var_orig, var_clip, ksize, k )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    integer,  intent(in)  :: ksize
    real(RP), intent(in)  :: var_orig(ADM_gall   ,ksize,ADM_lall)
    real(DP), intent(out) :: var_clip(ADM_gall_in,      ADM_lall)
    integer,  intent(in)  :: k

    integer  :: gmin, gmax, lall, gall_1d
    integer  :: i, j, l, n, ij
    !---------------------------------------------------------------------------

    gmin    = ADM_gmin
    gmax    = ADM_gmax
    gall_1d = ADM_gall_1d
    lall    = ADM_lall

!OCL XFILL
    !$omp parallel do default(none),private(i,j,l,n,ij), &
    !$omp shared(k,gmin,gmax,gall_1d,lall,var_clip,var_orig), &
    !$omp collapse(2)
    do l = 1,    lall
    do j = gmin, gmax+1
    do i = gmin, gmax+1
       n  = (gall_1d-1) * (j-2) + (i-1)
       ij = (gall_1d  ) * (j-1) + (i  )

       var_clip(n,l) = var_orig(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine GTL_clip_region_1layer_k_DP

end module mod_gtl
