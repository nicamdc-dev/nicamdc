!-------------------------------------------------------------------------------
!> Module generic tool
!!
!! @par Description
!!         This module is for the generic subroutine, e.g., global mean.
!!
!! @author NICAM developers, Team SCALE
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

  public :: GTL_mk_rigidrotation

  public :: GTL_clip_region
  public :: GTL_clip_region_1layer
  public :: GTL_clip_region_1layer_k

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
  function GTL_global_sum( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
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

    real(RP), intent(in) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP)             :: sum_g

    real(RP) :: sum
    integer  :: i, j, k, l
    !---------------------------------------------------------------------------

    sum = 0.0_RP
    do l = 1,        ADM_lall
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
  end function GTL_global_sum

  !-----------------------------------------------------------------------------
  function GTL_global_sum_srf( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
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

    real(RP), intent(in) :: var   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP)             :: sum_g

    real(RP) :: sum
    integer  :: i, j, l
    !---------------------------------------------------------------------------

    sum = 0.0_RP
    do l = 1, ADM_lall
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
  end function GTL_global_sum_srf

  !-----------------------------------------------------------------------------
  subroutine GTL_global_sum_eachlayer( var, var_pl, sum_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
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
    !---------------------------------------------------------------------------

    sum(:) = 0.0_RP
    do l = 1, ADM_lall
    do k = 1, ADM_kall
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
  function GTL_max( var, var_pl, kdim, kstart, kend ) result( vmax_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_comm, only: &
       COMM_Stat_max
    implicit none

    integer,  intent(in) :: kdim
    integer,  intent(in) :: kstart
    integer,  intent(in) :: kend
    real(RP), intent(in) :: var   (ADM_gall,   kdim,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,kdim,ADM_lall_pl)
    real(RP)             :: vmax_g

    real(RP) :: vmax
    integer  :: i, j, k, l
    !---------------------------------------------------------------------------

    vmax = -1.E+30_RP
    do l = 1,      ADM_lall
    do k = kstart, kend
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
  end function GTL_max

  !-----------------------------------------------------------------------------
  function GTL_max_k( var, var_pl, k ) result( vmax_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_comm, only: &
       COMM_Stat_max
    implicit none

    real(RP), intent(in) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer,  intent(in) :: k
    real(RP)             :: vmax_g

    real(RP) :: vmax
    integer  :: i, j, l
    !---------------------------------------------------------------------------

    vmax = -1.E+30_RP
    do l = 1,        ADM_lall
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
  function GTL_min( var, var_pl, kdim, kstart, kend, nonzero ) result( vmin_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_comm, only: &
       COMM_Stat_min
    implicit none

    integer,  intent(in) :: kdim
    integer,  intent(in) :: kstart
    integer,  intent(in) :: kend
    real(RP), intent(in) :: var   (ADM_gall,   kdim,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,kdim,ADM_lall_pl)
    real(RP)             :: vmin_g

    logical,  intent(in), optional :: nonzero

    real(RP) :: vmin
    integer  :: i, j, k, l
    !---------------------------------------------------------------------------

    if ( present(nonzero) ) then
       if ( nonzero ) then

          vmin = 1.E+30_RP
          do l = 1,      ADM_lall
          do k = kstart, kend
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             if (       var(suf(i,j),k,l) > 0.0_RP &
                  .AND. var(suf(i,j),k,l) < vmin ) then

                vmin = var(suf(i,j),k,l)

             endif
          enddo
          enddo
          enddo
          enddo

       if ( ADM_have_pl ) then
          do l = 1,      ADM_lall_pl
          do k = kstart, kend
             if (       var_pl(ADM_gslf_pl,k,l) > 0.0_RP &
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

    vmin = 1.E+30_RP
    do l = 1,      ADM_lall
    do k = kstart, kend
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
  end function GTL_min

  !-----------------------------------------------------------------------------
  function GTL_min_k( var, var_pl, k ) result( vmin_g )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl
    use mod_comm, only: &
       COMM_Stat_min
    implicit none

    real(RP), intent(in) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer,  intent(in) :: k
    real(RP)             :: vmin_g

    real(RP) :: vmin
    integer  :: i, j, l
    !---------------------------------------------------------------------------

    vmin = +1.E+30_RP
    do l = 1,        ADM_lall
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
       GRD_LAT, &
       GRD_LON, &
       GRD_s,   &
       GRD_s_pl
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
       u =  vmax * ( cos(GRD_s(g,k0,l,GRD_LAT)) * cos(alpha) &
                   + sin(GRD_s(g,k0,l,GRD_LAT))              &
                   * cos(GRD_s(g,k0,l,GRD_LON)) * sin(alpha) )
       v = -vmax * ( sin(GRD_s(g,k0,l,GRD_LON)) * sin(alpha) )

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

          u =  vmax * ( cos(GRD_s_pl(g,k0,l,GRD_LAT)) * cos(alpha) &
                      + sin(GRD_s_pl(g,k0,l,GRD_LAT))              &
                      * cos(GRD_s_pl(g,k0,l,GRD_LON)) * sin(alpha) )
          v = -vmax * ( sin(GRD_s_pl(g,k0,l,GRD_LON)) * sin(alpha) )

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
  subroutine GTL_clip_region( v, v_clip, kmin, kmax )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(RP), intent(in)  :: v     (ADM_gall   ,ADM_kall       ,ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_gall_in,1:(kmax-kmin+1),ADM_lall)

    integer :: i, j, k, l, n
    !---------------------------------------------------------------------------

    do l = 1,    ADM_lall
    do k = kmin, kmax
       n = 1
       do j = ADM_gmin, ADM_gmax+1
       do i = ADM_gmin, ADM_gmax+1
          v_clip(n,k-kmin+1,l) = v(suf(i,j),k,l)

          n = n + 1
       enddo
       enddo
    enddo
    enddo

    return
  end subroutine GTL_clip_region

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer( v, v_clip )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    real(RP), intent(in)  :: v     (ADM_gall   ,ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_gall_in,ADM_lall)

    integer :: i, j, l, n
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       n = 1
       do j = ADM_gmin, ADM_gmax+1
       do i = ADM_gmin, ADM_gmax+1
          v_clip(n,l) = v(suf(i,j),l)

          n = n + 1
       enddo
       enddo
    enddo

    return
  end subroutine GTL_clip_region_1layer

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer_k(v,v_clip,ksize,k)
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_gmin,    &
       ADM_gmax
    implicit none

    integer,  intent(in)  :: ksize
    real(RP), intent(in)  :: v     (ADM_gall   ,ksize, ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_gall_in,ADM_lall)
    integer,  intent(in)  :: k

    integer :: i, j, l, n
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       n = 1
       do j = ADM_gmin, ADM_gmax+1
       do i = ADM_gmin, ADM_gmax+1
          v_clip(n,l) = v(suf(i,j),k,l)

          n = n + 1
       enddo
       enddo
    enddo

    return
  end subroutine GTL_clip_region_1layer_k

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_gtl
