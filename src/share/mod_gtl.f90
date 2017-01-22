!-------------------------------------------------------------------------------
!> Module generic tool
!!
!! @par Description
!!         This module is for the generic subroutine
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
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
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
