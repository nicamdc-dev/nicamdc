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
