!-------------------------------------------------------------------------------
!>
!! variable conversion module
!!
!! @par Description
!!         Conversion tools for prognostic variables
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)   Imported from igdc-4.34
!! @li      2009-07-10 (H.Tomita)   Change the cnvvar_rhokin, cnvvar_kin for the energy conservation.
!! @li      2011-07-22 (T.Ohno)     add subroutines for plane hgrid systems
!!
!<
module mod_cnvvar
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: cnvvar_d2p
  public :: cnvvar_p2d
  public :: cnvvar_rhokin
  public :: cnvvar_kin

  public :: cnvvar_rhokin_ijkl

  public :: cnvvar_rhokin_plane ! [add] T.Ohno 110722
  public :: cnvvar_kin_plane    ! [add] T.Ohno 110722
  public :: cnvvar_ve_plane     ! [add] T.Ohno 110722

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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
  subroutine cnvvar_d2p(&
       ijdim,           & !--- IN : number of horizontal grid
       pre,             &  !--- IN : pressure ( physical )
       tem,             &  !--- IN : temperature ( physical )
       vx,              &  !--- IN : Vx ( physical )
       vy,              &  !--- IN : Vy ( physical )
       vz,              &  !--- IN : Vz ( physical )
       w,               &  !--- IN : w  ( physical )
       q,               &  !--- IN : mass concentration ( physical )
       rhog,            &  !--- OUT : density ( gam2 X G^{1/2} )
       rhogvx,          &  !--- OUT : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy,          &  !--- OUT : rho*Vy  ( gam2 X G^{1/2} )
       rhogvz,          &  !--- OUT : rho*Vz  ( gam2 X G^{1/2} )
       rhogw,           &  !--- OUT : rho*w   ( gam2 X G^{1/2} )
       rhoge,           &  !--- OUT : rho*ein ( gam2 X G^{1/2} )
       rhogq,           &  !--- OUT : rho*q   ( gam2 X G^{1/2} )
       gsqrtgam2,       &  !--- IN : G^{1/2} at the cell center
       gsqrtgam2h       &  !--- IN : G^{1/2} at the cell wall
       )
    !------
    !------ Generation of prognostic variables ( e.g. rhogvx )
    !------     1. calculation domain : (:,:,:)
    !------
    !
    use mod_adm, only :       &
         kdim => ADM_kall,    &
         kmin => ADM_kmin,    &
         kmax => ADM_kmax
    use mod_grd, only :  &
         GRD_afac,       &
         GRD_bfac
    use mod_thrmdyn, only :    &
         thrmdyn_ein,      &
         thrmdyn_rho,      &
         thrmdyn_qd
    use mod_runconf, only : &
         nqmax => TRC_VMAX
    !
    implicit none
    !
    integer, intent(in) :: ijdim
    real(8), intent(in) :: pre(1:ijdim,1:kdim)
    real(8), intent(in) :: tem(1:ijdim,1:kdim)
    real(8), intent(in) :: vx(1:ijdim,1:kdim)
    real(8), intent(in) :: vy(1:ijdim,1:kdim)
    real(8), intent(in) :: vz(1:ijdim,1:kdim)
    real(8), intent(in) :: w(1:ijdim,1:kdim)
    real(8), intent(in) :: q(1:ijdim,1:kdim,1:nqmax)
    real(8), intent(in) :: gsqrtgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2h(1:ijdim,1:kdim)
    !
    real(8), intent(out) :: rhog(1:ijdim,1:kdim)
    real(8), intent(out) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(out) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(out) :: rhogvz(1:ijdim,1:kdim)
    real(8), intent(out) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(out) :: rhoge(1:ijdim,1:kdim)
    real(8), intent(out) :: rhogq(1:ijdim,1:kdim,1:nqmax)
    !
    real(8) :: grho_h(1:ijdim,1:kdim)
    real(8) :: ein(1:ijdim,1:kdim)
    real(8) :: qd(1:ijdim,1:kdim)
    real(8) :: rho(1:ijdim,1:kdim)
    !
    integer :: nq
    integer :: k
    !
    !--- calculation of dry mass concentration
    call thrmdyn_qd( &
         ijdim,      & !--- in
         qd,         & !--- out
         q )           !--- in
    !
    !--- calculation  of density 
    call thrmdyn_rho( &
         ijdim,       & !--- in
         rho,         & !--- OUT
         pre,         & !--- IN
         tem,         & !--- IN
         qd,          & !--- IN
         q  )           !--- IN
    !
    !--- calculation of internal energy
    call thrmdyn_ein( &
         ijdim,       & !--- in
         ein,         & !--- OUT
         tem,         & !--- IN
         qd,          & !--- IN
         q )            !--- IN
    !
    !--- calculation of rhoge, rhogvx, rhogvy, rhogvz
    rhog(:,:) = rho(:,:) * gsqrtgam2(:,:)
    rhoge(:,:) = ein(:,:) * rhog(:,:)
    rhogvx(:,:) = vx(:,:) * rhog(:,:)
    rhogvy(:,:) = vy(:,:) * rhog(:,:)
    rhogvz(:,:) = vz(:,:) * rhog(:,:)
    do nq = 1,nqmax
       rhogq(:,:,nq) = q(:,:,nq) * rhog(:,:)
    end do
    !
    !
    !--- calculation of rhogw
    !------ interpolation of grho_h
    do k =kmin,kmax+1
       grho_h(:,k) = 0.5D0 *           &
            ( GRD_afac(k) * rho(:,k  ) &
            + GRD_bfac(k) * rho(:,k-1) &
            ) * gsqrtgam2h(:,k)
    end do
    grho_h(:,kmin-1) = grho_h(:,kmin)
    !
    !------ grho_h is ( gam2 * G^{1/2} * rho ) at the cell wall.
    rhogw(:,:) = w(:,:)* grho_h(:,:)
    !
  end subroutine cnvvar_d2p
  !-----------------------------------------------------------------------------
  subroutine cnvvar_p2d(&
       ijdim,           &  !--- IN : number of horizontal grid
       rho,             &  !--- OUT : density ( physical )
       pre,             &  !--- OUT : pressure ( physical )
       tem,             &  !--- OUT : temperature ( physical )
       vx,              &  !--- OUT : Vx ( physical )
       vy,              &  !--- OUT : Vy ( physical )
       vz,              &  !--- OUT : Vz ( physical )
       w,               &  !--- OUT : w  ( physical )
       q,               &  !--- OUT : mass concentration ( physical )
       rhog,            &  !--- IN : density ( gam2 X G^{1/2} )
       rhogvx,          &  !--- IN : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy,          &  !--- IN : rho*Vy  ( gam2 X G^{1/2} )
       rhogvz,          &  !--- IN : rho*Vz  ( gam2 X G^{1/2} )
       rhogw,           &  !--- IN : rho*w   ( gam2 X G^{1/2} )
       rhoge,           &  !--- IN : rho*ein ( gam2 X G^{1/2} )
       rhogq,           &  !--- IN : rho*q   ( gam2 X G^{1/2} )
       gsqrtgam2,       &  !--- IN : G^{1/2} at the cell center
       gsqrtgam2h       &  !--- IN : G^{1/2} at the cell wall
       )
    !------
    !------ Generation of diagnostic variables ( e.g. vx )
    !------     1. calculation domain : (:,:)
    !------
    !
    use mod_adm, only :       &
         kdim => ADM_kall,    &
         kmin => ADM_kmin,    &
         kmax => ADM_kmax
    use mod_grd, only :  &
         GRD_afac,       &
         GRD_bfac
    use mod_thrmdyn, only :    &
         thrmdyn_tempre,   &
         thrmdyn_qd
    use mod_runconf, only : &
         nqmax => TRC_VMAX
    !
    implicit none
    !
    integer, intent(in) :: ijdim
    real(8), intent(out) :: rho(1:ijdim,1:kdim)
    real(8), intent(out) :: pre(1:ijdim,1:kdim)
    real(8), intent(out) :: tem(1:ijdim,1:kdim)
    real(8), intent(out) :: vx(1:ijdim,1:kdim)
    real(8), intent(out) :: vy(1:ijdim,1:kdim)
    real(8), intent(out) :: vz(1:ijdim,1:kdim)
    real(8), intent(out) :: w(1:ijdim,1:kdim)
    real(8), intent(out) :: q(1:ijdim,1:kdim,1:nqmax)
    !
    real(8), intent(in) :: gsqrtgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2h(1:ijdim,1:kdim)
    !
    real(8), intent(in) :: rhog(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvz(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(in) :: rhoge(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogq(1:ijdim,1:kdim,1:nqmax)
    !
    real(8) :: rrhog(1:ijdim,1:kdim)
!    real(8) :: grho_h(1:ijdim,1:kdim)
    real(8) :: ein(1:ijdim,1:kdim)
    real(8) :: qd(1:ijdim,1:kdim)
    !
    integer :: nq
    integer :: k,ij
    !
    !--- calculation of rgrho
!del
    do k = 1, kdim
!del
!del
       do ij = 1, ijdim
          rrhog(ij,k) = 1.0D0 / rhog(ij,k)
       end do
    end do
    !
    !--- calculation of ein vx,vy,vz
!del
    do k = 1, kdim
!del
!del
       do ij = 1, ijdim
          rho(ij,k) = rhog(ij,k) / gsqrtgam2(ij,k)
       end do
!del
!del
       do ij = 1, ijdim
          ein(ij,k) = rhoge(ij,k) * rrhog(ij,k)
       end do
!del
!del
       do ij = 1, ijdim
          vx(ij,k) = rhogvx(ij,k) * rrhog(ij,k)
       end do
!del
!del
       do ij = 1, ijdim
          vy(ij,k) = rhogvy(ij,k) * rrhog(ij,k)
       end do
!del
!del
       do ij = 1, ijdim
          vz(ij,k) = rhogvz(ij,k) * rrhog(ij,k)
       end do
    end do
!del
    do nq = 1,nqmax
!del
       do k = 1, kdim
!del
          do ij = 1, ijdim
             q(ij,k,nq) = rhogq(ij,k,nq) * rrhog(ij,k)
          end do
       end do
    end do
    !
    !--- calculation of dry mass concentration
    call thrmdyn_qd( &
         ijdim,      &
         qd,         & !--- OUT  : dry mass concentration
         q )           !--- IN  : mass concentration
    !
    !
    !--- calculation of tem, pre, and th
    call thrmdyn_tempre( &
         ijdim,          &
         tem,            &  !--- out
         pre,            &  !--- out
         ein,            &  !--- in
         rho,            &  !--- in
         qd,             &  !--- in
         q )                !--- in
    !
    !---- Vertical velocity ( w )
    !------- calculation of grho at the call wall
!del
    do k = kmin, kmax+1
!del
       do ij = 1, ijdim
          w(ij,k) = rhogw(ij,k) / ( 0.5D0 * &
               ( GRD_afac(k) * rho(ij,k  )  &
               + GRD_bfac(k) * rho(ij,k-1)  &
               ) * gsqrtgam2h(ij,k) )
       end do
    end do
    !
!del
!del
    do ij = 1, ijdim
       w(ij,kmin-1) = rhogw(ij,kmin-1) / ( 0.5D0 * &
            ( GRD_afac(kmin) * rho(ij,kmin  )  &
            + GRD_bfac(kmin) * rho(ij,kmin-1)  &
            ) * gsqrtgam2h(ij,kmin) )
    end do
    !
  end subroutine cnvvar_p2d
  !-----------------------------------------------------------------------------
  subroutine cnvvar_rhokin(&
       ijdim,              & !--- IN : number of horizontal grid
       rhog,               &  !--- IN : rho     ( gam2 X G^{1/2}  )
       rhogvx,             &  !--- IN : rho*Vx  ( gam2 X G^{1/2}  )
       rhogvy,             &  !--- IN : rho*Vy  ( gam2 X G^{1/2}  )
       rhogvz,             &  !--- IN : rho*Vz  ( gam2 X G^{1/2}  )
       rhogw,              &  !--- IN : rho*w   ( gam2 X G^{1/2}  )
       gsqrtgam2,          &  !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,         &  !--- IN : G^{1/2} at the cell wall
       rhogkin,            &  !--- OUT : 1/2 rho*v^2  ( gam2 X G^{1/2}  )
       rhogkin_h,          &
       rhogkin_v           &
       )      
    !------ 
    !------ Calculation of kinetic energy.
    !------ 
    !
    use mod_adm, only :       &
         kdim => ADM_kall,    &
         kmin => ADM_kmin,    &
         kmax => ADM_kmax
    use mod_grd, only :      &
         GRD_afac, GRD_bfac,         GRD_cfac, GRD_dfac
    !
    implicit none
    integer, intent(in) :: ijdim
    real(8), intent(in) :: rhog(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvz(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2h(1:ijdim,1:kdim)
    real(8), intent(out) :: rhogkin(1:ijdim,1:kdim)
    !
    real(8) :: rhog_h(1:ijdim,1:kdim)
    real(8), intent(out), optional :: rhogkin_h(1:ijdim,1:kdim)
    real(8), intent(out), optional :: rhogkin_v(1:ijdim,1:kdim)
    real(8) :: rhogkin_h0(1:ijdim,1:kdim)
    real(8) :: rhogkin_v0(1:ijdim,1:kdim)
    !
    integer :: k,ij





    !--- rhogkin = gamma^2 * g_sqrt * rho * kin
    !
    do ij = 1, ijdim
       rhogkin_v0(ij,kmin) = 0.0D0
    end do
!del
!del
!del
    do ij = 1, ijdim
       rhogkin_v0(ij,kmax+1) = 0.0D0
    end do

    !--- horizontal kinetic energy
!del
    do k = kmin, kmax
!del
       do ij = 1, ijdim
          rhogkin_h0(ij,k)  &
               =((rhogvx(ij,k)/rhog(ij,k))**2 &
               + (rhogvy(ij,k)/rhog(ij,k))**2 &
               + (rhogvz(ij,k)/rhog(ij,k))**2)&
               * rhog(ij,k)*0.5D0
       end do
    end do
    !
    !--- rhog at the half level
!del
    do k = kmin, kmax+1
!del
       do ij = 1, ijdim
          rhog_h(ij,k) = 0.5D0 *            &
               ( GRD_afac(k) * rhog(ij,k  ) / gsqrtgam2(ij,k  ) &
               + GRD_bfac(k) * rhog(ij,k-1) / gsqrtgam2(ij,k-1) &
               ) * gsqrtgam2h(ij,k)
       end do
    end do
    !
    !--- vertical kinetic energy
!del
    do k = kmin+1, kmax
!del
       do ij = 1, ijdim
          rhogkin_v0(ij,k) &
               =((rhogw(ij,k)/rhog_h(ij,k))**2) &
               * rhog_h(ij,k)*0.5D0
       end do
    end do
    !
!del
    do k = kmin, kmax
!del
       do ij = 1, ijdim
          rhogkin(ij,k) = rhogkin_h0(ij,k)&
               + (GRD_dfac(k) * rhogkin_v0(ij,k+1)&
               +  GRD_cfac(k) * rhogkin_v0(ij,k  ) )*0.5D0
       end do
    end do
    !
    if(present(rhogkin_h)) then
!del
       do k = 1, kdim
!del
          do ij = 1, ijdim
             rhogkin_h(ij,k) = rhogkin_h0(ij,k)
          end do
       end do
    end if
    if(present(rhogkin_v)) then
!del
       do k = 1, kdim
!del
          do ij = 1, ijdim
             rhogkin_v(ij,k) = rhogkin_v0(ij,k)
          end do
       end do
    end if





  end subroutine cnvvar_rhokin
  !-----------------------------------------------------------------------------
  subroutine cnvvar_kin(&
       ijdim,           &  !--- IN : number of horizontal grid
       rhog,            &  !--- IN : density ( gam2 X G^{1/2} )
       rhogvx,          &  !--- IN : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy,          &  !--- IN : rho*Vy  ( gam2 X G^{1/2} )
       rhogvz,          &  !--- IN : rho*Vz  ( gam2 X G^{1/2} )
       rhogw,           &  !--- IN : rho*w   ( gam2 X G^{1/2} )
       gsqrtgam2,       &  !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,      &  !--- IN : G^{1/2} at the cell wall
       kin              &  !--- OUT : 1/2 v^2 ( physical )
       )
    !------ 
    !------ Calculation of kinetic energy.
    !------ 
    !
    use mod_adm, only :       &
         kdim => ADM_kall
    implicit none

    integer, intent(in) :: ijdim
    real(8), intent(in) :: rhog(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvz(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2h(1:ijdim,1:kdim)
    real(8), intent(out) :: kin(1:ijdim,1:kdim)

    call cnvvar_rhokin(&
         ijdim,        &
         rhog,         &
         rhogvx,       &
         rhogvy,       &
         rhogvz,       &
         rhogw,        &
         gsqrtgam2,    &
         gsqrtgam2h,   &
         kin           )      

    kin(:,:) = kin(:,:) / rhog(:,:)

  end subroutine cnvvar_kin

  !-----------------------------------------------------------------------------
  subroutine cnvvar_rhokin_plane(&
       ijdim,              & !--- IN : number of horizontal grid
       rhog,               &  !--- IN : rho     ( gam2 X G^{1/2}  )
       rhogvx,             &  !--- IN : rho*Vx  ( gam2 X G^{1/2}  )
       rhogvy,             &  !--- IN : rho*Vy  ( gam2 X G^{1/2}  )
       rhogw,              &  !--- IN : rho*w   ( gam2 X G^{1/2}  )
       gsqrtgam2,          &  !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,         &  !--- IN : G^{1/2} at the cell wall
       rhogkin,            &  !--- OUT : 1/2 rho*v^2  ( gam2 X G^{1/2}  )
       rhogkin_h,          &
       rhogkin_v           &
       )      
    ! [add] T.Ohno 110722
    !------ 
    !------ Calculation of kinetic energy.
    !------ 
    !
    use mod_adm, only :       &
         kdim => ADM_kall,    &
         kmin => ADM_kmin,    &
         kmax => ADM_kmax,    &
         ADM_VMISS
    use mod_grd, only :      &
         GRD_afac, GRD_bfac,         GRD_cfac, GRD_dfac
    !
    implicit none
    integer, intent(in) :: ijdim
    real(8), intent(in) :: rhog(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2h(1:ijdim,1:kdim)
    real(8), intent(out) :: rhogkin(1:ijdim,1:kdim)
    !
    real(8) :: rhog_h(1:ijdim,1:kdim)
    real(8), intent(out), optional :: rhogkin_h(1:ijdim,1:kdim)
    real(8), intent(out), optional :: rhogkin_v(1:ijdim,1:kdim)
    real(8) :: rhogkin_h0(1:ijdim,1:kdim)
    real(8) :: rhogkin_v0(1:ijdim,1:kdim)
    !
    integer :: k
    !
    !--- rhogkin = gamma^2 * g_sqrt * rho * kin
    !
    !--- horizontal kinetic energy
    rhogkin_h0(:,kmin-1) = ADM_VMISS
    rhogkin_h0(:,kmax+1) = ADM_VMISS
    do k = kmin, kmax
       rhogkin_h0(:,k)  &
            =((rhogvx(:,k)/rhog(:,k))**2 &
            + (rhogvy(:,k)/rhog(:,k))**2)&
            * rhog(:,k)*0.5D0
    end do
    !
    !--- rhog at the half level
    do k =kmin,kmax+1
       rhog_h(:,k) = 0.5D0 *            &
            ( GRD_afac(k) * rhog(:,k  )/gsqrtgam2(:,k  ) &
            + GRD_bfac(k) * rhog(:,k-1)/gsqrtgam2(:,k-1) &
            ) * gsqrtgam2h(:,k)
    end do
    !
    !--- vertical kinetic energy
    rhogkin_v0(:,kmin-1) = ADM_VMISS
    rhogkin_v0(:,kmin) = 0.0D0
    rhogkin_v0(:,kmax+1) = 0.0D0
    do k = kmin+1, kmax
       rhogkin_v0(:,k) &
            =((rhogw(:,k)/rhog_h(:,k))**2) &
            * rhog_h(:,k)*0.5D0
    end do
    !
    rhogkin(:,kmin-1) = ADM_VMISS
    rhogkin(:,kmax+1) = ADM_VMISS
    do k = kmin, kmax
       rhogkin(:,k) = rhogkin_h0(:,k)&
            + (GRD_dfac(k  )*rhogkin_v0(:,k+1)&
            +  GRD_cfac(k  )*rhogkin_v0(:,k  ) )*0.5D0
    end do
    !
    if(present(rhogkin_h)) then
       rhogkin_h = rhogkin_h0
    end if
    if(present(rhogkin_v)) then
       rhogkin_v = rhogkin_v0
    end if
    !
  end subroutine cnvvar_rhokin_plane

  !-----------------------------------------------------------------------------
  ! [add] T.Ohno 110722
  subroutine cnvvar_kin_plane(&
       ijdim,           &  !--- IN : number of horizontal grid
       rhog,            &  !--- IN : density ( gam2 X G^{1/2} )
       rhogvx,          &  !--- IN : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy,          &  !--- IN : rho*Vy  ( gam2 X G^{1/2} )
       rhogw,           &  !--- IN : rho*w   ( gam2 X G^{1/2} )
       gsqrtgam2,       &  !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,      &  !--- IN : G^{1/2} at the cell wall
       kin              &  !--- OUT : 1/2 v^2 ( physical )
       )
    use mod_adm, only :       &
         kdim => ADM_kall
    implicit none

    integer, intent(in) :: ijdim
    real(8), intent(in) :: rhog(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2h(1:ijdim,1:kdim)
    real(8), intent(out) :: kin(1:ijdim,1:kdim)

    call cnvvar_rhokin_plane(&
         ijdim,        &
         rhog,         &
         rhogvx,       &
         rhogvy,       &
         rhogw,        &
         gsqrtgam2,    &
         gsqrtgam2h,   &
         kin           &
       )      
    kin(:,:) = kin(:,:) / rhog(:,:)
    !
  end subroutine cnvvar_kin_plane
  !-----------------------------------------------------------------------------
  subroutine cnvvar_ve_plane( & 
       ijdim,           & !--- IN : number of horizontal grid
       rhog,            & !--- INOUT : density ( gam2 X G^{1/2} )
       rhogvx,          & !--- INOUT : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy,          & !--- INOUT : rho*Vy  ( gam2 X G^{1/2} )
       rhogw,           & !--- INOUT : rho*w   ( gam2 X G^{1/2} )
       rhoge,           & !--- INOUT : rho*ein ( gam2 X G^{1/2} )
       rhogq,           & !--- INOUT : rho*q   ( gam2 X G^{1/2} )
       rho,             & !--- OUT : density( physical )
       vx,              & !--- OUT : Vx ( physical )
       vy,              & !--- OUT : Vy ( physical )
       w,               & !--- OUT : w  ( physical )
       ein,             & !--- OUT : int. enrg. ( physical )
       q,               & !--- OUT : mass concentration. ( physical )
       pre,             & !--- OUT : pressure   ( physical )
       tem,             & !--- OUT : temperature( physical )
       rhod,            & !--- OUT : perturb dens. ( gam2 X G^{1/2} )
       pred,            & !--- OUT : perturb pres. ( gam2 X G^{1/2} )
       temd,            & !--- OUT : perturb temp. ( physical )
       pre_bs,          & !--- IN : pressure basicstate
       tem_bs,          & !--- IN : temperature basicstate
       rho_bs,          & !--- IN : density basicstate
       phi,             & !--- IN : geoptential
       gsqrtgam2,       & !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,      & !--- IN : G^{1/2} at the cell wall
       g3xh,            & !--- IN : G3X at the cell wall
       g3yh,            & !--- IN : G3Y at the cell wall
       bc_flag)
    ! [add] T.Ohno 110722
    !------
    !------ Generation of the phyisical variables from prognositc variables.
    !------    1. with boundary condision setting for all variables.
    !------    2. perturbation quantities are also calculated.
    !------
    !
    use mod_adm, only :       &
         kdim => ADM_kall,    &
         kmin => ADM_kmin,    &
         kmax => ADM_kmax
    use mod_grd, only :  &
         GRD_afac,       &
         GRD_bfac
    use mod_bndcnd, only : &
         bndcnd_all_plane
    use mod_thrmdyn, only :    &
         thrmdyn_tempre,   &
         thrmdyn_qd
    use mod_runconf, only : &
         nqmax => TRC_VMAX
    !
    implicit none
    !
    integer, intent(in) :: ijdim
    real(8), intent(inout) :: rhog(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhoge(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhogq(1:ijdim,1:kdim,1:nqmax)
    !
    real(8), intent(out) :: rho(1:ijdim,1:kdim)
    real(8), intent(out) :: pre(1:ijdim,1:kdim)
    real(8), intent(out) :: tem(1:ijdim,1:kdim)
    real(8), intent(out) :: rhod(1:ijdim,1:kdim)
    real(8), intent(out) :: pred(1:ijdim,1:kdim)
    real(8), intent(out) :: temd(1:ijdim,1:kdim)
    real(8), intent(out) :: ein(1:ijdim,1:kdim)
    real(8), intent(out) :: q(1:ijdim,1:kdim,1:nqmax)
    real(8), intent(out) :: vx(1:ijdim,1:kdim)
    real(8), intent(out) :: vy(1:ijdim,1:kdim)
    real(8), intent(out) :: w(1:ijdim,1:kdim)
    !
    real(8), intent(in)  :: phi(1:ijdim,1:kdim)
    real(8), intent(in)  :: pre_bs(1:ijdim,1:kdim)
    real(8), intent(in)  :: tem_bs(1:ijdim,1:kdim)
    real(8), intent(in)  :: rho_bs(1:ijdim,1:kdim)
    !
    real(8), intent(in)  :: gsqrtgam2(1:ijdim,1:kdim)
    real(8), intent(in)  :: gsqrtgam2h(1:ijdim,1:kdim)
    real(8), intent(in)  :: g3xh(1:ijdim,1:kdim)
    real(8), intent(in)  :: g3yh(1:ijdim,1:kdim)
    !
    logical, intent(in) ::  bc_flag
    !
    real(8) :: rhog_h(1:ijdim,1:kdim)
    real(8) :: qd(1:ijdim,1:kdim)
    real(8) :: rrhog(1:ijdim,1:kdim)
    !
    integer :: nq
    integer :: k
    !
    !--- calculation of rrhog
    rrhog(:,:) = 1.0D0/rhog(:,:)
    !
    !--- calculation of rho, ein vx,vy,vz, q
    rho(:,:) = rhog(:,:) / gsqrtgam2(:,:)
    ein(:,:) = rhoge(:,:) * rrhog(:,:)
    vx(:,:) = rhogvx(:,:) * rrhog(:,:)
    vy(:,:) = rhogvy(:,:) * rrhog(:,:)
    do nq=1,nqmax
       q(:,:,nq) = rhogq(:,:,nq) * rrhog(:,:)
    end do
    !
    !--- calculation of dry mass concentration
    call thrmdyn_qd( &
         ijdim,      & !--- in
         qd,         & !--- out
         q )           !--- in
    !
    !--- calculation of tem, pre, and th
    call thrmdyn_tempre( &
         ijdim,          &
         tem,            &  !--- out
         pre,            &  !--- out
         ein,            &  !--- in
         rho,            &  !--- in
         qd,             &  !--- in
         q )                !--- in
    !
    !---- Vertical velocity ( w )
    do k =kmin,kmax+1
       rhog_h(:,k) = 0.5D0 *           &
            ( GRD_afac(k) * rho(:,k  ) &
            + GRD_bfac(k) * rho(:,k-1) &
            ) * gsqrtgam2h(:,k)
    end do
    rhog_h(:,kmin-1) = rhog_h(:,kmin)
    !
    !------ except boundaries of [kmin] and [kmax+1]
    w(:,kmin+1:kmax)            &
         = rhogw(:,kmin+1:kmax) &
         / rhog_h(:,kmin+1:kmax)
    !
    !--- Setting the boundary conditions
    if(bc_flag) then
       call bndcnd_all_plane(&
            ijdim,     &  !--- in
            vx,        &  !--- inout
            vy,        &  !--- inout
            w,         &  !--- inout
            tem,       &  !--- inout
            rho,       &  !--- inout
            pre,       &  !--- inout
            ein,       &  !--- inout
            rhog,      &  !--- inout
            rhogvx,    &  !--- inout
            rhogvy,    &  !--- inout
            rhogw,     &  !--- inout
            rhoge,     &  !--- inout
            phi,       &  !--- in
            gsqrtgam2, &  !--- in
            gsqrtgam2h,&  !--- in
            g3xh,      &  !--- in
            g3yh       &  !--- in
            )
    end if
    !
    !--- Calculation of perturbations ( pred, rhod, temd )
    pred(:,:) = ( pre(:,:) - pre_bs(:,:) ) * gsqrtgam2(:,:)
    rhod(:,:) = ( rho(:,:) - rho_bs(:,:) ) * gsqrtgam2(:,:)
    temd(:,:) = tem(:,:) - tem_bs(:,:)
    !
  end subroutine cnvvar_ve_plane

  !-----------------------------------------------------------------------------
  subroutine cnvvar_rhokin_ijkl( &
       rhog,    rhog_pl,   &
       rhogvx,  rhogvx_pl, &
       rhogvy,  rhogvy_pl, &
       rhogvz,  rhogvz_pl, &
       rhogw,   rhogw_pl,  &
       rhogkin, rhogkin_pl )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac, &
       GRD_cfac, &
       GRD_dfac
    use mod_vmtr, only: &
       VMTR_GSGAM2H,     &
       VMTR_GSGAM2H_pl,  &
       VMTR_RGSGAM2,     &
       VMTR_RGSGAM2_pl
    implicit none

    real(8), intent(in)  :: rhog      (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 )
    real(8), intent(in)  :: rhog_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X vx
    real(8), intent(in)  :: rhogvx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X vy
    real(8), intent(in)  :: rhogvy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X vz
    real(8), intent(in)  :: rhogvz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw     (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X w
    real(8), intent(in)  :: rhogw_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogkin   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X kin
    real(8), intent(out) :: rhogkin_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: rhogkin_h   (ADM_gall,   ADM_kall) ! rho X ( G^{1/2} X gamma2 ) X kin (horizontal)
    real(8) :: rhogkin_h_pl(ADM_gall_pl,ADM_kall)
    real(8) :: rhogkin_v   (ADM_gall,   ADM_kall) ! rho X ( G^{1/2} X gamma2 ) X kin (vertical)
    real(8) :: rhogkin_v_pl(ADM_gall_pl,ADM_kall)

    real(8) :: rhog_h   (ADM_gall,   ADM_kall)
    real(8) :: rhog_h_pl(ADM_gall_pl,ADM_kall)

    integer :: ij, k, l
    !---------------------------------------------------------------------------

!del
    do l = 1, ADM_lall
!del
       do k = ADM_kmin+1, ADM_kmax
          !--- rhog at the half level
!del
          do ij = 1, ADM_gall
             rhog_h(ij,k) = 0.5D0 * ( GRD_afac(k) * VMTR_RGSGAM2(ij,k  ,l) * rhog(ij,k  ,l) &
                                    + GRD_bfac(k) * VMTR_RGSGAM2(ij,k-1,l) * rhog(ij,k-1,l) &
                                    ) * VMTR_GSGAM2H(ij,k,l)
          enddo
       enddo
!del
       do k = ADM_kmin+1, ADM_kmax
          !--- vertical kinetic energy
          do ij = 1, ADM_gall
             rhogkin_v(ij,k) = 0.5D0 * ( rhogw(ij,k,l) * rhogw(ij,k,l) ) / rhog_h(ij,k)
          enddo
       enddo
!del
       do ij = 1, ADM_gall
          rhogkin_v(ij,ADM_kmin  ) = 0.D0
          rhogkin_v(ij,ADM_kmax+1) = 0.D0
       enddo

!del
       do k = ADM_kmin, ADM_kmax
          !--- horizontal kinetic energy
          do ij = 1, ADM_gall
             rhogkin_h(ij,k) = 0.5D0 * ( rhogvx(ij,k,l) * rhogvx(ij,k,l) &
                                       + rhogvy(ij,k,l) * rhogvy(ij,k,l) &
                                       + rhogvz(ij,k,l) * rhogvz(ij,k,l) ) / rhog(ij,k,l)
          enddo
       enddo

!del
       do k = ADM_kmin, ADM_kmax
          !--- total kinetic energy
!del
          do ij = 1, ADM_gall
             rhogkin(ij,k,l) = rhogkin_h(ij,k)                           &
                             + 0.5D0 * ( GRD_dfac(k) * rhogkin_v(ij,k+1) &
                                       + GRD_cfac(k) * rhogkin_v(ij,k  ) )
          enddo
       enddo
!del
       do ij = 1, ADM_gall
          rhogkin(ij,ADM_kmin-1,l) = 0.D0
          rhogkin(ij,ADM_kmax+1,l) = 0.D0
       enddo

    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
!del
       do l = 1, ADM_lall_pl
!del
          do k = ADM_kmin+1, ADM_kmax
             !--- rhog at the half level
!del
             do ij = 1, ADM_gall_pl
                rhog_h_pl(ij,k) = 0.5D0 * ( GRD_afac(k) * VMTR_RGSGAM2_pl(ij,k  ,l) * rhog_pl(ij,k  ,l) &
                                          + GRD_bfac(k) * VMTR_RGSGAM2_pl(ij,k-1,l) * rhog_pl(ij,k-1,l) &
                                          ) * VMTR_GSGAM2H_pl(ij,k,l)
             enddo
          enddo
!del
          do k = ADM_kmin+1, ADM_kmax
             !--- vertical kinetic energy
             do ij = 1, ADM_gall_pl
                rhogkin_v_pl(ij,k) = 0.5D0 * ( rhogw_pl(ij,k,l) * rhogw_pl(ij,k,l) ) / rhog_h_pl(ij,k)
             enddo
          enddo
!del
          do ij = 1, ADM_gall_pl
             rhogkin_v_pl(ij,ADM_kmin  ) = 0.D0
             rhogkin_v_pl(ij,ADM_kmax+1) = 0.D0
          enddo
!del
          do k = ADM_kmin, ADM_kmax
             !--- horizontal kinetic energy
             do ij = 1, ADM_gall_pl
                rhogkin_h_pl(ij,k) = 0.5D0 * ( rhogvx_pl(ij,k,l) * rhogvx_pl(ij,k,l) &
                                             + rhogvy_pl(ij,k,l) * rhogvy_pl(ij,k,l) &
                                             + rhogvz_pl(ij,k,l) * rhogvz_pl(ij,k,l) ) / rhog_pl(ij,k,l)
             enddo
          enddo

!del
          do k = ADM_kmin, ADM_kmax
             !--- total kinetic energy
             do ij = 1, ADM_gall_pl
                rhogkin_pl(ij,k,l) = rhogkin_h_pl(ij,k)                           &
                                   + 0.5D0 * ( GRD_dfac(k) * rhogkin_v_pl(ij,k+1) &
                                             + GRD_cfac(k) * rhogkin_v_pl(ij,k  ) )
             enddo
          enddo
          do ij = 1, ADM_gall_pl
             rhogkin_pl(ij,ADM_kmin-1,l) = 0.D0
             rhogkin_pl(ij,ADM_kmax+1,l) = 0.D0
          enddo

       enddo
    endif

    return
  end subroutine cnvvar_rhokin_ijkl

end module mod_cnvvar
!-------------------------------------------------------------------------------
