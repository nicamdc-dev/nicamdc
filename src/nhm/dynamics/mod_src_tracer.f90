!-------------------------------------------------------------------------------
!>
!! Tracer advection module
!!
!! @par Description
!!         This module contains subroutines for tracer advection
!!
!! @author  H.Tomita, Y.Niwa
!!
!! @par History
!! @li      2008-01-24 (Y.Niwa   ) Imported from mod_src and mod_oprt
!! @li      2013-11-08 (H.Yashiro) Re-Arrange
!!
!<
module mod_src_tracer
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID,    &
     TI  => ADM_TI,  &
     TJ  => ADM_TJ,  &
     AI  => ADM_AI,  &
     AIJ => ADM_AIJ, &
     AJ  => ADM_AJ,  &
     K0  => ADM_KNONE
  use mod_grd, only: &
     XDIR => GRD_XDIR, &
     YDIR => GRD_YDIR, &
     ZDIR => GRD_ZDIR
  use mod_gmtr, only: &
     P_RAREA => GMTR_P_RAREA, &
     T_RAREA => GMTR_T_RAREA, &
     W1      => GMTR_T_W1,    &
     W2      => GMTR_T_W2,    &
     W3      => GMTR_T_W3,    &
     HNX     => GMTR_A_HNX,   &
     HNY     => GMTR_A_HNY,   &
     HNZ     => GMTR_A_HNZ,   &
     HTX     => GMTR_A_HTX,   &
     HTY     => GMTR_A_HTY,   &
     HTZ     => GMTR_A_HTZ,   &
     TNX     => GMTR_A_TNX,   &
     TNY     => GMTR_A_TNY,   &
     TNZ     => GMTR_A_TNZ,   &
     TN2X    => GMTR_A_TN2X,  &
     TN2Y    => GMTR_A_TN2Y,  &
     TN2Z    => GMTR_A_TN2Z
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: src_tracer_advection
  public :: vertical_limiter_thuburn
  public :: horizontal_flux
  public :: horizontal_limiter_thuburn

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
  real(8), private, allocatable, save :: local_t_var   (:,:,:,:,:)
  real(8), private, allocatable, save :: local_t_var_pl(:,:,:,:)

  real(8), private, allocatable, save :: GRD_xr   (:,:,:,:,:)
  real(8), private, allocatable, save :: GRD_xr_pl(:,:,:,:)

  !-----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------------
  ! Y.Niwa add 080124
  ! This routine is revised version of src_tracer_advection
  subroutine src_tracer_advection( &
       vmax,                        & !--- IN    : number of tracers
       rhogq,       rhogq_pl,       & !--- INOUT : rhogq   ( gam2 X G^{1/2} )
       rhog_in,     rhog_in_pl,     & !--- IN    : rho(old)( gam2 X G^{1/2} )
       rhog_mean,   rhog_mean_pl,   & !--- IN    : rho     ( gam2 X G^{1/2} )
       rhogvx_mean, rhogvx_mean_pl, & !--- IN    : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy_mean, rhogvy_mean_pl, & !--- IN    : rho*Vy  ( gam2 X G^{1/2} )
       rhogvz_mean, rhogvz_mean_pl, & !--- IN    : rho*Vz  ( gam2 X G^{1/2} )
       rhogw_mean,  rhogw_mean_pl,  & !--- IN    : rho*w   ( gam2 X G^{1/2} )
       frhog,       frhog_pl,       & !--- IN    : hyperviscosity tendency for rhog
       dt,                          & !--- IN    : delta t
       thubern_lim                  ) !--- IN    : switch of thubern limiter
    use mod_adm, only :  &
       ADM_have_pl,    &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl,    &
       ADM_kmin,       &
       ADM_kmax
    use mod_cnst, only: &
       CNST_EPS_ZERO
    use mod_grd, only: &
       GRD_rdgz, &
       GRD_afac, &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_C2Wfact_Gz,    &
       VMTR_C2Wfact_Gz_pl, &
       VMTR_RGSH,          &
       VMTR_RGSH_pl,       &
       VMTR_RGAM,          &
       VMTR_RGAM_pl,       &
       VMTR_RGAMH,         &
       VMTR_RGAMH_pl
    implicit none

    integer, intent(in)    :: vmax
    real(8), intent(inout) :: rhogq         (ADM_gall,   ADM_kall,ADM_lall,   vmax)
    real(8), intent(inout) :: rhogq_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax)
    real(8), intent(in)    :: rhog_in       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhog_in_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhog_mean     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhog_mean_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogvx_mean   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogvx_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogvy_mean   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogvy_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogvz_mean   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogvz_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogw_mean    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogw_mean_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: frhog         (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: frhog_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: dt
    logical, intent(in)    :: thubern_lim  ![add] 20130613 R.Yoshida

    real(8) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rrhog    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rrhog_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: q        (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: d        (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: d_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: q_h      (ADM_gall,   ADM_kall,ADM_lall   )   ! q at layer face
    real(8) :: q_h_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: flx_v    (ADM_gall,   ADM_kall,ADM_lall   )   ! mass flux
    real(8) :: flx_v_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ck       (ADM_gall,   ADM_kall,ADM_lall,   2) ! Courant number
    real(8) :: ck_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(8) :: q_a      (6,ADM_gall,   ADM_kall,ADM_lall   ) ! q at cell face
    real(8) :: q_a_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: flx_h    (6,ADM_gall,   ADM_kall,ADM_lall   ) ! mass flux
    real(8) :: flx_h_pl (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ch       (6,ADM_gall,   ADM_kall,ADM_lall   ) ! Courant number
    real(8) :: ch_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: cmask    (6,ADM_gall,   ADM_kall,ADM_lall   ) ! upwind direction mask
    real(8) :: cmask_pl (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: GRD_xc   (ADM_gall,   ADM_kall,ADM_lall,   AI:AJ,XDIR:ZDIR) ! mass centroid position
    real(8) :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,XDIR:ZDIR)

    real(8), parameter :: b1 = 0.D0
    real(8), parameter :: b2 = 1.D0
    real(8), parameter :: b3 = 1.D0 - (b1+b2)

    integer :: nstart, nend
    integer :: n, k, l, v, iq

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 1st
    !---------------------------------------------------------------------------
    call DEBUG_rapstart('++++Vertical Adv.')

    do l = 1, ADM_lall
       rrhog(:,:,l) = 1.D0 / rhog_in(:,:,l)

       d(:,:,l) = b1 * frhog(:,:,l) * rrhog(:,:,l) * dt

       do k = ADM_kmin+1, ADM_kmax
       do n = 1, ADM_gall
          flx_v(n,k,l) = ( ( VMTR_C2Wfact_Gz(1,n,k,l) * rhogvx_mean(n,k  ,l) &
                           + VMTR_C2Wfact_Gz(2,n,k,l) * rhogvx_mean(n,k-1,l) &
                           + VMTR_C2Wfact_Gz(3,n,k,l) * rhogvy_mean(n,k  ,l) &
                           + VMTR_C2Wfact_Gz(4,n,k,l) * rhogvy_mean(n,k-1,l) &
                           + VMTR_C2Wfact_Gz(5,n,k,l) * rhogvz_mean(n,k  ,l) &
                           + VMTR_C2Wfact_Gz(6,n,k,l) * rhogvz_mean(n,k-1,l) &
                           ) * VMTR_RGAMH(n,k,l)                             & ! horizontal contribution
                         + rhogw_mean(n,k,l) * VMTR_RGSH(n,k,l)              & ! vertical   contribution
                         ) * 0.5D0 * dt
       enddo
       enddo
       flx_v(:,ADM_kmin,  l) = 0.D0
       flx_v(:,ADM_kmax+1,l) = 0.D0

       do k = ADM_kmin, ADM_kmax
          ck(:,k,l,1) = -flx_v(:,k  ,l) * rrhog(:,k,l) * GRD_rdgz(k)
          ck(:,k,l,2) =  flx_v(:,k+1,l) * rrhog(:,k,l) * GRD_rdgz(k)
       enddo
       ck(:,ADM_kmin-1,l,1) = 0.D0
       ck(:,ADM_kmin-1,l,2) = 0.D0
       ck(:,ADM_kmax+1,l,1) = 0.D0
       ck(:,ADM_kmax+1,l,2) = 0.D0
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          rrhog_pl(:,:,l) = 1.D0 / rhog_in_pl(:,:,l)

          d_pl(:,:,l) = b1 * frhog_pl(:,:,l) * rrhog_pl(:,:,l) * dt

          do k = ADM_kmin+1, ADM_kmax
          do n = 1, ADM_gall_pl
             flx_v_pl(n,k,l) = ( ( VMTR_C2Wfact_Gz_pl(1,n,k,l) * rhogvx_mean_pl(n,k  ,l) &
                                 + VMTR_C2Wfact_Gz_pl(2,n,k,l) * rhogvx_mean_pl(n,k-1,l) &
                                 + VMTR_C2Wfact_Gz_pl(3,n,k,l) * rhogvy_mean_pl(n,k  ,l) &
                                 + VMTR_C2Wfact_Gz_pl(4,n,k,l) * rhogvy_mean_pl(n,k-1,l) &
                                 + VMTR_C2Wfact_Gz_pl(5,n,k,l) * rhogvz_mean_pl(n,k  ,l) &
                                 + VMTR_C2Wfact_Gz_pl(6,n,k,l) * rhogvz_mean_pl(n,k-1,l) &
                                 ) * VMTR_RGAMH_pl(n,k,l)                                & ! horizontal contribution
                               + rhogw_mean_pl(n,k,l) * VMTR_RGSH_pl(n,k,l)              & ! vertical   contribution
                               ) * 0.5D0 * dt
          enddo
          enddo
          flx_v_pl(:,ADM_kmin,  l) = 0.D0
          flx_v_pl(:,ADM_kmax+1,l) = 0.D0

          do k = ADM_kmin, ADM_kmax
             ck_pl(:,k,l,1) = -flx_v_pl(:,k  ,l) * rrhog_pl(:,k,l) * GRD_rdgz(k)
             ck_pl(:,k,l,2) =  flx_v_pl(:,k+1,l) * rrhog_pl(:,k,l) * GRD_rdgz(k)
          enddo
          ck_pl(:,ADM_kmin-1,l,1) = 0.D0
          ck_pl(:,ADM_kmin-1,l,2) = 0.D0
          ck_pl(:,ADM_kmax+1,l,1) = 0.D0
          ck_pl(:,ADM_kmax+1,l,2) = 0.D0
       enddo
    endif

    !--- vertical advection: 2nd-order centered difference
    do iq = 1, vmax

       do l = 1, ADM_lall
          q(:,:,l) = rhogq(:,:,l,iq) * rrhog(:,:,l)

          do k = ADM_kmin, ADM_kmax+1
             q_h(:,k,l) = 0.5D0 * ( GRD_afac(k) * q(:,k,  l) &
                                  + GRD_bfac(k) * q(:,k-1,l) )
          enddo
          q_h(:,ADM_kmin-1,l) = 0.D0
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             q_pl(:,:,l) = rhogq_pl(:,:,l,iq) * rrhog_pl(:,:,l)

             do k = ADM_kmin, ADM_kmax+1
                q_h_pl(:,k,l) = 0.5D0 * ( GRD_afac(k) * q_pl(:,k,  l) &
                                        + GRD_bfac(k) * q_pl(:,k-1,l) )
             enddo
             q_h_pl(:,ADM_kmin-1,l) = 0.D0
          enddo
       endif

       if ( thubern_lim ) then
          call vertical_limiter_thuburn( q_h(:,:,:),   q_h_pl(:,:,:),  & ! [INOUT]
                                         q  (:,:,:),   q_pl  (:,:,:),  & ! [IN]
                                         d  (:,:,:),   d_pl  (:,:,:),  & ! [IN]
                                         ck (:,:,:,:), ck_pl (:,:,:,:) ) ! [IN]
       endif

       !--- update rhogq
       do l = 1, ADM_lall
          q_h(:,ADM_kmin  ,l) = 0.D0
          q_h(:,ADM_kmax+1,l) = 0.D0

          do k = ADM_kmin, ADM_kmax
          do n = 1, ADM_gall
             rhogq(n,k,l,iq) = rhogq(n,k,l,iq) - ( flx_v(n,k+1,l)*q_h(n,k+1,l) &
                                                 - flx_v(n,k,  l)*q_h(n,k,  l) &
                                                 ) * GRD_rdgz(k)
          enddo
          enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             q_h_pl(:,ADM_kmin  ,l) = 0.D0
             q_h_pl(:,ADM_kmax+1,l) = 0.D0

             do k = ADM_kmin, ADM_kmax
             do n = 1, ADM_gall_pl
                rhogq_pl(n,k,l,iq) = rhogq_pl(n,k,l,iq) - ( flx_v_pl(n,k+1,l)*q_h_pl(n,k+1,l) &
                                                          - flx_v_pl(n,k,  l)*q_h_pl(n,k,  l) &
                                                          ) * GRD_rdgz(k)
             enddo
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

    !--- update rhog
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
       do n = 1, ADM_gall
          rhog(n,k,l) = rhog_in(n,k,l) - ( flx_v(n,k+1,l) &
                                         - flx_v(n,k  ,l) &
                                         ) * GRD_rdgz(k)  &
                                       + b1 * frhog(n,k,l) * dt
       enddo
       enddo

       do k = ADM_kmin, ADM_kmax
          rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
          rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
       enddo
       rhog(:,ADM_kmin-1,l) = rhog_in(:,ADM_kmin,l)
       rhog(:,ADM_kmax+1,l) = rhog_in(:,ADM_kmax,l)
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
          do n = 1, ADM_gall_pl
             rhog_pl(n,k,l) = rhog_in_pl(n,k,l) - ( flx_v_pl(n,k+1,l) &
                                                  - flx_v_pl(n,k  ,l) &
                                                  ) * GRD_rdgz(k)     &
                                                + b1 * frhog_pl(n,k,l) * dt
          enddo
          enddo

          rhog_pl(:,ADM_kmin-1,l) = rhog_in_pl(:,ADM_kmin,l)
          rhog_pl(:,ADM_kmax+1,l) = rhog_in_pl(:,ADM_kmax,l)
       enddo
    endif

    call DEBUG_rapend('++++Vertical Adv.')
    !---------------------------------------------------------------------------
    ! Horizontal advection by MIURA scheme
    !---------------------------------------------------------------------------
    call DEBUG_rapstart('++++Horizontal Adv.')

    do l = 1, ADM_lall
       rrhog(:,:,l) = 1.D0 / rhog(:,:,l)

       d(:,:,l) = b2 * frhog(:,:,l) * rrhog(:,:,l) * dt

       rhogvx(:,:,l) = rhogvx_mean(:,:,l) * VMTR_RGAM(:,:,l)
       rhogvy(:,:,l) = rhogvy_mean(:,:,l) * VMTR_RGAM(:,:,l)
       rhogvz(:,:,l) = rhogvz_mean(:,:,l) * VMTR_RGAM(:,:,l)
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          rrhog_pl(:,:,l) = 1.D0 / rhog_pl(:,:,l)

          d_pl(:,:,l) = b2 * frhog_pl(:,:,l) * rrhog_pl(:,:,l) * dt

          rhogvx_pl(:,:,l) = rhogvx_mean_pl(:,:,l) * VMTR_RGAM_pl(:,:,l)
          rhogvy_pl(:,:,l) = rhogvy_mean_pl(:,:,l) * VMTR_RGAM_pl(:,:,l)
          rhogvz_pl(:,:,l) = rhogvz_mean_pl(:,:,l) * VMTR_RGAM_pl(:,:,l)
       enddo
    endif

    call horizontal_flux( flx_h    (:,:,:,:),   flx_h_pl    (:,:,:),   & ! [OUT]
                          GRD_xc   (:,:,:,:,:), GRD_xc_pl   (:,:,:,:), & ! [OUT]
                          rhog_mean(:,:,:),     rhog_mean_pl(:,:,:),   & ! [IN]
                          rhogvx   (:,:,:),     rhogvx_pl   (:,:,:),   & ! [IN]
                          rhogvy   (:,:,:),     rhogvy_pl   (:,:,:),   & ! [IN]
                          rhogvz   (:,:,:),     rhogvz_pl   (:,:,:),   & ! [IN]
                          dt                                           ) ! [IN]

    !--- Courant number
    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_gall
       ch(1,n,k,l) = flx_h(1,n,k,l) * rrhog(n,k,l)
       ch(2,n,k,l) = flx_h(2,n,k,l) * rrhog(n,k,l)
       ch(3,n,k,l) = flx_h(3,n,k,l) * rrhog(n,k,l)
       ch(4,n,k,l) = flx_h(4,n,k,l) * rrhog(n,k,l)
       ch(5,n,k,l) = flx_h(5,n,k,l) * rrhog(n,k,l)
       ch(6,n,k,l) = flx_h(6,n,k,l) * rrhog(n,k,l)

       ! c <= 0(incoming), cmask = 1
       cmask(1,n,k,l) = 0.5D0 - sign(0.5D0,ch(1,n,k,l)-CNST_EPS_ZERO)
       cmask(2,n,k,l) = 0.5D0 - sign(0.5D0,ch(2,n,k,l)-CNST_EPS_ZERO)
       cmask(3,n,k,l) = 0.5D0 - sign(0.5D0,ch(3,n,k,l)-CNST_EPS_ZERO)
       cmask(4,n,k,l) = 0.5D0 - sign(0.5D0,ch(4,n,k,l)-CNST_EPS_ZERO)
       cmask(5,n,k,l) = 0.5D0 - sign(0.5D0,ch(5,n,k,l)-CNST_EPS_ZERO)
       cmask(6,n,k,l) = 0.5D0 - sign(0.5D0,ch(6,n,k,l)-CNST_EPS_ZERO)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do v = ADM_gmin_pl, ADM_gmax_pl
          ch_pl(v,k,l) = flx_h_pl(v,k,l) * rrhog_pl(n,k,l)

          cmask_pl(v,k,l) = 0.5D0 - sign(0.5D0,ch_pl(v,k,l)-CNST_EPS_ZERO)
       enddo
       enddo
       enddo
    endif

    do iq = 1, vmax

       q(:,:,:) = rhogq(:,:,:,iq) * rrhog(:,:,:)
       if ( ADM_have_pl ) then
          q_pl(:,:,:) = rhogq_pl(:,:,:,iq) * rrhog_pl(:,:,:)
       endif

       ! calculate q at cell face, upwind side
       call horizontal_remap( q_a   (:,:,:,:),   q_a_pl   (:,:,:),   & ! [OUT]
                              q     (:,:,:),     q_pl     (:,:,:),   & ! [IN]
                              cmask (:,:,:,:),   cmask_pl (:,:,:),   & ! [IN]
                              GRD_xc(:,:,:,:,:), GRD_xc_pl(:,:,:,:)  ) ! [IN]

       ! apply flux limiter
       if ( thubern_lim ) then
          call horizontal_limiter_thuburn( q_a  (:,:,:,:),   q_a_pl  (:,:,:), & ! [INOUT]
                                           q    (:,:,:),     q_pl    (:,:,:), & ! [IN]
                                           d    (:,:,:),     d_pl    (:,:,:), & ! [IN]
                                           ch   (:,:,:,:),   ch_pl   (:,:,:), & ! [IN]
                                           cmask(:,:,:,:),   cmask_pl(:,:,:)  ) ! [IN]
       endif

       !--- update rhogq
       do l = 1, ADM_lall
          nstart = suf(ADM_gmin,ADM_gmin)
          nend   = suf(ADM_gmax,ADM_gmax)

          do k = 1, ADM_kall
          do n = nstart, nend
             rhogq(n,k,l,iq) = rhogq(n,k,l,iq) - ( flx_h(1,n,k,l) * q_a(1,n,k,l) &
                                                 + flx_h(2,n,k,l) * q_a(2,n,k,l) &
                                                 + flx_h(3,n,k,l) * q_a(3,n,k,l) &
                                                 + flx_h(4,n,k,l) * q_a(4,n,k,l) &
                                                 + flx_h(5,n,k,l) * q_a(5,n,k,l) &
                                                 + flx_h(6,n,k,l) * q_a(6,n,k,l) )
          enddo
          enddo

          do k = 1, ADM_kall
             rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
             rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
          enddo
       enddo

       if ( ADM_have_pl ) then
          n = ADM_gslf_pl

          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do v = ADM_gmin_pl, ADM_gmax_pl
             rhogq_pl(n,k,l,iq) = rhogq_pl(n,k,l,iq) - flx_h_pl(v,k,l) * q_a_pl(v,k,l)
          enddo
          enddo
          enddo
       endif

    enddo ! tracer q LOOP

    !--- update rhog
    do l = 1, ADM_lall
       nstart = suf(ADM_gmin,ADM_gmin)
       nend   = suf(ADM_gmax,ADM_gmax)

       do k = 1, ADM_kall
       do n = nstart, nend
          rhog(n,k,l)= rhog(n,k,l) - ( flx_h(1,n,k,l) &
                                     + flx_h(2,n,k,l) &
                                     + flx_h(3,n,k,l) &
                                     + flx_h(4,n,k,l) &
                                     + flx_h(5,n,k,l) &
                                     + flx_h(6,n,k,l) ) + b2 * frhog(n,k,l) * dt
       enddo
       enddo

       do k = 1, ADM_kall
          rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
          rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
       enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          do v = ADM_gmin_pl, ADM_gmax_pl
             rhog_pl(n,k,l)= rhog_pl(n,k,l) - flx_h_pl(v,k,l)
          enddo
          rhog_pl(n,k,l)= rhog_pl(n,k,l) + b2 * frhog_pl(n,k,l) * dt
       enddo
       enddo
    endif

    call DEBUG_rapend('++++Horizontal Adv.')
    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 2nd
    !---------------------------------------------------------------------------
    call DEBUG_rapstart('++++Vertical Adv.')

    do l = 1, ADM_lall

       rrhog(:,:,l) = 1.D0 / rhog(:,:,l)

       d(:,:,l) = b3 * frhog(:,:,l) * dt * rrhog(:,:,l)

       do k = ADM_kmin, ADM_kmax
          ck(:,k,l,1) = -flx_v(:,k  ,l) * rrhog(:,k,l) * GRD_rdgz(k)
          ck(:,k,l,2) =  flx_v(:,k+1,l) * rrhog(:,k,l) * GRD_rdgz(k)
       enddo
       ck(:,ADM_kmin-1,l,1) = 0.D0
       ck(:,ADM_kmin-1,l,2) = 0.D0
       ck(:,ADM_kmax+1,l,1) = 0.D0
       ck(:,ADM_kmax+1,l,2) = 0.D0
    enddo ! l LOOP

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl

          rrhog_pl(:,:,l) = 1.D0 / rhog_pl(:,:,l)

          d_pl(:,:,l) = b3 * frhog_pl(:,:,l) * dt * rrhog_pl(:,:,l)

          do k = ADM_kmin, ADM_kmax
             ck_pl(:,k,l,1) = -flx_v_pl(:,k  ,l) * rrhog_pl(:,k,l) * GRD_rdgz(k)
             ck_pl(:,k,l,2) =  flx_v_pl(:,k+1,l) * rrhog_pl(:,k,l) * GRD_rdgz(k)
          enddo
          ck_pl(:,ADM_kmin-1,l,1) = 0.D0
          ck_pl(:,ADM_kmin-1,l,2) = 0.D0
          ck_pl(:,ADM_kmax+1,l,1) = 0.D0
          ck_pl(:,ADM_kmax+1,l,2) = 0.D0
       enddo
    endif

    !--- basic scheme ( 2nd-order centered difference )
    do iq = 1, vmax

       do l = 1, ADM_lall
          q(:,:,l) = rhogq(:,:,l,iq) * rrhog(:,:,l)

          do k = ADM_kmin, ADM_kmax+1
             q_h(:,k,l) = 0.5D0 * ( GRD_afac(k) * q(:,k,  l) &
                                  + GRD_bfac(k) * q(:,k-1,l) )
          enddo
          q_h(:,ADM_kmin-1,l) = 0.D0
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             q_pl(:,:,l) = rhogq_pl(:,:,l,iq) * rrhog_pl(:,:,l)

             do k = ADM_kmin, ADM_kmax+1
                q_h_pl(:,k,l) = 0.5D0 * ( GRD_afac(k) * q_pl(:,k,  l) &
                                        + GRD_bfac(k) * q_pl(:,k-1,l) )
             enddo
             q_h_pl(:,ADM_kmin-1,l) = 0.D0
          enddo
       endif

       if ( thubern_lim ) then
          call vertical_limiter_thuburn( q_h(:,:,:),   q_h_pl(:,:,:),   & ! [INOUT]
                                         q  (:,:,:),   q_pl  (:,:,:),   & ! [IN]
                                         d  (:,:,:),   d_pl  (:,:,:),   & ! [IN]
                                         ck (:,:,:,:), ck_pl (:,:,:,:)  ) ! [IN]
       endif

       !--- update rhogq
       do l = 1, ADM_lall
          q_h(:,ADM_kmin  ,l) = 0.D0
          q_h(:,ADM_kmax+1,l) = 0.D0

          do k = ADM_kmin, ADM_kmax
          do n = 1, ADM_gall
             rhogq(n,k,l,iq) = rhogq(n,k,l,iq) - ( flx_v(n,k+1,l)*q_h(n,k+1,l) &
                                                 - flx_v(n,k,  l)*q_h(n,k,  l) &
                                                 ) * GRD_rdgz(k)
          enddo
          enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             q_h_pl(:,ADM_kmin  ,l) = 0.D0
             q_h_pl(:,ADM_kmax+1,l) = 0.D0

             do k = ADM_kmin, ADM_kmax
             do n = 1, ADM_gall_pl
                rhogq_pl(n,k,l,iq) = rhogq_pl(n,k,l,iq) - ( flx_v_pl(n,k+1,l)*q_h_pl(n,k+1,l) &
                                                          - flx_v_pl(n,k,  l)*q_h_pl(n,k,  l) &
                                                          ) * GRD_rdgz(k)
             enddo
             enddo
          enddo
       endif
    enddo ! tracer q LOOP

    call DEBUG_rapend('++++Vertical Adv.')

    return
  end subroutine src_tracer_advection

  !-----------------------------------------------------------------------------
  !> prepare horizontal advection trem: mass flux, GRD_xc
  subroutine horizontal_flux( &
       flx_h,  flx_h_pl,  &
       GRD_xc, GRD_xc_pl, &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )
    use mod_adm, only: &
       ADM_have_pl,    &
       ADM_have_sgp,   &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl,    &
       ADM_kmin,       &
       ADM_kmax
    use mod_cnst, only: &
       CNST_UNDEF, &
       CNST_EPS_ZERO
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_xt,   &
       GRD_xt_pl
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    use mod_oprt, only: &
       cinterp_HN,    &
       cinterp_PRA
    implicit none

    real(8), intent(out) :: flx_h    (6,ADM_gall   ,ADM_kall,ADM_lall   )               ! horizontal mass flux
    real(8), intent(out) :: flx_h_pl (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,AI:AJ,XDIR:ZDIR) ! mass centroid position
    real(8), intent(out) :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)
    real(8), intent(in)  :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )                 ! rho at cell center
    real(8), intent(in)  :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: dt

    real(8) :: rhot     (ADM_gall   ,TI:TJ) ! rho at cell vertex
    real(8) :: rhot_pl  (ADM_gall_pl)
    real(8) :: rhovxt   (ADM_gall   ,TI:TJ)
    real(8) :: rhovxt_pl(ADM_gall_pl)
    real(8) :: rhovyt   (ADM_gall   ,TI:TJ)
    real(8) :: rhovyt_pl(ADM_gall_pl)
    real(8) :: rhovzt   (ADM_gall   ,TI:TJ)
    real(8) :: rhovzt_pl(ADM_gall_pl)

    real(8) :: flux
    real(8) :: rrhoa2

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1

    integer :: nstart,nend
    integer :: n, k, l, v

    logical, save :: first = .true.  ! Y.Niwa add 080130

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    if ( first ) then
       first = .false.

       allocate( local_t_var   (ADM_gall   ,k0,ADM_lall   ,TI:TJ,W1:W3) )
       allocate( local_t_var_pl(ADM_gall_pl,k0,ADM_lall_pl,      W1:W3) )
       local_t_var   (:,:,:,:,:) = CNST_UNDEF
       local_t_var_pl(:,:,:,:)   = CNST_UNDEF

       allocate( GRD_xr   (ADM_gall   ,k0,ADM_lall   ,AI:AJ,XDIR:ZDIR) )
       allocate( GRD_xr_pl(ADM_gall_pl,k0,ADM_lall_pl,      XDIR:ZDIR) )
       GRD_xr   (:,:,:,:,:) = CNST_UNDEF
       GRD_xr_pl(:,:,:,:)   = CNST_UNDEF

       do l = 1, ADM_lall
          do n = 1, ADM_gall
             local_t_var(n,k0,l,TI,W1) = GMTR_t_var(n,k0,l,TI,W1)
             local_t_var(n,k0,l,TI,W2) = GMTR_t_var(n,k0,l,TI,W2)
             local_t_var(n,k0,l,TI,W3) = GMTR_t_var(n,k0,l,TI,W3)

             local_t_var(n,k0,l,TJ,W1) = GMTR_t_var(n,k0,l,TJ,W1)
             local_t_var(n,k0,l,TJ,W2) = GMTR_t_var(n,k0,l,TJ,W2)
             local_t_var(n,k0,l,TJ,W3) = GMTR_t_var(n,k0,l,TJ,W3)
          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )

          do n = nstart, nend
             ij     = n
             ijm1   = n     - ADM_gall_1d

             GRD_xr(n,k0,l,AI ,XDIR) = 0.5D0 * ( GRD_xt(ijm1,k0,l,TJ,XDIR) + GRD_xt(ij,k0,l,TI,XDIR) )
             GRD_xr(n,k0,l,AI ,YDIR) = 0.5D0 * ( GRD_xt(ijm1,k0,l,TJ,YDIR) + GRD_xt(ij,k0,l,TI,YDIR) )
             GRD_xr(n,k0,l,AI ,ZDIR) = 0.5D0 * ( GRD_xt(ijm1,k0,l,TJ,ZDIR) + GRD_xt(ij,k0,l,TI,ZDIR) )
          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )

          do n = nstart, nend
             ij     = n

             GRD_xr(n,k0,l,AIJ,XDIR) = 0.5D0 * ( GRD_xt(ij,k0,l,TI,XDIR) + GRD_xt(ij,k0,l,TJ,XDIR) )
             GRD_xr(n,k0,l,AIJ,YDIR) = 0.5D0 * ( GRD_xt(ij,k0,l,TI,YDIR) + GRD_xt(ij,k0,l,TJ,YDIR) )
             GRD_xr(n,k0,l,AIJ,ZDIR) = 0.5D0 * ( GRD_xt(ij,k0,l,TI,ZDIR) + GRD_xt(ij,k0,l,TJ,ZDIR) )
          enddo

          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )

          do n = nstart, nend
             ij     = n
             im1j   = n - 1

             GRD_xr(n,k0,l,AJ ,XDIR) = 0.5D0 * ( GRD_xt(ij,k0,l,TJ,XDIR) + GRD_xt(im1j,k0,l,TI,XDIR) )
             GRD_xr(n,k0,l,AJ ,YDIR) = 0.5D0 * ( GRD_xt(ij,k0,l,TJ,YDIR) + GRD_xt(im1j,k0,l,TI,YDIR) )
             GRD_xr(n,k0,l,AJ ,ZDIR) = 0.5D0 * ( GRD_xt(ij,k0,l,TJ,ZDIR) + GRD_xt(im1j,k0,l,TI,ZDIR) )
          enddo
       enddo

       if ( ADM_have_pl ) then
          n = ADM_gslf_pl

          do l = 1, ADM_lall_pl
             do v = ADM_gmin_pl, ADM_gmax_pl
                local_t_var_pl(v,k0,l,W1) = GMTR_t_var_pl(v,k0,l,W1)
                local_t_var_pl(v,k0,l,W2) = GMTR_t_var_pl(v,k0,l,W2)
                local_t_var_pl(v,k0,l,W3) = GMTR_t_var_pl(v,k0,l,W3)
             enddo

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijm1 = v - 1
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                GRD_xr_pl(v,k0,l,XDIR) = 0.5D0 * (GRD_xt_pl(ijm1,k0,l,XDIR)+GRD_xt_pl(ij,k0,l,XDIR))
                GRD_xr_pl(v,k0,l,YDIR) = 0.5D0 * (GRD_xt_pl(ijm1,k0,l,YDIR)+GRD_xt_pl(ij,k0,l,YDIR))
                GRD_xr_pl(v,k0,l,ZDIR) = 0.5D0 * (GRD_xt_pl(ijm1,k0,l,ZDIR)+GRD_xt_pl(ij,k0,l,ZDIR))
             enddo
          enddo
       endif
    endif

    do l = 1, ADM_lall
    do k = 1, ADM_kall

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d

          rhot  (n,TI) = rho  (ij    ,k,l) * local_t_var(n,k0,l,TI,W1) &
                       + rho  (ip1j  ,k,l) * local_t_var(n,k0,l,TI,W2) &
                       + rho  (ip1jp1,k,l) * local_t_var(n,k0,l,TI,W3)
          rhovxt(n,TI) = rhovx(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                       + rhovx(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                       + rhovx(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
          rhovyt(n,TI) = rhovy(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                       + rhovy(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                       + rhovy(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
          rhovzt(n,TI) = rhovz(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                       + rhovz(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                       + rhovz(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
       enddo

       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          rhot  (n,TJ) = rho  (ij    ,k,l) * local_t_var(n,k0,l,TJ,W1) &
                       + rho  (ip1jp1,k,l) * local_t_var(n,k0,l,TJ,W2) &
                       + rho  (ijp1  ,k,l) * local_t_var(n,k0,l,TJ,W3)
          rhovxt(n,TJ) = rhovx(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                       + rhovx(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                       + rhovx(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
          rhovyt(n,TJ) = rhovy(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                       + rhovy(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                       + rhovy(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
          rhovzt(n,TJ) = rhovz(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                       + rhovz(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                       + rhovz(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
       enddo

       if ( ADM_have_sgp(l) ) then
          rhot  (suf(ADM_gmin-1,ADM_gmin-1),TI) = rhot  (suf(ADM_gmin,ADM_gmin-1),TJ)
          rhovxt(suf(ADM_gmin-1,ADM_gmin-1),TI) = rhovxt(suf(ADM_gmin,ADM_gmin-1),TJ)
          rhovyt(suf(ADM_gmin-1,ADM_gmin-1),TI) = rhovyt(suf(ADM_gmin,ADM_gmin-1),TJ)
          rhovzt(suf(ADM_gmin-1,ADM_gmin-1),TI) = rhovzt(suf(ADM_gmin,ADM_gmin-1),TJ)
       endif

       !--- calculate flux and mass centroid position

       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ijm1   = n     - ADM_gall_1d
          ip1j   = n + 1

          flux = 0.5D0 * ( (rhovxt(ijm1,TJ)+rhovxt(ij,TI)) * cinterp_HN(AI ,1,ij,l) &
                         + (rhovyt(ijm1,TJ)+rhovyt(ij,TI)) * cinterp_HN(AI ,2,ij,l) &
                         + (rhovzt(ijm1,TJ)+rhovzt(ij,TI)) * cinterp_HN(AI ,3,ij,l) )

          flx_h(1,ij  ,k,l) =  flux * cinterp_PRA(ij  ,l) * dt
          flx_h(4,ip1j,k,l) = -flux * cinterp_PRA(ip1j,l) * dt
       enddo

       do n = nstart, nend
          ij     = n
          ijm1   = n     - ADM_gall_1d

          rrhoa2 = 1.D0 / max( rhot(ijm1,TJ) + rhot(ij,TI), CNST_EPS_ZERO ) ! doubled

          GRD_xc(n,k,l,AI,XDIR) = GRD_xr(n,k0,l,AI,XDIR) - (rhovxt(ijm1,TJ)+rhovxt(ij,TI)) * rrhoa2 * dt * 0.5D0
          GRD_xc(n,k,l,AI,YDIR) = GRD_xr(n,k0,l,AI,YDIR) - (rhovyt(ijm1,TJ)+rhovyt(ij,TI)) * rrhoa2 * dt * 0.5D0
          GRD_xc(n,k,l,AI,ZDIR) = GRD_xr(n,k0,l,AI,ZDIR) - (rhovzt(ijm1,TJ)+rhovzt(ij,TI)) * rrhoa2 * dt * 0.5D0
       enddo

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ip1jp1 = n + 1 + ADM_gall_1d

          flux = 0.5D0 * ( (rhovxt(ij,TI)+rhovxt(ij,TJ)) * cinterp_HN(AIJ,1,ij,l) &
                         + (rhovyt(ij,TI)+rhovyt(ij,TJ)) * cinterp_HN(AIJ,2,ij,l) &
                         + (rhovzt(ij,TI)+rhovzt(ij,TJ)) * cinterp_HN(AIJ,3,ij,l) )

          flx_h(2,ij    ,k,l) =  flux * cinterp_PRA(ij    ,l) * dt
          flx_h(5,ip1jp1,k,l) = -flux * cinterp_PRA(ip1jp1,l) * dt
       enddo

       do n = nstart, nend
          ij     = n

          rrhoa2 = 1.D0 / max( rhot(ij,TI) + rhot(ij,TJ), CNST_EPS_ZERO ) ! doubled

          GRD_xc(n,k,l,AIJ,XDIR) = GRD_xr(n,k0,l,AIJ,XDIR) - (rhovxt(ij,TI)+rhovxt(ij,TJ)) * rrhoa2 * dt * 0.5D0
          GRD_xc(n,k,l,AIJ,YDIR) = GRD_xr(n,k0,l,AIJ,YDIR) - (rhovyt(ij,TI)+rhovyt(ij,TJ)) * rrhoa2 * dt * 0.5D0
          GRD_xc(n,k,l,AIJ,ZDIR) = GRD_xr(n,k0,l,AIJ,ZDIR) - (rhovzt(ij,TI)+rhovzt(ij,TJ)) * rrhoa2 * dt * 0.5D0
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          im1j   = n - 1
          ijp1   = n     + ADM_gall_1d

          flux = 0.5D0 * ( (rhovxt(ij,TJ)+rhovxt(im1j,TI)) * cinterp_HN(AJ ,1,ij,l) &
                         + (rhovyt(ij,TJ)+rhovyt(im1j,TI)) * cinterp_HN(AJ ,2,ij,l) &
                         + (rhovzt(ij,TJ)+rhovzt(im1j,TI)) * cinterp_HN(AJ ,3,ij,l) )

          flx_h(3,ij  ,k,l) =  flux * cinterp_PRA(ij  ,l) * dt
          flx_h(6,ijp1,k,l) = -flux * cinterp_PRA(ijp1,l) * dt
       enddo

       do n = nstart, nend
          ij     = n
          im1j   = n - 1

          rrhoa2 = 1.D0 / max( rhot(ij,TJ) + rhot(im1j,TI), CNST_EPS_ZERO ) ! doubled

          GRD_xc(n,k,l,AJ,XDIR) = GRD_xr(n,k0,l,AJ,XDIR) - (rhovxt(ij,TJ)+rhovxt(im1j,TI)) * rrhoa2 * dt * 0.5D0
          GRD_xc(n,k,l,AJ,YDIR) = GRD_xr(n,k0,l,AJ,YDIR) - (rhovyt(ij,TJ)+rhovyt(im1j,TI)) * rrhoa2 * dt * 0.5D0
          GRD_xc(n,k,l,AJ,ZDIR) = GRD_xr(n,k0,l,AJ,ZDIR) - (rhovzt(ij,TJ)+rhovzt(im1j,TI)) * rrhoa2 * dt * 0.5D0
       enddo

       if ( ADM_have_sgp(l) ) then
          flx_h(6,suf(ADM_gmin,ADM_gmin),k,l) = 0.D0
       endif

    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl

             rhot_pl  (v) = rho_pl  (n   ,k,l) * local_T_var_pl(ij,k0,l,W1) &
                          + rho_pl  (ij  ,k,l) * local_T_var_pl(ij,k0,l,W2) &
                          + rho_pl  (ijp1,k,l) * local_T_var_pl(ij,k0,l,W3)
             rhovxt_pl(v) = rhovx_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                          + rhovx_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                          + rhovx_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)
             rhovyt_pl(v) = rhovy_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                          + rhovy_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                          + rhovy_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)
             rhovzt_pl(v) = rhovz_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                          + rhovz_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                          + rhovz_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             flux = 0.5D0 * ( (rhovxt_pl(ijm1)+rhovxt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HNX) &
                            + (rhovyt_pl(ijm1)+rhovyt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HNY) &
                            + (rhovzt_pl(ijm1)+rhovzt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HNZ) )

             flx_h_pl(v,k,l) = flux * GMTR_P_var_pl(n,k0,l,P_RAREA) * dt
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             rrhoa2 = 1.D0 / max( rhot_pl(ijm1) + rhot_pl(ij), CNST_EPS_ZERO ) ! doubled

             GRD_xc_pl(v,k,l,XDIR) = GRD_xr_pl(v,k0,l,XDIR) - (rhovxt_pl(ijm1)+rhovxt_pl(ij)) * rrhoa2 * dt * 0.5D0
             GRD_xc_pl(v,k,l,YDIR) = GRD_xr_pl(v,k0,l,YDIR) - (rhovyt_pl(ijm1)+rhovyt_pl(ij)) * rrhoa2 * dt * 0.5D0
             GRD_xc_pl(v,k,l,ZDIR) = GRD_xr_pl(v,k0,l,ZDIR) - (rhovzt_pl(ijm1)+rhovzt_pl(ij)) * rrhoa2 * dt * 0.5D0
          enddo

       enddo
       enddo
    endif

    return
  end subroutine horizontal_flux

  !-----------------------------------------------------------------------------
  !> local linear approximation of q (Miura, 2007)
  subroutine horizontal_remap( &
       q_a,    q_a_pl,   &
       q,      q_pl,     &
       cmask,  cmask_pl, &
       GRD_xc, GRD_xc_pl )
    use mod_adm, only: &
       ADM_have_pl,    &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_x,   &
       GRD_x_pl
    use mod_oprt, only: &
       OPRT_gradient
    implicit none

    real(8), intent(out) :: q_a      (6,ADM_gall   ,ADM_kall,ADM_lall   )               ! q at cell face
    real(8), intent(out) :: q_a_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: q        (  ADM_gall   ,ADM_kall,ADM_lall   )               ! q at cell center
    real(8), intent(in)  :: q_pl     (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: cmask    (6,ADM_gall   ,ADM_kall,ADM_lall   )               ! upwind direction mask
    real(8), intent(in)  :: cmask_pl (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,AI:AJ,XDIR:ZDIR) ! position of the mass centroid
    real(8), intent(in)  :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)

    real(8)  :: gradq   (ADM_gall   ,ADM_kall,ADM_lall   ,XDIR:ZDIR) ! grad(q)
    real(8)  :: gradq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,XDIR:ZDIR)

    real(8) :: q_ap, q_am

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: nstart,nend
    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call OPRT_gradient( q    (:,:,:),   q_pl    (:,:,:),  & ![IN]
                        gradq(:,:,:,:), gradq_pl(:,:,:,:) ) ![OUT]

    call COMM_data_transfer( gradq(:,:,:,:), gradq_pl(:,:,:,:) )

    ! interpolated Q at cell arc
    do l = 1, ADM_lall
    do k = 1, ADM_kall
       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ip1j   = n + 1

          q_ap = q(ij    ,k,l) + gradq(ij    ,k,l,XDIR) * ( GRD_xc(ij    ,k,l,AI ,XDIR) - GRD_x(ij    ,k0,l,XDIR) ) &
                               + gradq(ij    ,k,l,YDIR) * ( GRD_xc(ij    ,k,l,AI ,YDIR) - GRD_x(ij    ,k0,l,YDIR) ) &
                               + gradq(ij    ,k,l,ZDIR) * ( GRD_xc(ij    ,k,l,AI ,ZDIR) - GRD_x(ij    ,k0,l,ZDIR) )

          q_am = q(ip1j  ,k,l) + gradq(ip1j  ,k,l,XDIR) * ( GRD_xc(ij    ,k,l,AI ,XDIR) - GRD_x(ip1j  ,k0,l,XDIR) ) &
                               + gradq(ip1j  ,k,l,YDIR) * ( GRD_xc(ij    ,k,l,AI ,YDIR) - GRD_x(ip1j  ,k0,l,YDIR) ) &
                               + gradq(ip1j  ,k,l,ZDIR) * ( GRD_xc(ij    ,k,l,AI ,ZDIR) - GRD_x(ip1j  ,k0,l,ZDIR) )

          q_a(1,n,k,l) = (      cmask(1,n,k,l) ) * q_am &
                       + ( 1.D0-cmask(1,n,k,l) ) * q_ap
       enddo

       do n = nstart, nend
          ij     = n
          ip1jp1 = n + 1 + ADM_gall_1d

          q_ap = q(ij    ,k,l) + gradq(ij    ,k,l,XDIR) * ( GRD_xc(ij    ,k,l,AIJ,XDIR) - GRD_x(ij    ,k0,l,XDIR) ) &
                               + gradq(ij    ,k,l,YDIR) * ( GRD_xc(ij    ,k,l,AIJ,YDIR) - GRD_x(ij    ,k0,l,YDIR) ) &
                               + gradq(ij    ,k,l,ZDIR) * ( GRD_xc(ij    ,k,l,AIJ,ZDIR) - GRD_x(ij    ,k0,l,ZDIR) )

          q_am = q(ip1jp1,k,l) + gradq(ip1jp1,k,l,XDIR) * ( GRD_xc(ij    ,k,l,AIJ,XDIR) - GRD_x(ip1jp1,k0,l,XDIR) ) &
                               + gradq(ip1jp1,k,l,YDIR) * ( GRD_xc(ij    ,k,l,AIJ,YDIR) - GRD_x(ip1jp1,k0,l,YDIR) ) &
                               + gradq(ip1jp1,k,l,ZDIR) * ( GRD_xc(ij    ,k,l,AIJ,ZDIR) - GRD_x(ip1jp1,k0,l,ZDIR) )

          q_a(2,n,k,l) = (      cmask(2,n,k,l) ) * q_am &
                       + ( 1.D0-cmask(2,n,k,l) ) * q_ap
       enddo

       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d

          q_ap = q(ij    ,k,l) + gradq(ij    ,k,l,XDIR) * ( GRD_xc(ij    ,k,l,AJ ,XDIR) - GRD_x(ij    ,k0,l,XDIR) ) &
                               + gradq(ij    ,k,l,YDIR) * ( GRD_xc(ij    ,k,l,AJ ,YDIR) - GRD_x(ij    ,k0,l,YDIR) ) &
                               + gradq(ij    ,k,l,ZDIR) * ( GRD_xc(ij    ,k,l,AJ ,ZDIR) - GRD_x(ij    ,k0,l,ZDIR) )

          q_am = q(ijp1  ,k,l) + gradq(ijp1  ,k,l,XDIR) * ( GRD_xc(ij    ,k,l,AJ ,XDIR) - GRD_x(ijp1  ,k0,l,XDIR) ) &
                               + gradq(ijp1  ,k,l,YDIR) * ( GRD_xc(ij    ,k,l,AJ ,YDIR) - GRD_x(ijp1  ,k0,l,YDIR) ) &
                               + gradq(ijp1  ,k,l,ZDIR) * ( GRD_xc(ij    ,k,l,AJ ,ZDIR) - GRD_x(ijp1  ,k0,l,ZDIR) )

          q_a(3,n,k,l) = (      cmask(3,n,k,l) ) * q_am &
                       + ( 1.D0-cmask(3,n,k,l) ) * q_ap
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          im1j   = n - 1

          q_ap = q(im1j  ,k,l) + gradq(im1j  ,k,l,XDIR) * ( GRD_xc(im1j  ,k,l,AI ,XDIR) - GRD_x(im1j  ,k0,l,XDIR) ) &
                               + gradq(im1j  ,k,l,YDIR) * ( GRD_xc(im1j  ,k,l,AI ,YDIR) - GRD_x(im1j  ,k0,l,YDIR) ) &
                               + gradq(im1j  ,k,l,ZDIR) * ( GRD_xc(im1j  ,k,l,AI ,ZDIR) - GRD_x(im1j  ,k0,l,ZDIR) )

          q_am = q(ij    ,k,l) + gradq(ij    ,k,l,XDIR) * ( GRD_xc(im1j  ,k,l,AI ,XDIR) - GRD_x(ij    ,k0,l,XDIR) ) &
                               + gradq(ij    ,k,l,YDIR) * ( GRD_xc(im1j  ,k,l,AI ,YDIR) - GRD_x(ij    ,k0,l,YDIR) ) &
                               + gradq(ij    ,k,l,ZDIR) * ( GRD_xc(im1j  ,k,l,AI ,ZDIR) - GRD_x(ij    ,k0,l,ZDIR) )

          q_a(4,n,k,l) = (      cmask(4,n,k,l) ) * q_am &
                       + ( 1.D0-cmask(4,n,k,l) ) * q_ap
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          im1jm1 = n - 1 - ADM_gall_1d

          q_ap = q(im1jm1,k,l) + gradq(im1jm1,k,l,XDIR) * ( GRD_xc(im1jm1,k,l,AIJ,XDIR) - GRD_x(im1jm1,k0,l,XDIR) ) &
                               + gradq(im1jm1,k,l,YDIR) * ( GRD_xc(im1jm1,k,l,AIJ,YDIR) - GRD_x(im1jm1,k0,l,YDIR) ) &
                               + gradq(im1jm1,k,l,ZDIR) * ( GRD_xc(im1jm1,k,l,AIJ,ZDIR) - GRD_x(im1jm1,k0,l,ZDIR) )

          q_am = q(ij    ,k,l) + gradq(ij    ,k,l,XDIR) * ( GRD_xc(im1jm1,k,l,AIJ,XDIR) - GRD_x(ij    ,k0,l,XDIR) ) &
                               + gradq(ij    ,k,l,YDIR) * ( GRD_xc(im1jm1,k,l,AIJ,YDIR) - GRD_x(ij    ,k0,l,YDIR) ) &
                               + gradq(ij    ,k,l,ZDIR) * ( GRD_xc(im1jm1,k,l,AIJ,ZDIR) - GRD_x(ij    ,k0,l,ZDIR) )

          q_a(5,n,k,l) = (      cmask(5,n,k,l) ) * q_am &
                       + ( 1.D0-cmask(5,n,k,l) ) * q_ap
       enddo

       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ijm1   = n     - ADM_gall_1d

          q_ap = q(ijm1  ,k,l) + gradq(ijm1  ,k,l,XDIR) * ( GRD_xc(ijm1,k,l,AJ ,XDIR) - GRD_x(ijm1  ,k0,l,XDIR) ) &
                               + gradq(ijm1  ,k,l,YDIR) * ( GRD_xc(ijm1,k,l,AJ ,YDIR) - GRD_x(ijm1  ,k0,l,YDIR) ) &
                               + gradq(ijm1  ,k,l,ZDIR) * ( GRD_xc(ijm1,k,l,AJ ,ZDIR) - GRD_x(ijm1  ,k0,l,ZDIR) )

          q_am = q(ij    ,k,l) + gradq(ij    ,k,l,XDIR) * ( GRD_xc(ijm1,k,l,AJ ,XDIR) - GRD_x(ij    ,k0,l,XDIR) ) &
                               + gradq(ij    ,k,l,YDIR) * ( GRD_xc(ijm1,k,l,AJ ,YDIR) - GRD_x(ij    ,k0,l,YDIR) ) &
                               + gradq(ij    ,k,l,ZDIR) * ( GRD_xc(ijm1,k,l,AJ ,ZDIR) - GRD_x(ij    ,k0,l,ZDIR) )

          q_a(6,n,k,l) = (      cmask(6,n,k,l) ) * q_am &
                       + ( 1.D0-cmask(6,n,k,l) ) * q_ap
       enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do v = ADM_gmin_pl, ADM_gmax_pl
          q_ap = q_pl(n,k,l) + gradq_pl(n,k,l,XDIR) * ( GRD_xc_pl(v,k,l,XDIR) - GRD_x_pl(n,k0,l,XDIR) ) &
                             + gradq_pl(n,k,l,YDIR) * ( GRD_xc_pl(v,k,l,YDIR) - GRD_x_pl(n,k0,l,YDIR) ) &
                             + gradq_pl(n,k,l,ZDIR) * ( GRD_xc_pl(v,k,l,ZDIR) - GRD_x_pl(n,k0,l,ZDIR) )

          q_am = q_pl(v,k,l) + gradq_pl(v,k,l,XDIR) * ( GRD_xc_pl(v,k,l,XDIR) - GRD_x_pl(v,k0,l,XDIR) ) &
                             + gradq_pl(v,k,l,YDIR) * ( GRD_xc_pl(v,k,l,YDIR) - GRD_x_pl(v,k0,l,YDIR) ) &
                             + gradq_pl(v,k,l,ZDIR) * ( GRD_xc_pl(v,k,l,ZDIR) - GRD_x_pl(v,k0,l,ZDIR) )

          q_a_pl(v,k,l) = (      cmask_pl(v,k,l) ) * q_am &
                        + ( 1.D0-cmask_pl(v,k,l) ) * q_ap
       enddo
       enddo
       enddo
    endif

    return
  end subroutine horizontal_remap

  !-----------------------------------------------------------------------------
  subroutine vertical_limiter_thuburn( &
       q_h, q_h_pl, &
       q,   q_pl,   &
       d,   d_pl,   &
       ck,  ck_pl   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_MAX_REAL, &
       CNST_EPS_ZERO
    implicit none

    real(8), intent(inout) :: q_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: q     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: d     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: ck    (ADM_gall,   ADM_kall,ADM_lall   ,2)
    real(8), intent(in)    :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(8) :: Qin_min    (ADM_gall,   ADM_kall,2)
    real(8) :: Qin_min_pl (ADM_gall_pl,ADM_kall,2)
    real(8) :: Qin_max    (ADM_gall,   ADM_kall,2)
    real(8) :: Qin_max_pl (ADM_gall_pl,ADM_kall,2)

    real(8) :: Qout_min   (ADM_gall,   ADM_kall)
    real(8) :: Qout_min_pl(ADM_gall_pl,ADM_kall)
    real(8) :: Qout_max   (ADM_gall,   ADM_kall)
    real(8) :: Qout_max_pl(ADM_gall_pl,ADM_kall)

    real(8) :: qnext_min, qnext_min_pl
    real(8) :: qnext_max, qnext_max_pl
    real(8) :: Cin,       Cin_pl
    real(8) :: Cout,      Cout_pl
    real(8) :: CQin_min,  CQin_min_pl
    real(8) :: CQin_max,  CQin_max_pl

    real(8) :: mask, mask1, mask2
    real(8) :: zerosw

    integer :: n, k, l
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       Qin_min(:,:,:) =  CNST_MAX_REAL
       Qin_max(:,:,:) = -CNST_MAX_REAL

       do k = ADM_kmin, ADM_kmax+1
       do n = 1,ADM_gall
          if ( ck(n,k,l,1) > 0.D0 ) then ! incoming
             Qin_min(n,k-1,2) = min( q(n,k,l), q(n,k-1,l) )
             Qin_max(n,k-1,2) = max( q(n,k,l), q(n,k-1,l) )
          else ! outgoing
             Qin_min(n,k  ,1) = min( q(n,k,l), q(n,k-1,l) )
             Qin_max(n,k  ,1) = max( q(n,k,l), q(n,k-1,l) )
          endif
       enddo
       enddo

       do k = 1, ADM_kall
       do n = 1, ADM_gall

          qnext_min = min( Qin_min(n,k,1), Qin_min(n,k,2) )
          if( qnext_min ==  CNST_MAX_REAL ) qnext_min = q(n,k,l)

          qnext_max = max( Qin_max(n,k,1), Qin_max(n,k,2) )
          if( qnext_max == -CNST_MAX_REAL ) qnext_max = q(n,k,l)

          mask1 = 0.5D0 - sign(0.5D0,ck(n,k,l,1)-CNST_EPS_ZERO)
          mask2 = 0.5D0 - sign(0.5D0,ck(n,k,l,2)-CNST_EPS_ZERO)

          Cin  = (      mask1 ) * ck(n,k,l,1) &
               + (      mask2 ) * ck(n,k,l,2)
          Cout = ( 1.D0-mask1 ) * ck(n,k,l,1) &
               + ( 1.D0-mask2 ) * ck(n,k,l,2)

          CQin_max = mask1 * ( ck(n,k,l,1) * Qin_max(n,k,1) ) &
                   + mask2 * ( ck(n,k,l,2) * Qin_max(n,k,2) )
          CQin_min = mask1 * ( ck(n,k,l,1) * Qin_min(n,k,1) ) &
                   + mask2 * ( ck(n,k,l,2) * Qin_min(n,k,2) )

          zerosw = 0.5D0 - sign(0.5D0,abs(Cout)-CNST_EPS_ZERO) ! if Cout = 0, sw = 1

          Qout_min(n,k) = ( q(n,k,l) - CQin_max - qnext_max*(1.D0-Cin-Cout+d(n,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.D0 - zerosw )                        &
                        + q(n,k,l) * zerosw
          Qout_max(n,k) = ( q(n,k,l) - CQin_min - qnext_min*(1.D0-Cin-Cout+d(n,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.D0 - zerosw )                        &
                        + q(n,k,l) * zerosw

       enddo
       enddo

       do k = ADM_kmin, ADM_kmax+1
       do n = 1, ADM_gall
          mask = 0.5D0 - sign(0.5D0,ck(n,k,l,1)-CNST_EPS_ZERO)

          q_h(n,k,l) = (      mask ) * min( max( q_h(n,k,l), Qin_min(n,k  ,1) ), Qin_max(n,k  ,1) ) &
                     + ( 1.D0-mask ) * min( max( q_h(n,k,l), Qin_min(n,k-1,2) ), Qin_max(n,k-1,2) )

          q_h(n,k,l) = (      mask ) * max( min( q_h(n,k,l), Qout_max(n,k-1) ), Qout_min(n,k-1) ) &
                     + ( 1.D0-mask ) * max( min( q_h(n,k,l), Qout_max(n,k  ) ), Qout_min(n,k  ) )
       enddo
       enddo

    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          Qin_min_pl(:,:,:) =  CNST_MAX_REAL
          Qin_max_pl(:,:,:) = -CNST_MAX_REAL

          do k = ADM_kmin, ADM_kmax+1
          do n = 1, ADM_gall_pl
             if ( ck_pl(n,k,l,1) > 0.D0 ) then
                Qin_min_pl(n,k-1,2) = min( q_pl(n,k,l), q_pl(n,k-1,l) )
                Qin_max_pl(n,k-1,2) = max( q_pl(n,k,l), q_pl(n,k-1,l) )
             else
                Qin_min_pl(n,k,  1) = min( q_pl(n,k,l), q_pl(n,k-1,l) )
                Qin_max_pl(n,k,  1) = max( q_pl(n,k,l), q_pl(n,k-1,l) )
             endif
          enddo
          enddo

          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl

             qnext_min_pl = min( Qin_min_pl(n,k,1), Qin_min_pl(n,k,2) )
             if( qnext_min_pl ==  CNST_MAX_REAL ) qnext_min_pl = q_pl(n,k,l)
             qnext_min_pl = max( 0.D0, qnext_min_pl)

             qnext_max_pl = max( Qin_max_pl(n,k,1), Qin_max_pl(n,k,2) )
             if( qnext_max_pl == -CNST_MAX_REAL ) qnext_max_pl = q_pl(n,k,l)

             mask1 = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,1) )
             mask2 = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,2) )

             Cin_pl  = (      mask1 ) * ck_pl(n,k,l,1) &
                     + (      mask2 ) * ck_pl(n,k,l,2)
             Cout_pl = ( 1.D0-mask1 ) * ck_pl(n,k,l,1) &
                     + ( 1.D0-mask2 ) * ck_pl(n,k,l,2)

             CQin_max_pl = mask1 * ( ck_pl(n,k,l,1) * Qin_max_pl(n,k,1) ) &
                         + mask2 * ( ck_pl(n,k,l,2) * Qin_max_pl(n,k,2) )
             CQin_min_pl = mask1 * ( ck_pl(n,k,l,1) * Qin_min_pl(n,k,1) ) &
                         + mask2 * ( ck_pl(n,k,l,2) * Qin_min_pl(n,k,2) )

             zerosw = 0.5D0 - sign(0.5D0,abs(Cout_pl)-CNST_EPS_ZERO) ! if Cout_pl = 0, sw = 1

             Qout_min_pl(n,k) = ( q_pl(n,k,l) - CQin_max_pl - qnext_max_pl*(1.D0-Cin_pl-Cout_pl+d_pl(n,k,l)) ) &
                              / ( Cout_pl + zerosw ) * ( 1.D0 - zerosw )                                       &
                              + q_pl(n,k,l) * zerosw
             Qout_max_pl(n,k) = ( q_pl(n,k,l) - CQin_min_pl - qnext_min_pl*(1.D0-Cin_pl-Cout_pl+d_pl(n,k,l)) ) &
                              / ( Cout_pl + zerosw ) * ( 1.D0 - zerosw )                                       &
                              + q_pl(n,k,l) * zerosw
          enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
          do n = 1, ADM_gall_pl
             mask = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,1) )

             q_h_pl(n,k,l) = (      mask ) * min( max( q_h_pl(n,k,l), Qin_min_pl(n,k  ,1) ), Qin_max_pl(n,k  ,1) ) &
                           + ( 1.D0-mask ) * min( max( q_h_pl(n,k,l), Qin_min_pl(n,k-1,2) ), Qin_max_pl(n,k-1,2) )

             q_h_pl(n,k,l) = (      mask ) * max( min( q_h_pl(n,k,l), Qout_max_pl(n,k-1) ), Qout_min_pl(n,k-1) ) &
                           + ( 1.D0-mask ) * max( min( q_h_pl(n,k,l), Qout_max_pl(n,k  ) ), Qout_min_pl(n,k  ) )
          enddo
          enddo

       enddo
    endif

    return
  end subroutine vertical_limiter_thuburn

  !-----------------------------------------------------------------------------
  !> Miura(2004)'s scheme with Thuburn(1996) limiter
  subroutine horizontal_limiter_thuburn( &
       q_a,    q_a_pl,  &
       q,      q_pl,    &
       d,      d_pl,    &
       ch,     ch_pl,   &
       cmask,  cmask_pl )
    use mod_adm, only: &
       ADM_have_pl,    &
       ADM_have_sgp,   &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl
    use mod_cnst, only: &
       CNST_MAX_REAL, &
       CNST_EPS_ZERO
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(8), intent(inout) :: q_a     (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: q_a_pl  (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: q       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: q_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: d       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: d_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: ch      (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: ch_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: cmask   (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: cmask_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: q_min_AI, q_min_AIJ, q_min_AJ, q_min_pl
    real(8) :: q_max_AI, q_max_AIJ, q_max_AJ, q_max_pl

    real(8) :: qnext_min   (ADM_gall), qnext_min_pl
    real(8) :: qnext_max   (ADM_gall), qnext_max_pl
    real(8) :: Cin_sum     (ADM_gall), Cin_sum_pl
    real(8) :: Cout_sum    (ADM_gall), Cout_sum_pl
    real(8) :: CQin_max_sum(ADM_gall), CQin_max_sum_pl
    real(8) :: CQin_min_sum(ADM_gall), CQin_min_sum_pl

    integer, parameter :: I_min = 1
    integer, parameter :: I_max = 2
    real(8) :: Qin    (6,ADM_gall   ,ADM_kall,ADM_lall   ,2)
    real(8) :: Qin_pl (2,ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(8) :: Qout   (  ADM_gall   ,ADM_kall,ADM_lall   ,2)
    real(8) :: Qout_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(8) :: zerosw

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1, ip2jp1
    integer :: im1j, ijm1

    integer :: nstart, nend
    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    !---< (i) define inflow bounds, eq.(32)&(33) >---
    do l = 1, ADM_lall
    do k = 1, ADM_kall
       nstart = suf(ADM_gmin,ADM_gmin)
       nend   = suf(ADM_gmax,ADM_gmax)

       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d

          q_min_AI  = min( q(ij,k,l), q(ijm1,k,l), q(ip1j,k,l), q(ip1jp1,k,l) )
          q_max_AI  = max( q(ij,k,l), q(ijm1,k,l), q(ip1j,k,l), q(ip1jp1,k,l) )
          q_min_AIJ = min( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijp1,k,l) )
          q_max_AIJ = max( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijp1,k,l) )
          q_min_AJ  = min( q(ij,k,l), q(ip1jp1,k,l), q(ijp1,k,l), q(im1j,k,l) )
          q_max_AJ  = max( q(ij,k,l), q(ip1jp1,k,l), q(ijp1,k,l), q(im1j,k,l) )

          Qin(1,ij,    k,l,I_min) = (      cmask(1,n,k,l) ) * q_min_AI         &
                                  + ( 1.D0-cmask(1,n,k,l) ) * CNST_MAX_REAL
          Qin(4,ip1j,  k,l,I_min) = (      cmask(1,n,k,l) ) * CNST_MAX_REAL    &
                                  + ( 1.D0-cmask(1,n,k,l) ) * q_min_AI
          Qin(1,ij,    k,l,I_max) = (      cmask(1,n,k,l) ) * q_max_AI         &
                                  + ( 1.D0-cmask(1,n,k,l) ) * (-CNST_MAX_REAL)
          Qin(4,ip1j,  k,l,I_max) = (      cmask(1,n,k,l) ) * (-CNST_MAX_REAL) &
                                  + ( 1.D0-cmask(1,n,k,l) ) * q_max_AI

          Qin(2,ij,    k,l,I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
                                  + ( 1.D0-cmask(2,n,k,l) ) * CNST_MAX_REAL
          Qin(5,ip1jp1,k,l,I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
                                  + ( 1.D0-cmask(2,n,k,l) ) * q_min_AIJ
          Qin(2,ij,    k,l,I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
                                  + ( 1.D0-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
          Qin(5,ip1jp1,k,l,I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
                                  + ( 1.D0-cmask(2,n,k,l) ) * q_max_AIJ

          Qin(3,ij,    k,l,I_min) = (      cmask(3,n,k,l) ) * q_min_AJ         &
                                  + ( 1.D0-cmask(3,n,k,l) ) * CNST_MAX_REAL
          Qin(6,ijp1,  k,l,I_min) = (      cmask(3,n,k,l) ) * CNST_MAX_REAL    &
                                  + ( 1.D0-cmask(3,n,k,l) ) * q_min_AJ
          Qin(3,ij,    k,l,I_max) = (      cmask(3,n,k,l) ) * q_max_AJ         &
                                  + ( 1.D0-cmask(3,n,k,l) ) * (-CNST_MAX_REAL)
          Qin(6,ijp1,  k,l,I_max) = (      cmask(3,n,k,l) ) * (-CNST_MAX_REAL) &
                                  + ( 1.D0-cmask(3,n,k,l) ) * q_max_AJ
       enddo

       ! peeling
       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmin-1,ADM_gmin  )

       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          q_min_AI  = min( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijm1,k,l) )
          q_max_AI  = max( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijm1,k,l) )

          Qin(1,ij,    k,l,I_min) = (      cmask(1,n,k,l) ) * q_min_AI         &
                                  + ( 1.D0-cmask(1,n,k,l) ) * CNST_MAX_REAL
          Qin(4,ip1j,  k,l,I_min) = (      cmask(1,n,k,l) ) * CNST_MAX_REAL    &
                                  + ( 1.D0-cmask(1,n,k,l) ) * q_min_AI
          Qin(1,ij,    k,l,I_max) = (      cmask(1,n,k,l) ) * q_max_AI         &
                                  + ( 1.D0-cmask(1,n,k,l) ) * (-CNST_MAX_REAL)
          Qin(4,ip1j,  k,l,I_max) = (      cmask(1,n,k,l) ) * (-CNST_MAX_REAL) &
                                  + ( 1.D0-cmask(1,n,k,l) ) * q_max_AI
       enddo

       ! peeling
       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmin-1,ADM_gmin  )

       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          q_min_AIJ = min( q(ij,k,l), q(ip1jp1,k,l), q(ip1j,k,l), q(ijp1,k,l) )
          q_max_AIJ = max( q(ij,k,l), q(ip1jp1,k,l), q(ip1j,k,l), q(ijp1,k,l) )

          Qin(2,ij,    k,l,I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
                                  + ( 1.D0-cmask(2,n,k,l) ) * CNST_MAX_REAL
          Qin(5,ip1jp1,k,l,I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
                                  + ( 1.D0-cmask(2,n,k,l) ) * q_min_AIJ
          Qin(2,ij,    k,l,I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
                                  + ( 1.D0-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
          Qin(5,ip1jp1,k,l,I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
                                  + ( 1.D0-cmask(2,n,k,l) ) * q_max_AIJ
       enddo

       ! peeling
       nstart = suf(ADM_gmin,  ADM_gmin-1)
       nend   = suf(ADM_gmin-1,ADM_gmin  )

       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1

          q_min_AJ  = min( q(ij,k,l), q(ijp1,k,l), q(ip1jp1,k,l), q(im1j,k,l) )
          q_max_AJ  = max( q(ij,k,l), q(ijp1,k,l), q(ip1jp1,k,l), q(im1j,k,l) )

          Qin(3,ij,    k,l,I_min) = (      cmask(3,n,k,l) ) * q_min_AJ         &
                                  + ( 1.D0-cmask(3,n,k,l) ) * CNST_MAX_REAL
          Qin(6,ijp1,  k,l,I_min) = (      cmask(3,n,k,l) ) * CNST_MAX_REAL    &
                                  + ( 1.D0-cmask(3,n,k,l) ) * q_min_AJ
          Qin(3,ij,    k,l,I_max) = (      cmask(3,n,k,l) ) * q_max_AJ         &
                                  + ( 1.D0-cmask(3,n,k,l) ) * (-CNST_MAX_REAL)
          Qin(6,ijp1,  k,l,I_max) = (      cmask(3,n,k,l) ) * (-CNST_MAX_REAL) &
                                  + ( 1.D0-cmask(3,n,k,l) ) * q_max_AJ
       enddo

       if ( ADM_have_sgp(l) ) then
          n = suf(ADM_gmin-1,ADM_gmin-1)
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          ip2jp1 = n + 2 + ADM_gall_1d

          q_min_AIJ = min( q(ij,k,l), q(ip1jp1,k,l), q(ip2jp1,k,l), q(ijp1,k,l) )
          q_max_AIJ = max( q(ij,k,l), q(ip1jp1,k,l), q(ip2jp1,k,l), q(ijp1,k,l) )

          Qin(2,ij,    k,l,I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
                                  + ( 1.D0-cmask(2,n,k,l) ) * CNST_MAX_REAL
          Qin(5,ip1jp1,k,l,I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
                                  + ( 1.D0-cmask(2,n,k,l) ) * q_min_AIJ
          Qin(2,ij,    k,l,I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
                                  + ( 1.D0-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
          Qin(5,ip1jp1,k,l,I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
                                  + ( 1.D0-cmask(2,n,k,l) ) * q_max_AIJ
       endif

    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do v = ADM_gmin_pl, ADM_gmax_pl
          ij   = v
          ijp1 = v + 1
          ijm1 = v - 1
          if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl
          if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl

          q_min_pl = min( q_pl(n,k,l), q_pl(ij,k,l), q_pl(ijm1,k,l), q_pl(ijp1,k,l) )
          q_max_pl = max( q_pl(n,k,l), q_pl(ij,k,l), q_pl(ijm1,k,l), q_pl(ijp1,k,l) )

          Qin_pl(1,ij,k,l,I_min) = (      cmask_pl(ij,k,l) ) * q_min_pl         &
                                 + ( 1.D0-cmask_pl(ij,k,l) ) * CNST_MAX_REAL
          Qin_pl(2,ij,k,l,I_min) = (      cmask_pl(ij,k,l) ) * CNST_MAX_REAL    &
                                 + ( 1.D0-cmask_pl(ij,k,l) ) * q_min_pl
          Qin_pl(1,ij,k,l,I_max) = (      cmask_pl(ij,k,l) ) * q_max_pl         &
                                 + ( 1.D0-cmask_pl(ij,k,l) ) * (-CNST_MAX_REAL)
          Qin_pl(2,ij,k,l,I_max) = (      cmask_pl(ij,k,l) ) * (-CNST_MAX_REAL) &
                                 + ( 1.D0-cmask_pl(ij,k,l) ) * q_max_pl
       enddo
       enddo
       enddo
    endif

    !---< (iii) define allowable range of q at next step, eq.(42)&(43) >---
    nstart = suf(ADM_gmin, ADM_gmin )
    nend   = suf(ADM_gmax ,ADM_gmax )

    do l = 1, ADM_lall
    do k = 1, ADM_kall

       do n = nstart, nend
          qnext_min(n) = minval( Qin(1:6,n,k,l,I_min) )
          if( qnext_min(n) ==  CNST_MAX_REAL ) qnext_min(n) = q(n,k,l)
       enddo

       do n = nstart, nend
          qnext_max(n) = maxval( Qin(1:6,n,k,l,I_max) )
          if( qnext_max(n) == -CNST_MAX_REAL ) qnext_max(n) = q(n,k,l)
       enddo

       do n = nstart, nend
          Cin_sum(n)  = (      cmask(1,n,k,l) ) * ch(1,n,k,l) &
                      + (      cmask(2,n,k,l) ) * ch(2,n,k,l) &
                      + (      cmask(3,n,k,l) ) * ch(3,n,k,l) &
                      + (      cmask(4,n,k,l) ) * ch(4,n,k,l) &
                      + (      cmask(5,n,k,l) ) * ch(5,n,k,l) &
                      + (      cmask(6,n,k,l) ) * ch(6,n,k,l)

          Cout_sum(n) = ( 1.D0-cmask(1,n,k,l) ) * ch(1,n,k,l) &
                      + ( 1.D0-cmask(2,n,k,l) ) * ch(2,n,k,l) &
                      + ( 1.D0-cmask(3,n,k,l) ) * ch(3,n,k,l) &
                      + ( 1.D0-cmask(4,n,k,l) ) * ch(4,n,k,l) &
                      + ( 1.D0-cmask(5,n,k,l) ) * ch(5,n,k,l) &
                      + ( 1.D0-cmask(6,n,k,l) ) * ch(6,n,k,l)
       enddo

       do n = nstart, nend
          CQin_min_sum(n) = cmask(1,n,k,l) * ch(1,n,k,l) * Qin(1,n,k,l,I_min) &
                          + cmask(2,n,k,l) * ch(2,n,k,l) * Qin(2,n,k,l,I_min) &
                          + cmask(3,n,k,l) * ch(3,n,k,l) * Qin(3,n,k,l,I_min) &
                          + cmask(4,n,k,l) * ch(4,n,k,l) * Qin(4,n,k,l,I_min) &
                          + cmask(5,n,k,l) * ch(5,n,k,l) * Qin(5,n,k,l,I_min) &
                          + cmask(6,n,k,l) * ch(6,n,k,l) * Qin(6,n,k,l,I_min)

          CQin_max_sum(n) = cmask(1,n,k,l) * ch(1,n,k,l) * Qin(1,n,k,l,I_max) &
                          + cmask(2,n,k,l) * ch(2,n,k,l) * Qin(2,n,k,l,I_max) &
                          + cmask(3,n,k,l) * ch(3,n,k,l) * Qin(3,n,k,l,I_max) &
                          + cmask(4,n,k,l) * ch(4,n,k,l) * Qin(4,n,k,l,I_max) &
                          + cmask(5,n,k,l) * ch(5,n,k,l) * Qin(5,n,k,l,I_max) &
                          + cmask(6,n,k,l) * ch(6,n,k,l) * Qin(6,n,k,l,I_max)
       enddo

       do n = nstart, nend
          zerosw = 0.5D0 - sign(0.5D0,abs(Cout_sum(n))-CNST_EPS_ZERO) ! if Cout_sum = 0, sw = 1

          Qout(n,k,l,I_min) = ( q(n,k,l) - CQin_max_sum(n) - qnext_max(n)*(1.D0-Cin_sum(n)-Cout_sum(n)+d(n,k,l)) ) &
                            / ( Cout_sum(n) + zerosw ) * ( 1.D0 - zerosw )                                         &
                            + q(n,k,l) * zerosw
          Qout(n,k,l,I_max) = ( q(n,k,l) - CQin_min_sum(n) - qnext_min(n)*(1.D0-Cin_sum(n)-Cout_sum(n)+d(n,k,l)) ) &
                            / ( Cout_sum(n) + zerosw ) * ( 1.D0 - zerosw )                                         &
                            + q(n,k,l) * zerosw
       enddo ! n loop
    enddo ! k loop
    enddo ! l loop

    Qout(     1:nstart-1,:,:,I_min) = q(     1:nstart-1,:,:)
    Qout(nend+1:ADM_gall,:,:,I_min) = q(nend+1:ADM_gall,:,:)
    Qout(     1:nstart-1,:,:,I_max) = q(     1:nstart-1,:,:)
    Qout(nend+1:ADM_gall,:,:,I_max) = q(nend+1:ADM_gall,:,:)

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          qnext_min_pl = minval( Qin_pl(1,ADM_gmin_pl:ADM_gmax_pl,k,l,I_min) )
          if( qnext_min_pl ==  CNST_MAX_REAL ) qnext_min_pl = q_pl(n,k,l)

          qnext_max_pl = maxval( Qin_pl(1,ADM_gmin_pl:ADM_gmax_pl,k,l,I_max) )
          if( qnext_max_pl == -CNST_MAX_REAL ) qnext_max_pl = q_pl(n,k,l)

          Cin_sum_pl      = 0.D0
          Cout_sum_pl     = 0.D0
          CQin_max_sum_pl = 0.D0
          CQin_min_sum_pl = 0.D0
          do v = ADM_gmin_pl, ADM_gmax_pl
             Cin_sum_pl      = Cin_sum_pl      + (      cmask_pl(v,k,l) ) * ch_pl(v,k,l)
             Cout_sum_pl     = Cout_sum_pl     + ( 1.D0-cmask_pl(v,k,l) ) * ch_pl(v,k,l)
             CQin_min_sum_pl = CQin_min_sum_pl + (      cmask_pl(v,k,l) ) * ch_pl(v,k,l) * Qin_pl(1,v,k,l,I_min)
             CQin_max_sum_pl = CQin_max_sum_pl + (      cmask_pl(v,k,l) ) * ch_pl(v,k,l) * Qin_pl(1,v,k,l,I_max)
          enddo

          zerosw = 0.5D0 - sign(0.5D0,abs(Cout_sum_pl)-CNST_EPS_ZERO) ! if Cout_sum_pl = 0, sw = 1

          Qout_pl(n,k,l,I_min) = ( q_pl(n,k,l) - CQin_max_sum_pl - qnext_max_pl*(1.D0-Cin_sum_pl-Cout_sum_pl+d_pl(n,k,l)) ) &
                               / ( Cout_sum_pl + zerosw ) * ( 1.D0 - zerosw )                                               &
                               + q_pl(n,k,l) * zerosw
          Qout_pl(n,k,l,I_max) = ( q_pl(n,k,l) - CQin_min_sum_pl - qnext_min_pl*(1.D0-Cin_sum_pl-Cout_sum_pl+d_pl(n,k,l)) ) &
                               / ( Cout_sum_pl + zerosw ) * ( 1.D0 - zerosw )                                               &
                               + q_pl(n,k,l) * zerosw
       enddo
       enddo
    endif

    call COMM_data_transfer( Qout(:,:,:,:), Qout_pl(:,:,:,:) )

    !---- apply inflow/outflow limiter
    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax  ,ADM_gmax  )

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          q_a(1,n,k,l) = (      cmask(1,n,k,l) ) * min(max(q_a(1,n,k,l), Qin (1,ij    ,k,l,I_min)), Qin (1,ij    ,k,l,I_max)) &
                       + ( 1.D0-cmask(1,n,k,l) ) * min(max(q_a(1,n,k,l), Qin (4,ip1j  ,k,l,I_min)), Qin (4,ip1j  ,k,l,I_max))
          q_a(1,n,k,l) = (      cmask(1,n,k,l) ) * max(min(q_a(1,n,k,l), Qout(  ip1j  ,k,l,I_max)), Qout(  ip1j  ,k,l,I_min)) &
                       + ( 1.D0-cmask(1,n,k,l) ) * max(min(q_a(1,n,k,l), Qout(  ij    ,k,l,I_max)), Qout(  ij    ,k,l,I_min))
          q_a(4,ip1j,k,l) = q_a(1,n,k,l)

          q_a(2,n,k,l) = (      cmask(2,n,k,l) ) * min(max(q_a(2,n,k,l), Qin (2,ij    ,k,l,I_min)), Qin (2,ij    ,k,l,I_max)) &
                       + ( 1.D0-cmask(2,n,k,l) ) * min(max(q_a(2,n,k,l), Qin (5,ip1jp1,k,l,I_min)), Qin (5,ip1jp1,k,l,I_max))
          q_a(2,n,k,l) = (      cmask(2,n,k,l) ) * max(min(q_a(2,n,k,l), Qout(  ip1jp1,k,l,I_max)), Qout(  ip1jp1,k,l,I_min)) &
                       + ( 1.D0-cmask(2,n,k,l) ) * max(min(q_a(2,n,k,l), Qout(  ij    ,k,l,I_max)), Qout(  ij    ,k,l,I_min))
          q_a(5,ip1jp1,k,l) = q_a(2,n,k,l)

          q_a(3,n,k,l) = (      cmask(3,n,k,l) ) * min(max(q_a(3,n,k,l), Qin (3,ij    ,k,l,I_min)), Qin (3,ij    ,k,l,I_max)) &
                       + ( 1.D0-cmask(3,n,k,l) ) * min(max(q_a(3,n,k,l), Qin (6,ijp1  ,k,l,I_min)), Qin (6,ijp1  ,k,l,I_max))
          q_a(3,n,k,l) = (      cmask(3,n,k,l) ) * max(min(q_a(3,n,k,l), Qout(  ijp1  ,k,l,I_max)), Qout(  ijp1  ,k,l,I_min)) &
                       + ( 1.D0-cmask(3,n,k,l) ) * max(min(q_a(3,n,k,l), Qout(  ij    ,k,l,I_max)), Qout(  ij    ,k,l,I_min))
          q_a(6,ijp1,k,l) = q_a(3,n,k,l)
       enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do v = ADM_gmin_pl, ADM_gmax_pl
          q_a_pl(v,k,l) = (      cmask_pl(v,k,l) ) * min(max(q_a_pl(v,k,l), Qin_pl (1,v,k,l,I_min)), Qin_pl (1,v,k,l,I_max)) &
                        + ( 1.D0-cmask_pl(v,k,l) ) * min(max(q_a_pl(v,k,l), Qin_pl (2,v,k,l,I_min)), Qin_pl (2,v,k,l,I_max))
          q_a_pl(v,k,l) = (      cmask_pl(v,k,l) ) * max(min(q_a_pl(v,k,l), Qout_pl(  v,k,l,I_max)), Qout_pl(  v,k,l,I_min)) &
                        + ( 1.D0-cmask_pl(v,k,l) ) * max(min(q_a_pl(v,k,l), Qout_pl(  n,k,l,I_max)), Qout_pl(  n,k,l,I_min))
       enddo
       enddo
       enddo
    endif

    return
  end subroutine horizontal_limiter_thuburn

end module mod_src_tracer
!-------------------------------------------------------------------------------------
