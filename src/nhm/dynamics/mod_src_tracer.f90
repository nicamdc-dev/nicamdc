!-------------------------------------------------------------------------------
!> Module tracer advection
!!
!! @par Description
!!         This module contains subroutines for tracer advection
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_src_tracer
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof

  use mod_adm, only: &
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
     P_RAREA => GMTR_p_RAREA, &
     T_RAREA => GMTR_t_RAREA, &
     W1      => GMTR_t_W1,    &
     W2      => GMTR_t_W2,    &
     W3      => GMTR_t_W3,    &
     HNX     => GMTR_A_HNX,   &
     HNY     => GMTR_A_HNY,   &
     HNZ     => GMTR_A_HNZ,   &
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

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: horizontal_flux
  private :: horizontal_remap
  private :: vertical_limiter_thuburn
  private :: horizontal_limiter_thuburn

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------------
  subroutine src_tracer_advection( &
       vmax,                        &
       rhogq,       rhogq_pl,       &
       rhog_in,     rhog_in_pl,     &
       rhog_mean,   rhog_mean_pl,   &
       rhogvx_mean, rhogvx_mean_pl, &
       rhogvy_mean, rhogvy_mean_pl, &
       rhogvz_mean, rhogvz_mean_pl, &
       rhogw_mean,  rhogw_mean_pl,  &
       frhog,       frhog_pl,       &
       dt,                          &
       thubern_lim                  )
    use mod_const, only: &
       CONST_EPS
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
       ADM_gmin_pl, &
       ADM_gmax_pl, &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgz, &
       GRD_afac, &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_RGSQRTH,      &
       VMTR_RGSQRTH_pl,   &
       VMTR_RGAM,         &
       VMTR_RGAM_pl,      &
       VMTR_RGAMH,        &
       VMTR_RGAMH_pl,     &
       VMTR_C2WfactGz,    &
       VMTR_C2WfactGz_pl
    implicit none

    integer,  intent(in)    :: vmax                                                  ! number of tracers
    real(RP), intent(inout) :: rhogq         (ADM_gall   ,ADM_kall,ADM_lall   ,vmax) ! rhogq   ( G^1/2 x gam2 )
    real(RP), intent(inout) :: rhogq_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax)
    real(RP), intent(in)    :: rhog_in       (ADM_gall   ,ADM_kall,ADM_lall   )      ! rho(old)( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog_in_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog_mean     (ADM_gall   ,ADM_kall,ADM_lall   )      ! rho     ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog_mean_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhogvx_mean   (ADM_gall   ,ADM_kall,ADM_lall   )      ! rho*Vx  ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogvx_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhogvy_mean   (ADM_gall   ,ADM_kall,ADM_lall   )      ! rho*Vy  ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogvy_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhogvz_mean   (ADM_gall   ,ADM_kall,ADM_lall   )      ! rho*Vz  ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogvz_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhogw_mean    (ADM_gall   ,ADM_kall,ADM_lall   )      ! rho*w   ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogw_mean_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: frhog         (ADM_gall   ,ADM_kall,ADM_lall   )      ! hyperviscosity tendency for rhog
    real(RP), intent(in)    :: frhog_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: dt                                                    ! delta t
    logical,  intent(in)    :: thubern_lim                                           ! switch of thubern limiter [add] 20130613 R.Yoshida

    real(RP) :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: q        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: d        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: d_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: q_h      (ADM_gall   ,ADM_kall,ADM_lall   )                 ! q at layer face
    real(RP) :: q_h_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: flx_v    (ADM_gall   ,ADM_kall,ADM_lall   )                 ! mass flux
    real(RP) :: flx_v_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ck       (ADM_gall   ,ADM_kall,ADM_lall   ,2)               ! Courant number
    real(RP) :: ck_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(RP) :: q_a      (ADM_gall   ,ADM_kall,ADM_lall   ,6)               ! q at cell face
    real(RP) :: q_a_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP) :: flx_h    (ADM_gall   ,ADM_kall,ADM_lall   ,6)               ! mass flux
    real(RP) :: flx_h_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP) :: ch       (ADM_gall   ,ADM_kall,ADM_lall   ,6)               ! Courant number
    real(RP) :: ch_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP) :: cmask    (ADM_gall   ,ADM_kall,ADM_lall   ,6)               ! upwind direction mask
    real(RP) :: cmask_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP) :: GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,AI:AJ,XDIR:ZDIR) ! mass centroid position
    real(RP) :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)

    real(RP), parameter :: b1 = 0.0_RP
    real(RP), parameter :: b2 = 1.0_RP
    real(RP), parameter :: b3 = 1.0_RP - (b1+b2)

    integer  :: gall, kall, kmin, kmax
    real(RP) :: EPS

    integer  :: nstart, nend
    integer  :: g, k, l, v, iq
    !---------------------------------------------------------------------------

    gall = ADM_gall
    kall = ADM_kall
    kmin = ADM_kmin
    kmax = ADM_kmax

    EPS = CONST_EPS

    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 1st
    !---------------------------------------------------------------------------
    call PROF_rapstart('____vertical_adv',2)

    do l = 1, ADM_lall
       !$omp parallel default(none),private(g,k), &
       !$omp shared(l,gall,kall,kmin,kmax,flx_v,rhogvx_mean,rhogvy_mean,rhogvz_mean,rhogw_mean,dt, &
       !$omp        d,ck,frhog,rhog_in,GRD_rdgz,VMTR_RGAMH,VMTR_RGSQRTH,VMTR_C2WfactGz)

       !$omp do
       do k = kmin+1, kmax
       do g = 1, gall
          flx_v(g,k,l) = ( ( VMTR_C2WfactGz(g,k,1,l) * rhogvx_mean(g,k  ,l) &
                           + VMTR_C2WfactGz(g,k,2,l) * rhogvx_mean(g,k-1,l) &
                           + VMTR_C2WfactGz(g,k,3,l) * rhogvy_mean(g,k  ,l) &
                           + VMTR_C2WfactGz(g,k,4,l) * rhogvy_mean(g,k-1,l) &
                           + VMTR_C2WfactGz(g,k,5,l) * rhogvz_mean(g,k  ,l) &
                           + VMTR_C2WfactGz(g,k,6,l) * rhogvz_mean(g,k-1,l) &
                           ) * VMTR_RGAMH(g,k,l)                            & ! horizontal contribution
                         + rhogw_mean(g,k,l) * VMTR_RGSQRTH(g,k,l)          & ! vertical   contribution
                         ) * 0.5_RP * dt
       enddo
       enddo
       !$omp end do nowait

       !$omp do
       do g = 1, gall
          flx_v(g,kmin,  l) = 0.0_RP
          flx_v(g,kmax+1,l) = 0.0_RP
       enddo
       !$omp end do

       !$omp do
       do k = 1, kall
       do g = 1, gall
          d(g,k,l) = b1 * frhog(g,k,l) / rhog_in(g,k,l) * dt
       enddo
       enddo
       !$omp end do nowait

       !$omp do
       do k = kmin, kmax
       do g = 1, gall
          ck(g,k,l,1) = -flx_v(g,k  ,l) / rhog_in(g,k,l) * GRD_rdgz(k)
          ck(g,k,l,2) =  flx_v(g,k+1,l) / rhog_in(g,k,l) * GRD_rdgz(k)
       enddo
       enddo
       !$omp end do nowait

       !$omp do
       do g = 1, gall
          ck(g,kmin-1,l,1) = 0.0_RP
          ck(g,kmin-1,l,2) = 0.0_RP
          ck(g,kmax+1,l,1) = 0.0_RP
          ck(g,kmax+1,l,2) = 0.0_RP
       enddo
       !$omp end do

       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             flx_v_pl(g,k,l) = ( ( VMTR_C2WfactGz_pl(g,k,1,l) * rhogvx_mean_pl(g,k  ,l) &
                                 + VMTR_C2WfactGz_pl(g,k,2,l) * rhogvx_mean_pl(g,k-1,l) &
                                 + VMTR_C2WfactGz_pl(g,k,3,l) * rhogvy_mean_pl(g,k  ,l) &
                                 + VMTR_C2WfactGz_pl(g,k,4,l) * rhogvy_mean_pl(g,k-1,l) &
                                 + VMTR_C2WfactGz_pl(g,k,5,l) * rhogvz_mean_pl(g,k  ,l) &
                                 + VMTR_C2WfactGz_pl(g,k,6,l) * rhogvz_mean_pl(g,k-1,l) &
                                 ) * VMTR_RGAMH_pl(g,k,l)                               & ! horizontal contribution
                               + rhogw_mean_pl(g,k,l) * VMTR_RGSQRTH_pl(g,k,l)          & ! vertical   contribution
                               ) * 0.5_RP * dt
          enddo
          enddo
          do g = 1, ADM_gall_pl
             flx_v_pl(g,ADM_kmin,  l) = 0.0_RP
             flx_v_pl(g,ADM_kmax+1,l) = 0.0_RP
          enddo

          d_pl(:,:,l) = b1 * frhog_pl(:,:,l) / rhog_in_pl(:,:,l) * dt

          do k = ADM_kmin, ADM_kmax
             ck_pl(:,k,l,1) = -flx_v_pl(:,k  ,l) / rhog_in_pl(:,k,l) * GRD_rdgz(k)
             ck_pl(:,k,l,2) =  flx_v_pl(:,k+1,l) / rhog_in_pl(:,k,l) * GRD_rdgz(k)
          enddo
          ck_pl(:,ADM_kmin-1,l,1) = 0.0_RP
          ck_pl(:,ADM_kmin-1,l,2) = 0.0_RP
          ck_pl(:,ADM_kmax+1,l,1) = 0.0_RP
          ck_pl(:,ADM_kmax+1,l,2) = 0.0_RP
       enddo
    endif

    !--- vertical advection: 2nd-order centered difference
    do iq = 1, vmax

       do l = 1, ADM_lall
          !$omp parallel default(none),private(g,k), &
          !$omp shared(l,iq,gall,kall,kmin,kmax,q_h,q,rhogq,rhog_in,GRD_afac,GRD_bfac)

          !$omp do
          do k = 1, kall
          do g = 1, gall
             q(g,k,l) = rhogq(g,k,l,iq) / rhog_in(g,k,l)
          enddo
          enddo
          !$omp end do

          !$omp do
          do k = kmin, kmax+1
          do g = 1, gall
             q_h(g,k,l) = GRD_afac(k) * q(g,k,  l) &
                        + GRD_bfac(k) * q(g,k-1,l)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do g = 1, gall
             q_h(g,kmin-1,l) = 0.0_RP
          enddo
          !$omp end do

          !$omp end parallel
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             q_pl(:,:,l) = rhogq_pl(:,:,l,iq) / rhog_in_pl(:,:,l)

             do k = ADM_kmin, ADM_kmax+1
                q_h_pl(:,k,l) = GRD_afac(k) * q_pl(:,k,  l) &
                              + GRD_bfac(k) * q_pl(:,k-1,l)
             enddo
             q_h_pl(:,ADM_kmin-1,l) = 0.0_RP
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
          !$omp parallel default(none),private(g,k), &
          !$omp shared(l,iq,gall,kmin,kmax,rhogq,q_h,flx_v,GRD_rdgz)

          !$omp do
          do g = 1, gall
             q_h(g,kmin  ,l) = 0.0_RP
             q_h(g,kmax+1,l) = 0.0_RP
          enddo
          !$omp end do

          !$omp do
          do k = kmin, kmax
          do g = 1, gall
             rhogq(g,k,l,iq) = rhogq(g,k,l,iq) - ( flx_v(g,k+1,l) * q_h(g,k+1,l) &
                                                 - flx_v(g,k,  l) * q_h(g,k,  l) &
                                                 ) * GRD_rdgz(k)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do g = 1, gall
             rhogq(g,kmin-1,l,iq) = 0.0_RP
             rhogq(g,kmax+1,l,iq) = 0.0_RP
          enddo
          !$omp end do

          !$omp end parallel
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             q_h_pl(:,ADM_kmin  ,l) = 0.0_RP
             q_h_pl(:,ADM_kmax+1,l) = 0.0_RP

             do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                rhogq_pl(g,k,l,iq) = rhogq_pl(g,k,l,iq) - ( flx_v_pl(g,k+1,l) * q_h_pl(g,k+1,l) &
                                                          - flx_v_pl(g,k,  l) * q_h_pl(g,k,  l) &
                                                          ) * GRD_rdgz(k)
             enddo
             enddo

             do g = 1, ADM_gall_pl
                rhogq_pl(g,ADM_kmin-1,l,iq) = 0.0_RP
                rhogq_pl(g,ADM_kmax+1,l,iq) = 0.0_RP
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

    !--- update rhog
    do l = 1, ADM_lall
       !$omp parallel default(none),private(g,k), &
       !$omp shared(l,gall,kmin,kmax,rhog,rhog_in,flx_v,GRD_rdgz,frhog,dt)

       !$omp do
       do k = kmin, kmax
       do g = 1, gall
          rhog(g,k,l) = rhog_in(g,k,l) - ( flx_v(g,k+1,l) &
                                         - flx_v(g,k  ,l) &
                                         ) * GRD_rdgz(k)  &
                                       + b1 * frhog(g,k,l) * dt
       enddo
       enddo
       !$omp end do

       !$omp do
       do g = 1, gall
          rhog(g,kmin-1,l) = rhog_in(g,kmin,l)
          rhog(g,kmax+1,l) = rhog_in(g,kmax,l)
       enddo
       !$omp end do

       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             rhog_pl(g,k,l) = rhog_in_pl(g,k,l) - ( flx_v_pl(g,k+1,l) &
                                                  - flx_v_pl(g,k  ,l) &
                                                  ) * GRD_rdgz(k)     &
                                                + b1 * frhog_pl(g,k,l) * dt
          enddo
          enddo

          rhog_pl(:,ADM_kmin-1,l) = rhog_in_pl(:,ADM_kmin,l)
          rhog_pl(:,ADM_kmax+1,l) = rhog_in_pl(:,ADM_kmax,l)
       enddo
    endif

    call PROF_rapend('____vertical_adv',2)
    !---------------------------------------------------------------------------
    ! Horizontal advection by MIURA scheme
    !---------------------------------------------------------------------------
    call PROF_rapstart('____horizontal_adv',2)

    d(:,:,:) = b2 * frhog(:,:,:) / rhog(:,:,:) * dt

    !$omp parallel workshare
    rhogvx(:,:,:) = rhogvx_mean(:,:,:) * VMTR_RGAM(:,:,:)
    rhogvy(:,:,:) = rhogvy_mean(:,:,:) * VMTR_RGAM(:,:,:)
    rhogvz(:,:,:) = rhogvz_mean(:,:,:) * VMTR_RGAM(:,:,:)
    !$omp end parallel workshare

    if ( ADM_have_pl ) then
       d_pl(:,:,:) = b2 * frhog_pl(:,:,:) / rhog_pl(:,:,:) * dt

       rhogvx_pl(:,:,:) = rhogvx_mean_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
       rhogvy_pl(:,:,:) = rhogvy_mean_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
       rhogvz_pl(:,:,:) = rhogvz_mean_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
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
       !$omp parallel do default(none),private(g,k), &
       !$omp shared(l,gall,kall,ch,cmask,flx_h,rhog,EPS)
       do k = 1, kall
       do g = 1, gall
          ch(g,k,l,1) = flx_h(g,k,l,1) / rhog(g,k,l)
          ch(g,k,l,2) = flx_h(g,k,l,2) / rhog(g,k,l)
          ch(g,k,l,3) = flx_h(g,k,l,3) / rhog(g,k,l)
          ch(g,k,l,4) = flx_h(g,k,l,4) / rhog(g,k,l)
          ch(g,k,l,5) = flx_h(g,k,l,5) / rhog(g,k,l)
          ch(g,k,l,6) = flx_h(g,k,l,6) / rhog(g,k,l)

          ! c <= 0(incoming), cmask = 1
          cmask(g,k,l,1) = 0.5_RP - sign(0.5_RP,ch(g,k,l,1)-EPS)
          cmask(g,k,l,2) = 0.5_RP - sign(0.5_RP,ch(g,k,l,2)-EPS)
          cmask(g,k,l,3) = 0.5_RP - sign(0.5_RP,ch(g,k,l,3)-EPS)
          cmask(g,k,l,4) = 0.5_RP - sign(0.5_RP,ch(g,k,l,4)-EPS)
          cmask(g,k,l,5) = 0.5_RP - sign(0.5_RP,ch(g,k,l,5)-EPS)
          cmask(g,k,l,6) = 0.5_RP - sign(0.5_RP,ch(g,k,l,6)-EPS)
       enddo
       enddo
       !$omp end parallel do
    enddo

    if ( ADM_have_pl ) then
       g = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do v = ADM_gmin_pl, ADM_gmax_pl
          ch_pl(v,k,l) = flx_h_pl(v,k,l) / rhog_pl(g,k,l)

          cmask_pl(v,k,l) = 0.5_RP - sign(0.5_RP,ch_pl(v,k,l)-EPS)
       enddo
       enddo
       enddo
    endif

    do iq = 1, vmax

       !$omp parallel workshare
       q(:,:,:) = rhogq(:,:,:,iq) / rhog(:,:,:)
       !$omp end parallel workshare

       if ( ADM_have_pl ) then
          q_pl(:,:,:) = rhogq_pl(:,:,:,iq) / rhog_pl(:,:,:)
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

          !$omp parallel do default(none),private(g,k), &
          !$omp shared(l,iq,nstart,nend,kall,rhogq,flx_h,q_a)
          do k = 1, kall
          do g = nstart, nend
             rhogq(g,k,l,iq) = rhogq(g,k,l,iq) - ( flx_h(g,k,l,1) * q_a(g,k,l,1) &
                                                 + flx_h(g,k,l,2) * q_a(g,k,l,2) &
                                                 + flx_h(g,k,l,3) * q_a(g,k,l,3) &
                                                 + flx_h(g,k,l,4) * q_a(g,k,l,4) &
                                                 + flx_h(g,k,l,5) * q_a(g,k,l,5) &
                                                 + flx_h(g,k,l,6) * q_a(g,k,l,6) )
          enddo
          enddo
          !$omp end parallel do
       enddo

       if ( ADM_have_pl ) then
          g = ADM_gslf_pl

          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do v = ADM_gmin_pl, ADM_gmax_pl
             rhogq_pl(g,k,l,iq) = rhogq_pl(g,k,l,iq) - flx_h_pl(v,k,l) * q_a_pl(v,k,l)
          enddo
          enddo
          enddo
       endif

    enddo ! tracer q LOOP

    !--- update rhog
    do l = 1, ADM_lall
       nstart = suf(ADM_gmin,ADM_gmin)
       nend   = suf(ADM_gmax,ADM_gmax)

       !$omp parallel do default(none),private(g,k), &
       !$omp shared(l,iq,nstart,nend,kall,rhog,flx_h,frhog,dt)
       do k = 1, kall
       do g = nstart, nend
          rhog(g,k,l)= rhog(g,k,l) - ( flx_h(g,k,l,1) &
                                     + flx_h(g,k,l,2) &
                                     + flx_h(g,k,l,3) &
                                     + flx_h(g,k,l,4) &
                                     + flx_h(g,k,l,5) &
                                     + flx_h(g,k,l,6) ) + b2 * frhog(g,k,l) * dt
       enddo
       enddo
       !$omp end parallel do
    enddo

    if ( ADM_have_pl ) then
       g = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          do v = ADM_gmin_pl, ADM_gmax_pl
             rhog_pl(g,k,l)= rhog_pl(g,k,l) - flx_h_pl(v,k,l)
          enddo
          rhog_pl(g,k,l)= rhog_pl(g,k,l) + b2 * frhog_pl(g,k,l) * dt
       enddo
       enddo
    endif

    call PROF_rapend('____horizontal_adv',2)
    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 2nd
    !---------------------------------------------------------------------------
    call PROF_rapstart('____vertical_adv',2)

    do l = 1, ADM_lall
       !$omp parallel default(none),private(g,k), &
       !$omp shared(l,gall,kall,kmin,kmax,ck,d,frhog,rhog,flx_v,GRD_rdgz,dt)

       !$omp do
       do k = 1, kall
       do g = 1, gall
          d(g,k,l) = b3 * frhog(g,k,l) / rhog(g,k,l) * dt
       enddo
       enddo
       !$omp end do nowait

       !$omp do
       do k = kmin, kmax
       do g = 1, gall
          ck(g,k,l,1) = -flx_v(g,k  ,l) / rhog(g,k,l) * GRD_rdgz(k)
          ck(g,k,l,2) =  flx_v(g,k+1,l) / rhog(g,k,l) * GRD_rdgz(k)
       enddo
       enddo
       !$omp end do nowait

       !$omp do
       do g = 1, gall
          ck(g,kmin-1,l,1) = 0.0_RP
          ck(g,kmin-1,l,2) = 0.0_RP
          ck(g,kmax+1,l,1) = 0.0_RP
          ck(g,kmax+1,l,2) = 0.0_RP
       enddo
       !$omp end do

       !$omp end parallel
    enddo ! l LOOP

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          d_pl(:,:,l) = b3 * frhog_pl(:,:,l) / rhog_pl(:,:,l) * dt

          do k = ADM_kmin, ADM_kmax
             ck_pl(:,k,l,1) = -flx_v_pl(:,k  ,l) / rhog_pl(:,k,l) * GRD_rdgz(k)
             ck_pl(:,k,l,2) =  flx_v_pl(:,k+1,l) / rhog_pl(:,k,l) * GRD_rdgz(k)
          enddo
          ck_pl(:,ADM_kmin-1,l,1) = 0.0_RP
          ck_pl(:,ADM_kmin-1,l,2) = 0.0_RP
          ck_pl(:,ADM_kmax+1,l,1) = 0.0_RP
          ck_pl(:,ADM_kmax+1,l,2) = 0.0_RP
       enddo
    endif

    !--- vertical advection: 2nd-order centered difference
    do iq = 1, vmax

       do l = 1, ADM_lall
          !$omp parallel default(none),private(g,k), &
          !$omp shared(l,iq,gall,kall,kmin,kmax,q_h,q,rhogq,rhog,GRD_afac,GRD_bfac)

          !$omp do
          do k = 1, kall
          do g = 1, gall
             q(g,k,l) = rhogq(g,k,l,iq) / rhog(g,k,l)
          enddo
          enddo
          !$omp end do

          !$omp do
          do k = kmin, kmax+1
          do g = 1, gall
             q_h(g,k,l) = GRD_afac(k) * q(g,k,  l) &
                        + GRD_bfac(k) * q(g,k-1,l)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do g = 1, gall
             q_h(g,kmin-1,l) = 0.0_RP
          enddo
          !$omp end do

          !$omp end parallel
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             q_pl(:,:,l) = rhogq_pl(:,:,l,iq) / rhog_pl(:,:,l)

             do k = ADM_kmin, ADM_kmax+1
                q_h_pl(:,k,l) = GRD_afac(k) * q_pl(:,k,  l) &
                              + GRD_bfac(k) * q_pl(:,k-1,l)
             enddo
             q_h_pl(:,ADM_kmin-1,l) = 0.0_RP
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
          !$omp parallel default(none),private(g,k), &
          !$omp shared(l,iq,gall,kmin,kmax,rhogq,q_h,flx_v,GRD_rdgz)

          !$omp do
          do g = 1, gall
             q_h(g,kmin  ,l) = 0.0_RP
             q_h(g,kmax+1,l) = 0.0_RP
          enddo
          !$omp end do

          !$omp do
          do k = kmin, kmax
          do g = 1, gall
             rhogq(g,k,l,iq) = rhogq(g,k,l,iq) - ( flx_v(g,k+1,l) * q_h(g,k+1,l) &
                                                 - flx_v(g,k,  l) * q_h(g,k,  l) &
                                                 ) * GRD_rdgz(k)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do g = 1, gall
             rhogq(g,kmin-1,l,iq) = 0.0_RP
             rhogq(g,kmax+1,l,iq) = 0.0_RP
          enddo
          !$omp end do

          !$omp end parallel
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             q_h_pl(:,ADM_kmin  ,l) = 0.0_RP
             q_h_pl(:,ADM_kmax+1,l) = 0.0_RP

             do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                rhogq_pl(g,k,l,iq) = rhogq_pl(g,k,l,iq) - ( flx_v_pl(g,k+1,l)*q_h_pl(g,k+1,l) &
                                                          - flx_v_pl(g,k,  l)*q_h_pl(g,k,  l) &
                                                          ) * GRD_rdgz(k)
             enddo
             enddo

             do g = 1, ADM_gall_pl
                rhogq_pl(g,ADM_kmin-1,l,iq) = 0.0_RP
                rhogq_pl(g,ADM_kmax+1,l,iq) = 0.0_RP
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

    call PROF_rapend('____vertical_adv',2)

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
    use mod_const, only: &
       CONST_EPS
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
    use mod_grd, only: &
       GRD_xr,   &
       GRD_xr_pl
    use mod_gmtr, only: &
       GMTR_p,    &
       GMTR_p_pl, &
       GMTR_t,    &
       GMTR_t_pl, &
       GMTR_a,    &
       GMTR_a_pl
    implicit none

    real(RP), intent(out) :: flx_h    (ADM_gall   ,ADM_kall,ADM_lall   ,6)               ! horizontal mass flux
    real(RP), intent(out) :: flx_h_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(out) :: GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,AI:AJ,XDIR:ZDIR) ! mass centroid position
    real(RP), intent(out) :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)
    real(RP), intent(in)  :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )                 ! rho at cell center
    real(RP), intent(in)  :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: dt

    real(RP) :: rhot_TI  (ADM_gall   ) ! rho at cell vertex
    real(RP) :: rhot_TJ  (ADM_gall   ) ! rho at cell vertex
    real(RP) :: rhot_pl  (ADM_gall_pl)
    real(RP) :: rhovxt_TI(ADM_gall   )
    real(RP) :: rhovxt_TJ(ADM_gall   )
    real(RP) :: rhovxt_pl(ADM_gall_pl)
    real(RP) :: rhovyt_TI(ADM_gall   )
    real(RP) :: rhovyt_TJ(ADM_gall   )
    real(RP) :: rhovyt_pl(ADM_gall_pl)
    real(RP) :: rhovzt_TI(ADM_gall   )
    real(RP) :: rhovzt_TJ(ADM_gall   )
    real(RP) :: rhovzt_pl(ADM_gall_pl)

    real(RP) :: rhovxt2
    real(RP) :: rhovyt2
    real(RP) :: rhovzt2
    real(RP) :: flux
    real(RP) :: rrhoa2

    integer  :: gmin, gmax, kall, iall
    real(RP) :: EPS

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('____horizontal_adv_flux',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    EPS  = CONST_EPS

    do l = 1, ADM_lall
       !$omp parallel default(none), &
       !$omp private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,                                        &
       !$omp         rrhoa2,rhovxt2,rhovyt2,rhovzt2,flux),                                       &
       !$omp shared(l,ADM_have_sgp,gmin,gmax,kall,iall,rho,rhovx,rhovy,rhovz,flx_h,dt,           &
       !$omp        rhot_TI,rhovxt_TI,rhovyt_TI,rhovzt_TI,rhot_TJ,rhovxt_TJ,rhovyt_TJ,rhovzt_TJ, &
       !$omp        GRD_xc,GRD_xr,GMTR_p,GMTR_t,GMTR_a,EPS)
       do k = 1, kall

          ! (i,j),(i+1,j)
          !$omp do
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1

             rhot_TI  (ij) = rho  (ij  ,k,l) * GMTR_t(ij,K0,l,TI,W1) &
                           + rho  (ip1j,k,l) * GMTR_t(ij,K0,l,TI,W2)
             rhovxt_TI(ij) = rhovx(ij  ,k,l) * GMTR_t(ij,K0,l,TI,W1) &
                           + rhovx(ip1j,k,l) * GMTR_t(ij,K0,l,TI,W2)
             rhovyt_TI(ij) = rhovy(ij  ,k,l) * GMTR_t(ij,K0,l,TI,W1) &
                           + rhovy(ip1j,k,l) * GMTR_t(ij,K0,l,TI,W2)
             rhovzt_TI(ij) = rhovz(ij  ,k,l) * GMTR_t(ij,K0,l,TI,W1) &
                           + rhovz(ip1j,k,l) * GMTR_t(ij,K0,l,TI,W2)

             rhot_TJ  (ij) = rho  (ij  ,k,l) * GMTR_t(ij,K0,l,TJ,W1)
             rhovxt_TJ(ij) = rhovx(ij  ,k,l) * GMTR_t(ij,K0,l,TJ,W1)
             rhovyt_TJ(ij) = rhovy(ij  ,k,l) * GMTR_t(ij,K0,l,TJ,W1)
             rhovzt_TJ(ij) = rhovz(ij  ,k,l) * GMTR_t(ij,K0,l,TJ,W1)
          enddo
          enddo
          !$omp end do

          ! (i,j+1),(i+1,j+1)
          !$omp do
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ijp1   = ij + iall
             ip1jp1 = ij + iall + 1

             rhot_TI  (ij) = rhot_TI  (ij) + rho  (ip1jp1,k,l) * GMTR_t(ij,K0,l,TI,W3)
             rhovxt_TI(ij) = rhovxt_TI(ij) + rhovx(ip1jp1,k,l) * GMTR_t(ij,K0,l,TI,W3)
             rhovyt_TI(ij) = rhovyt_TI(ij) + rhovy(ip1jp1,k,l) * GMTR_t(ij,K0,l,TI,W3)
             rhovzt_TI(ij) = rhovzt_TI(ij) + rhovz(ip1jp1,k,l) * GMTR_t(ij,K0,l,TI,W3)

             rhot_TJ  (ij) = rhot_TJ  (ij) + rho  (ip1jp1,k,l) * GMTR_t(ij,K0,l,TJ,W2) &
                                           + rho  (ijp1  ,k,l) * GMTR_t(ij,K0,l,TJ,W3)
             rhovxt_TJ(ij) = rhovxt_TJ(ij) + rhovx(ip1jp1,k,l) * GMTR_t(ij,K0,l,TJ,W2) &
                                           + rhovx(ijp1  ,k,l) * GMTR_t(ij,K0,l,TJ,W3)
             rhovyt_TJ(ij) = rhovyt_TJ(ij) + rhovy(ip1jp1,k,l) * GMTR_t(ij,K0,l,TJ,W2) &
                                           + rhovy(ijp1  ,k,l) * GMTR_t(ij,K0,l,TJ,W3)
             rhovzt_TJ(ij) = rhovzt_TJ(ij) + rhovz(ip1jp1,k,l) * GMTR_t(ij,K0,l,TJ,W2) &
                                           + rhovz(ijp1  ,k,l) * GMTR_t(ij,K0,l,TJ,W3)
          enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then
             !$omp master
             j = gmin-1
             i = gmin-1

             ij   = (j-1)*iall + i
             ip1j = ij + 1

             rhot_TI  (ij) = rhot_TJ  (ip1j)
             rhovxt_TI(ij) = rhovxt_TJ(ip1j)
             rhovyt_TI(ij) = rhovyt_TJ(ip1j)
             rhovzt_TI(ij) = rhovzt_TJ(ip1j)
             !$omp end master
          endif

          !--- calculate flux and mass centroid position

          !$omp do
          do j = 1, iall
          do i = 1, iall
             if (      i < gmin .OR. i > gmax &
                  .OR. j < gmin .OR. j > gmax ) then
                ij = (j-1)*iall + i

                flx_h(ij,k,l,1) = 0.0_RP
                flx_h(ij,k,l,2) = 0.0_RP
                flx_h(ij,k,l,3) = 0.0_RP
                flx_h(ij,k,l,4) = 0.0_RP
                flx_h(ij,k,l,5) = 0.0_RP
                flx_h(ij,k,l,6) = 0.0_RP

                GRD_xc(ij,k,l,AI ,XDIR) = 0.0_RP
                GRD_xc(ij,k,l,AI ,YDIR) = 0.0_RP
                GRD_xc(ij,k,l,AI ,ZDIR) = 0.0_RP
                GRD_xc(ij,k,l,AIJ,XDIR) = 0.0_RP
                GRD_xc(ij,k,l,AIJ,YDIR) = 0.0_RP
                GRD_xc(ij,k,l,AIJ,ZDIR) = 0.0_RP
                GRD_xc(ij,k,l,AJ ,XDIR) = 0.0_RP
                GRD_xc(ij,k,l,AJ ,YDIR) = 0.0_RP
                GRD_xc(ij,k,l,AJ ,ZDIR) = 0.0_RP
             endif
          enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin  , gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1
             ijm1   = ij - iall

             rrhoa2  = 1.0_RP / max( rhot_TJ(ijm1) + rhot_TI(ij), EPS ) ! doubled
             rhovxt2 = rhovxt_TJ(ijm1) + rhovxt_TI(ij)
             rhovyt2 = rhovyt_TJ(ijm1) + rhovyt_TI(ij)
             rhovzt2 = rhovzt_TJ(ijm1) + rhovzt_TI(ij)

             flux = 0.5_RP * ( rhovxt2 * GMTR_a(ij,k0,l,AI ,HNX) &
                             + rhovyt2 * GMTR_a(ij,k0,l,AI ,HNY) &
                             + rhovzt2 * GMTR_a(ij,k0,l,AI ,HNZ) )

             flx_h(ij  ,k,l,1) =  flux * GMTR_p(ij  ,k0,l,P_RAREA) * dt
             flx_h(ip1j,k,l,4) = -flux * GMTR_p(ip1j,k0,l,P_RAREA) * dt

             GRD_xc(ij,k,l,AI,XDIR) = GRD_xr(ij,K0,l,AI,XDIR) - rhovxt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AI,YDIR) = GRD_xr(ij,K0,l,AI,YDIR) - rhovyt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AI,ZDIR) = GRD_xr(ij,K0,l,AI,ZDIR) - rhovzt2 * rrhoa2 * dt * 0.5_RP
          enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1jp1 = ij + iall + 1

             rrhoa2  = 1.0_RP / max( rhot_TI(ij) + rhot_TJ(ij), EPS ) ! doubled
             rhovxt2 = rhovxt_TI(ij) + rhovxt_TJ(ij)
             rhovyt2 = rhovyt_TI(ij) + rhovyt_TJ(ij)
             rhovzt2 = rhovzt_TI(ij) + rhovzt_TJ(ij)

             flux = 0.5_RP * ( rhovxt2 * GMTR_a(ij,k0,l,AIJ,HNX) &
                             + rhovyt2 * GMTR_a(ij,k0,l,AIJ,HNY) &
                             + rhovzt2 * GMTR_a(ij,k0,l,AIJ,HNZ) )

             flx_h(ij    ,k,l,2) =  flux * GMTR_p(ij    ,k0,l,P_RAREA) * dt
             flx_h(ip1jp1,k,l,5) = -flux * GMTR_p(ip1jp1,k0,l,P_RAREA) * dt

             GRD_xc(ij,k,l,AIJ,XDIR) = GRD_xr(ij,K0,l,AIJ,XDIR) - rhovxt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AIJ,YDIR) = GRD_xr(ij,K0,l,AIJ,YDIR) - rhovyt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AIJ,ZDIR) = GRD_xr(ij,K0,l,AIJ,ZDIR) - rhovzt2 * rrhoa2 * dt * 0.5_RP
          enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin-1, gmax
          do i = gmin  , gmax
             ij     = (j-1)*iall + i
             ijp1   = ij + iall
             im1j   = ij - 1

             rrhoa2  = 1.0_RP / max( rhot_TJ(ij) + rhot_TI(im1j), EPS ) ! doubled
             rhovxt2 = rhovxt_TJ(ij) + rhovxt_TI(im1j)
             rhovyt2 = rhovyt_TJ(ij) + rhovyt_TI(im1j)
             rhovzt2 = rhovzt_TJ(ij) + rhovzt_TI(im1j)

             flux = 0.5_RP * ( rhovxt2 * GMTR_a(ij,k0,l,AJ ,HNX) &
                             + rhovyt2 * GMTR_a(ij,k0,l,AJ ,HNY) &
                             + rhovzt2 * GMTR_a(ij,k0,l,AJ ,HNZ) )

             flx_h(ij  ,k,l,3) =  flux * GMTR_p(ij  ,k0,l,P_RAREA) * dt
             flx_h(ijp1,k,l,6) = -flux * GMTR_p(ijp1,k0,l,P_RAREA) * dt

             GRD_xc(ij,k,l,AJ,XDIR) = GRD_xr(ij,K0,l,AJ,XDIR) - rhovxt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AJ,YDIR) = GRD_xr(ij,K0,l,AJ,YDIR) - rhovyt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc(ij,k,l,AJ,ZDIR) = GRD_xr(ij,K0,l,AJ,ZDIR) - rhovzt2 * rrhoa2 * dt * 0.5_RP
          enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then
             !$omp master
             j = gmin
             i = gmin

             ij = (j-1)*iall + i

             flx_h(ij,k,l,6) = 0.0_RP
             !$omp end master
          endif

       enddo
       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl

             rhot_pl  (v) = rho_pl  (n   ,k,l) * GMTR_t_pl(ij,K0,l,W1) &
                          + rho_pl  (ij  ,k,l) * GMTR_t_pl(ij,K0,l,W2) &
                          + rho_pl  (ijp1,k,l) * GMTR_t_pl(ij,K0,l,W3)
             rhovxt_pl(v) = rhovx_pl(n   ,k,l) * GMTR_t_pl(ij,K0,l,W1) &
                          + rhovx_pl(ij  ,k,l) * GMTR_t_pl(ij,K0,l,W2) &
                          + rhovx_pl(ijp1,k,l) * GMTR_t_pl(ij,K0,l,W3)
             rhovyt_pl(v) = rhovy_pl(n   ,k,l) * GMTR_t_pl(ij,K0,l,W1) &
                          + rhovy_pl(ij  ,k,l) * GMTR_t_pl(ij,K0,l,W2) &
                          + rhovy_pl(ijp1,k,l) * GMTR_t_pl(ij,K0,l,W3)
             rhovzt_pl(v) = rhovz_pl(n   ,k,l) * GMTR_t_pl(ij,K0,l,W1) &
                          + rhovz_pl(ij  ,k,l) * GMTR_t_pl(ij,K0,l,W2) &
                          + rhovz_pl(ijp1,k,l) * GMTR_t_pl(ij,K0,l,W3)
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             rrhoa2  = 1.0_RP / max( rhot_pl(ijm1) + rhot_pl(ij), EPS ) ! doubled
             rhovxt2 = rhovxt_pl(ijm1) + rhovxt_pl(ij)
             rhovyt2 = rhovyt_pl(ijm1) + rhovyt_pl(ij)
             rhovzt2 = rhovzt_pl(ijm1) + rhovzt_pl(ij)

             flux = 0.5_RP * ( rhovxt2 * GMTR_a_pl(ij,K0,l,HNX) &
                             + rhovyt2 * GMTR_a_pl(ij,K0,l,HNY) &
                             + rhovzt2 * GMTR_a_pl(ij,K0,l,HNZ) )

             flx_h_pl(v,k,l) = flux * GMTR_p_pl(n,K0,l,P_RAREA) * dt

             GRD_xc_pl(v,k,l,XDIR) = GRD_xr_pl(v,K0,l,XDIR) - rhovxt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc_pl(v,k,l,YDIR) = GRD_xr_pl(v,K0,l,YDIR) - rhovyt2 * rrhoa2 * dt * 0.5_RP
             GRD_xc_pl(v,k,l,ZDIR) = GRD_xr_pl(v,K0,l,ZDIR) - rhovzt2 * rrhoa2 * dt * 0.5_RP
          enddo

       enddo
       enddo
    endif

    call PROF_rapend  ('____horizontal_adv_flux',2)

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
       ADM_nxyz,    &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl, &
       ADM_gmin_pl, &
       ADM_gmax_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_x,   &
       GRD_x_pl
    use mod_oprt, only: &
       OPRT_gradient,    &
       OPRT_coef_grad,   &
       OPRT_coef_grad_pl
    implicit none

    real(RP), intent(out) :: q_a      (ADM_gall   ,ADM_kall,ADM_lall   ,6)               ! q at cell face
    real(RP), intent(out) :: q_a_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)  :: q        (ADM_gall   ,ADM_kall,ADM_lall     )               ! q at cell center
    real(RP), intent(in)  :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)  :: cmask    (ADM_gall   ,ADM_kall,ADM_lall   ,6)               ! upwind direction mask
    real(RP), intent(in)  :: cmask_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)  :: GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,AI:AJ,XDIR:ZDIR) ! position of the mass centroid
    real(RP), intent(in)  :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)

    real(RP)  :: gradq   (ADM_gall   ,ADM_kall,ADM_lall   ,ADM_nxyz) ! grad(q)
    real(RP)  :: gradq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,ADM_nxyz)

    real(RP) :: q_ap1(ADM_gall), q_am1(ADM_gall)
    real(RP) :: q_ap2(ADM_gall), q_am2(ADM_gall)
    real(RP) :: q_ap3(ADM_gall), q_am3(ADM_gall)
    real(RP) :: q_ap4(ADM_gall), q_am4(ADM_gall)
    real(RP) :: q_ap5(ADM_gall), q_am5(ADM_gall)
    real(RP) :: q_ap6(ADM_gall), q_am6(ADM_gall)
    real(RP) :: q_ap, q_am

    integer :: kall, iall

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: nstart1, nstart2, nstart3, nstart4, nend
    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('____horizontal_adv_remap',2)

    call OPRT_gradient( gradq         (:,:,:,:), gradq_pl         (:,:,:,:), & ! [OUT]
                        q             (:,:,:),   q_pl             (:,:,:),   & ! [IN]
                        OPRT_coef_grad(:,:,:,:), OPRT_coef_grad_pl(:,:,:)    ) ! [IN]

    call COMM_data_transfer( gradq(:,:,:,:), gradq_pl(:,:,:,:) )

    kall = ADM_kall
    iall = ADM_gall_1d

    nstart1 = suf(ADM_gmin-1,ADM_gmin-1)
    nstart2 = suf(ADM_gmin  ,ADM_gmin-1)
    nstart3 = suf(ADM_gmin  ,ADM_gmin  )
    nstart4 = suf(ADM_gmin-1,ADM_gmin  )
    nend    = suf(ADM_gmax  ,ADM_gmax  )

    ! interpolated Q at cell arc
    do l = 1, ADM_lall

       !$omp parallel default(none),private(n,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(q,gradq,GRD_xc,GRD_x,l,kall,iall,nstart1,nstart2,nstart3,nstart4,nend, &
       !$omp        q_a,cmask,q_ap1,q_am1,q_ap2,q_am2,q_ap3,q_am3,q_ap4,q_am4,q_ap5,q_am5,q_ap6,q_am6)
       do k = 1, kall

          !$omp do
          do n = nstart1, nend
             ij = n

             q_ap1(n) = q(ij,k,l) + gradq(ij,k,l,XDIR) * ( GRD_xc(ij,k,l,AI,XDIR) - GRD_x(ij,K0,l,XDIR) ) &
                                  + gradq(ij,k,l,YDIR) * ( GRD_xc(ij,k,l,AI,YDIR) - GRD_x(ij,K0,l,YDIR) ) &
                                  + gradq(ij,k,l,ZDIR) * ( GRD_xc(ij,k,l,AI,ZDIR) - GRD_x(ij,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij   = n
             ip1j = n + 1

             q_am1(n) = q(ip1j,k,l) + gradq(ip1j,k,l,XDIR) * ( GRD_xc(ij,k,l,AI,XDIR) - GRD_x(ip1j,K0,l,XDIR) ) &
                                    + gradq(ip1j,k,l,YDIR) * ( GRD_xc(ij,k,l,AI,YDIR) - GRD_x(ip1j,K0,l,YDIR) ) &
                                    + gradq(ip1j,k,l,ZDIR) * ( GRD_xc(ij,k,l,AI,ZDIR) - GRD_x(ip1j,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij = n

             q_ap2(n) = q(ij,k,l) + gradq(ij,k,l,XDIR) * ( GRD_xc(ij,k,l,AIJ,XDIR) - GRD_x(ij,K0,l,XDIR) ) &
                                  + gradq(ij,k,l,YDIR) * ( GRD_xc(ij,k,l,AIJ,YDIR) - GRD_x(ij,K0,l,YDIR) ) &
                                  + gradq(ij,k,l,ZDIR) * ( GRD_xc(ij,k,l,AIJ,ZDIR) - GRD_x(ij,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij     = n
             ip1jp1 = n + 1 + iall

             q_am2(n) = q(ip1jp1,k,l) + gradq(ip1jp1,k,l,XDIR) * ( GRD_xc(ij,k,l,AIJ,XDIR) - GRD_x(ip1jp1,K0,l,XDIR) ) &
                                      + gradq(ip1jp1,k,l,YDIR) * ( GRD_xc(ij,k,l,AIJ,YDIR) - GRD_x(ip1jp1,K0,l,YDIR) ) &
                                      + gradq(ip1jp1,k,l,ZDIR) * ( GRD_xc(ij,k,l,AIJ,ZDIR) - GRD_x(ip1jp1,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij = n

             q_ap3(n) = q(ij,k,l) + gradq(ij,k,l,XDIR) * ( GRD_xc(ij,k,l,AJ,XDIR) - GRD_x(ij,K0,l,XDIR) ) &
                                  + gradq(ij,k,l,YDIR) * ( GRD_xc(ij,k,l,AJ,YDIR) - GRD_x(ij,K0,l,YDIR) ) &
                                  + gradq(ij,k,l,ZDIR) * ( GRD_xc(ij,k,l,AJ,ZDIR) - GRD_x(ij,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij   = n
             ijp1 = n + iall

             q_am3(n) = q(ijp1,k,l) + gradq(ijp1,k,l,XDIR) * ( GRD_xc(ij,k,l,AJ,XDIR) - GRD_x(ijp1,K0,l,XDIR) ) &
                                    + gradq(ijp1,k,l,YDIR) * ( GRD_xc(ij,k,l,AJ,YDIR) - GRD_x(ijp1,K0,l,YDIR) ) &
                                    + gradq(ijp1,k,l,ZDIR) * ( GRD_xc(ij,k,l,AJ,ZDIR) - GRD_x(ijp1,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart2, nend
             ij   = n
             im1j = n - 1

             q_ap4(n) = q(im1j,k,l) + gradq(im1j,k,l,XDIR) * ( GRD_xc(im1j,k,l,AI,XDIR) - GRD_x(im1j,K0,l,XDIR) ) &
                                    + gradq(im1j,k,l,YDIR) * ( GRD_xc(im1j,k,l,AI,YDIR) - GRD_x(im1j,K0,l,YDIR) ) &
                                    + gradq(im1j,k,l,ZDIR) * ( GRD_xc(im1j,k,l,AI,ZDIR) - GRD_x(im1j,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart2, nend
             ij = n

             q_am4(n) = q(ij,k,l) + gradq(ij,k,l,XDIR) * ( GRD_xc(im1j,k,l,AI,XDIR) - GRD_x(ij,K0,l,XDIR) ) &
                                  + gradq(ij,k,l,YDIR) * ( GRD_xc(im1j,k,l,AI,YDIR) - GRD_x(ij,K0,l,YDIR) ) &
                                  + gradq(ij,k,l,ZDIR) * ( GRD_xc(im1j,k,l,AI,ZDIR) - GRD_x(ij,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart3, nend
             ij     = n
             im1jm1 = n - 1 - iall

             q_ap5(n) = q(im1jm1,k,l) + gradq(im1jm1,k,l,XDIR) * ( GRD_xc(im1jm1,k,l,AIJ,XDIR) - GRD_x(im1jm1,K0,l,XDIR) ) &
                                      + gradq(im1jm1,k,l,YDIR) * ( GRD_xc(im1jm1,k,l,AIJ,YDIR) - GRD_x(im1jm1,K0,l,YDIR) ) &
                                      + gradq(im1jm1,k,l,ZDIR) * ( GRD_xc(im1jm1,k,l,AIJ,ZDIR) - GRD_x(im1jm1,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart3, nend
             ij = n

             q_am5(n) = q(ij,k,l) + gradq(ij,k,l,XDIR) * ( GRD_xc(im1jm1,k,l,AIJ,XDIR) - GRD_x(ij,K0,l,XDIR) ) &
                                  + gradq(ij,k,l,YDIR) * ( GRD_xc(im1jm1,k,l,AIJ,YDIR) - GRD_x(ij,K0,l,YDIR) ) &
                                  + gradq(ij,k,l,ZDIR) * ( GRD_xc(im1jm1,k,l,AIJ,ZDIR) - GRD_x(ij,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart4, nend
             ij   = n
             ijm1 = n - iall

             q_ap6(n) = q(ijm1,k,l) + gradq(ijm1,k,l,XDIR) * ( GRD_xc(ijm1,k,l,AJ,XDIR) - GRD_x(ijm1,K0,l,XDIR) ) &
                                    + gradq(ijm1,k,l,YDIR) * ( GRD_xc(ijm1,k,l,AJ,YDIR) - GRD_x(ijm1,K0,l,YDIR) ) &
                                    + gradq(ijm1,k,l,ZDIR) * ( GRD_xc(ijm1,k,l,AJ,ZDIR) - GRD_x(ijm1,K0,l,ZDIR) )
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart4, nend
             ij = n

             q_am6(n) = q(ij,k,l) + gradq(ij,k,l,XDIR) * ( GRD_xc(ijm1,k,l,AJ,XDIR) - GRD_x(ij,K0,l,XDIR) ) &
                                  + gradq(ij,k,l,YDIR) * ( GRD_xc(ijm1,k,l,AJ,YDIR) - GRD_x(ij,K0,l,YDIR) ) &
                                  + gradq(ij,k,l,ZDIR) * ( GRD_xc(ijm1,k,l,AJ,ZDIR) - GRD_x(ij,K0,l,ZDIR) )
          enddo
          !$omp end do

          !$omp do
          do n = nstart1, nstart2
             q_ap4(n) = 0.0_RP
             q_am4(n) = 0.0_RP
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nstart3
             q_ap5(n) = 0.0_RP
             q_am5(n) = 0.0_RP
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nstart4
             q_ap6(n) = 0.0_RP
             q_am6(n) = 0.0_RP
          enddo
          !$omp end do

          !$omp do
          do n = nstart1, nend
             q_a(n,k,l,1) = (        cmask(n,k,l,1) ) * q_am1(n) &
                          + ( 1.0_RP-cmask(n,k,l,1) ) * q_ap1(n)
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             q_a(n,k,l,2) = (        cmask(n,k,l,2) ) * q_am2(n) &
                          + ( 1.0_RP-cmask(n,k,l,2) ) * q_ap2(n)
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             q_a(n,k,l,3) = (        cmask(n,k,l,3) ) * q_am3(n) &
                          + ( 1.0_RP-cmask(n,k,l,3) ) * q_ap3(n)
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             q_a(n,k,l,4) = (        cmask(n,k,l,4) ) * q_am4(n) &
                          + ( 1.0_RP-cmask(n,k,l,4) ) * q_ap4(n)
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             q_a(n,k,l,5) = (        cmask(n,k,l,5) ) * q_am5(n) &
                          + ( 1.0_RP-cmask(n,k,l,5) ) * q_ap5(n)
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             q_a(n,k,l,6) = (        cmask(n,k,l,6) ) * q_am6(n) &
                          + ( 1.0_RP-cmask(n,k,l,6) ) * q_ap6(n)
          enddo
          !$omp end do

       enddo
       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do v = ADM_gmin_pl, ADM_gmax_pl
          q_ap = q_pl(n,k,l) + gradq_pl(n,k,l,XDIR) * ( GRD_xc_pl(v,k,l,XDIR) - GRD_x_pl(n,K0,l,XDIR) ) &
                             + gradq_pl(n,k,l,YDIR) * ( GRD_xc_pl(v,k,l,YDIR) - GRD_x_pl(n,K0,l,YDIR) ) &
                             + gradq_pl(n,k,l,ZDIR) * ( GRD_xc_pl(v,k,l,ZDIR) - GRD_x_pl(n,K0,l,ZDIR) )

          q_am = q_pl(v,k,l) + gradq_pl(v,k,l,XDIR) * ( GRD_xc_pl(v,k,l,XDIR) - GRD_x_pl(v,K0,l,XDIR) ) &
                             + gradq_pl(v,k,l,YDIR) * ( GRD_xc_pl(v,k,l,YDIR) - GRD_x_pl(v,K0,l,YDIR) ) &
                             + gradq_pl(v,k,l,ZDIR) * ( GRD_xc_pl(v,k,l,ZDIR) - GRD_x_pl(v,K0,l,ZDIR) )

          q_a_pl(v,k,l) = (        cmask_pl(v,k,l) ) * q_am &
                        + ( 1.0_RP-cmask_pl(v,k,l) ) * q_ap
       enddo
       enddo
       enddo
    endif

    call PROF_rapend  ('____horizontal_adv_remap',2)

    return
  end subroutine horizontal_remap

  !-----------------------------------------------------------------------------
  subroutine vertical_limiter_thuburn( &
       q_h, q_h_pl, &
       q,   q_pl,   &
       d,   d_pl,   &
       ck,  ck_pl   )
    use mod_const, only: &
       CONST_HUGE, &
       CONST_EPS
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    real(RP), intent(inout) :: q_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: q     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: d     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: ck    (ADM_gall   ,ADM_kall,ADM_lall   ,2)
    real(RP), intent(in)    :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(RP) :: Qout_min_k
    real(RP) :: Qout_max_k
    real(RP) :: Qout_min_km1(ADM_gall)
    real(RP) :: Qout_max_km1(ADM_gall)
    real(RP) :: Qout_min_pl(ADM_gall_pl,ADM_kall)
    real(RP) :: Qout_max_pl(ADM_gall_pl,ADM_kall)

    real(RP) :: Qin_minL, Qin_maxL
    real(RP) :: Qin_minU, Qin_maxU
    real(RP) :: qnext_min, qnext_max
    real(RP) :: Cin, Cout
    real(RP) :: CQin_min, CQin_max
    real(RP) :: inflagL, inflagU
    real(RP) :: zerosw

    integer  :: gall, kmin, kmax
    real(RP) :: EPS, BIG

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call PROF_rapstart('____vertical_adv_limiter',2)

    gall = ADM_gall
    kmin = ADM_kmin
    kmax = ADM_kmax

    EPS  = CONST_EPS
    BIG  = CONST_HUGE

    do l = 1, ADM_lall
       !$omp parallel default(none), &
       !$omp private(g,k,zerosw,inflagL,inflagU,Qin_minL,Qin_minU,Qin_maxL,Qin_maxU, &
       !$omp         qnext_min,qnext_max,Cin,Cout,CQin_min,CQin_max,Qout_min_k,Qout_max_k),  &
       !$omp shared(l,gall,kmin,kmax,q_h,ck,q,d,Qout_min_km1,Qout_max_km1,EPS,BIG)

       !$omp do
       do g = 1, gall
          k = kmin ! peeling

          inflagL = 0.5_RP - sign(0.5_RP, ck(g,k  ,l,1)) ! incoming flux: flag=1
          inflagU = 0.5_RP - sign(0.5_RP,-ck(g,k+1,l,1)) ! incoming flux: flag=1

          Qin_minL = min( q(g,k,l), q(g,k-1,l) ) + ( 1.0_RP-inflagL ) * BIG
          Qin_minU = min( q(g,k,l), q(g,k+1,l) ) + ( 1.0_RP-inflagU ) * BIG
          Qin_maxL = max( q(g,k,l), q(g,k-1,l) ) - ( 1.0_RP-inflagL ) * BIG
          Qin_maxU = max( q(g,k,l), q(g,k+1,l) ) - ( 1.0_RP-inflagU ) * BIG

          qnext_min = min( Qin_minL, Qin_minU, q(g,k,l) )
          qnext_max = max( Qin_maxL, Qin_maxU, q(g,k,l) )

          Cin      = (        inflagL ) * ck(g,k,l,1) &
                   + (        inflagU ) * ck(g,k,l,2)
          Cout     = ( 1.0_RP-inflagL ) * ck(g,k,l,1) &
                   + ( 1.0_RP-inflagU ) * ck(g,k,l,2)

          CQin_min = (        inflagL ) * ck(g,k,l,1) * Qin_minL &
                   + (        inflagU ) * ck(g,k,l,2) * Qin_minU
          CQin_max = (        inflagL ) * ck(g,k,l,1) * Qin_maxL &
                   + (        inflagU ) * ck(g,k,l,2) * Qin_maxU

          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-EPS) ! if Cout = 0, sw = 1

          Qout_min_k = ( q(g,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(g,k,l)) ) &
                     / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                     + q(g,k,l) * zerosw
          Qout_max_k = ( q(g,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(g,k,l)) ) &
                     / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                     + q(g,k,l) * zerosw

          Qout_min_km1(g) = Qout_min_k
          Qout_max_km1(g) = Qout_max_k
       enddo
       !$omp end do

       do k = kmin+1, kmax
          !$omp do
          do g = 1, gall
             inflagL = 0.5_RP - sign(0.5_RP, ck(g,k  ,l,1)) ! incoming flux: flag=1
             inflagU = 0.5_RP - sign(0.5_RP,-ck(g,k+1,l,1)) ! incoming flux: flag=1

             Qin_minL = min( q(g,k,l), q(g,k-1,l) ) + ( 1.0_RP-inflagL ) * BIG
             Qin_minU = min( q(g,k,l), q(g,k+1,l) ) + ( 1.0_RP-inflagU ) * BIG
             Qin_maxL = max( q(g,k,l), q(g,k-1,l) ) - ( 1.0_RP-inflagL ) * BIG
             Qin_maxU = max( q(g,k,l), q(g,k+1,l) ) - ( 1.0_RP-inflagU ) * BIG

             qnext_min = min( Qin_minL, Qin_minU, q(g,k,l) )
             qnext_max = max( Qin_maxL, Qin_maxU, q(g,k,l) )

             Cin      = (        inflagL ) * ck(g,k,l,1) &
                      + (        inflagU ) * ck(g,k,l,2)
             Cout     = ( 1.0_RP-inflagL ) * ck(g,k,l,1) &
                      + ( 1.0_RP-inflagU ) * ck(g,k,l,2)

             CQin_min = (        inflagL ) * ck(g,k,l,1) * Qin_minL &
                      + (        inflagU ) * ck(g,k,l,2) * Qin_minU
             CQin_max = (        inflagL ) * ck(g,k,l,1) * Qin_maxL &
                      + (        inflagU ) * ck(g,k,l,2) * Qin_maxU

             zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-EPS) ! if Cout = 0, sw = 1

             Qout_min_k = ( q(g,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(g,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q(g,k,l) * zerosw
             Qout_max_k = ( q(g,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(g,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q(g,k,l) * zerosw

             q_h(g,k,l) = (        inflagL ) * max( min( q_h(g,k,l), Qout_max_km1(g) ), Qout_min_km1(g) ) &
                        + ( 1.0_RP-inflagL ) * max( min( q_h(g,k,l), Qout_max_k      ), Qout_min_k      )

             Qout_min_km1(g) = Qout_min_k
             Qout_max_km1(g) = Qout_max_k
          enddo
          !$omp end do
       enddo

       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl

          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             inflagL = 0.5_RP - sign(0.5_RP, ck_pl(g,k  ,l,1)) ! incoming flux: flag=1
             inflagU = 0.5_RP - sign(0.5_RP,-ck_pl(g,k+1,l,1)) ! incoming flux: flag=1

             Qin_minL = min( q_pl(g,k,l), q_pl(g,k-1,l) ) + ( 1.0_RP-inflagL ) * BIG
             Qin_minU = min( q_pl(g,k,l), q_pl(g,k+1,l) ) + ( 1.0_RP-inflagU ) * BIG
             Qin_maxL = max( q_pl(g,k,l), q_pl(g,k-1,l) ) - ( 1.0_RP-inflagL ) * BIG
             Qin_maxU = max( q_pl(g,k,l), q_pl(g,k+1,l) ) - ( 1.0_RP-inflagU ) * BIG

             qnext_min = min( Qin_minL, Qin_minU, q_pl(g,k,l) )
             qnext_max = max( Qin_maxL, Qin_maxU, q_pl(g,k,l) )

             Cin      = (        inflagL ) * ( ck_pl(g,k,l,1) ) &
                      + (        inflagU ) * ( ck_pl(g,k,l,2) )
             Cout     = ( 1.0_RP-inflagL ) * ( ck_pl(g,k,l,1) ) &
                      + ( 1.0_RP-inflagU ) * ( ck_pl(g,k,l,2) )

             CQin_max = (        inflagL ) * ( ck_pl(g,k,l,1) * Qin_maxL ) &
                      + (        inflagU ) * ( ck_pl(g,k,l,2) * Qin_maxU )
             CQin_min = (        inflagL ) * ( ck_pl(g,k,l,1) * Qin_minL ) &
                      + (        inflagU ) * ( ck_pl(g,k,l,2) * Qin_minU )

             zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-EPS) ! if Cout = 0, sw = 1

             Qout_min_pl(g,k) = ( q_pl(g,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d_pl(g,k,l)) ) &
                              / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                              &
                              + q_pl(g,k,l) * zerosw
             Qout_max_pl(g,k) = ( q_pl(g,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d_pl(g,k,l)) ) &
                              / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                              &
                              + q_pl(g,k,l) * zerosw
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             inflagL = 0.5_RP - sign(0.5_RP,ck_pl(g,k,l,1)) ! incoming flux: flag=1

             q_h_pl(g,k,l) = (        inflagL ) * max( min( q_h_pl(g,k,l), Qout_max_pl(g,k-1) ), Qout_min_pl(g,k-1) ) &
                           + ( 1.0_RP-inflagL ) * max( min( q_h_pl(g,k,l), Qout_max_pl(g,k  ) ), Qout_min_pl(g,k  ) )
          enddo
          enddo

       enddo
    endif

    call PROF_rapend  ('____vertical_adv_limiter',2)

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
    use mod_const, only: &
       CONST_HUGE, &
       CONST_EPS
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
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(RP), intent(inout) :: q_a     (ADM_gall   ,ADM_kall,ADM_lall   ,6)
    real(RP), intent(inout) :: q_a_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)    :: q       (ADM_gall   ,ADM_kall,ADM_lall     )
    real(RP), intent(in)    :: q_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)    :: d       (ADM_gall   ,ADM_kall,ADM_lall     )
    real(RP), intent(in)    :: d_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)    :: ch      (ADM_gall   ,ADM_kall,ADM_lall   ,6)
    real(RP), intent(in)    :: ch_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl  )
    real(RP), intent(in)    :: cmask   (ADM_gall   ,ADM_kall,ADM_lall   ,6)
    real(RP), intent(in)    :: cmask_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl  )

    real(RP) :: q_min_AI, q_min_AIJ, q_min_AJ, q_min_pl
    real(RP) :: q_max_AI, q_max_AIJ, q_max_AJ, q_max_pl

    real(RP) :: qnext_min   , qnext_min_pl
    real(RP) :: qnext_max   , qnext_max_pl
    real(RP) :: Cin_sum     , Cin_sum_pl
    real(RP) :: Cout_sum    , Cout_sum_pl
    real(RP) :: CQin_max_sum, CQin_max_sum_pl
    real(RP) :: CQin_min_sum, CQin_min_sum_pl

    integer, parameter :: I_min = 1
    integer, parameter :: I_max = 2
    real(RP) :: Qin    (ADM_gall   ,ADM_kall,ADM_lall   ,2,6)
    real(RP) :: Qin_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2,2)
    real(RP) :: Qout   (ADM_gall   ,ADM_kall,ADM_lall   ,2  )
    real(RP) :: Qout_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2  )

    real(RP) :: ch_masked1
    real(RP) :: ch_masked2
    real(RP) :: ch_masked3
    real(RP) :: ch_masked4
    real(RP) :: ch_masked5
    real(RP) :: ch_masked6
    real(RP) :: ch_masked
    real(RP) :: zerosw

    integer  :: gmin, gmax, kall, iall
    real(RP) :: EPS, BIG

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1, ip2jp1
    integer :: im1j, ijm1

    integer :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('____horizontal_adv_limiter',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    EPS  = CONST_EPS
    BIG  = CONST_HUGE

    do l = 1, ADM_lall
       !$omp parallel default(none), &
       !$omp private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,ip2jp1,                        &
       !$omp         q_min_AI,q_min_AIJ,q_min_AJ,q_max_AI,q_max_AIJ,q_max_AJ,zerosw,    &
       !$omp         ch_masked1,ch_masked2,ch_masked3,ch_masked4,ch_masked5,ch_masked6, &
       !$omp         qnext_min,qnext_max,Cin_sum,Cout_sum,CQin_min_sum,CQin_max_sum),   &
       !$omp shared(l,ADM_have_sgp,gmin,gmax,kall,iall,q,cmask,d,ch,Qin,Qout,EPS,BIG)
       do k = 1, kall
          !---< (i) define inflow bounds, eq.(32)&(33) >---
          !$omp do
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1
             ip1jp1 = ij + iall + 1
             ijp1   = ij + iall
             im1j   = ij - 1
             ijm1   = ij - iall

             im1j   = max( im1j  , 1 )
             ijm1   = max( ijm1  , 1 )

             q_min_AI  = min( q(ij,k,l), q(ijm1,k,l), q(ip1j,k,l), q(ip1jp1,k,l) )
             q_max_AI  = max( q(ij,k,l), q(ijm1,k,l), q(ip1j,k,l), q(ip1jp1,k,l) )
             q_min_AIJ = min( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijp1,k,l) )
             q_max_AIJ = max( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijp1,k,l) )
             q_min_AJ  = min( q(ij,k,l), q(ip1jp1,k,l), q(ijp1,k,l), q(im1j,k,l) )
             q_max_AJ  = max( q(ij,k,l), q(ip1jp1,k,l), q(ijp1,k,l), q(im1j,k,l) )

             Qin(ij,    k,l,I_min,1) = (        cmask(ij,k,l,1) ) * q_min_AI &
                                     + ( 1.0_RP-cmask(ij,k,l,1) ) * BIG
             Qin(ip1j,  k,l,I_min,4) = (        cmask(ij,k,l,1) ) * BIG      &
                                     + ( 1.0_RP-cmask(ij,k,l,1) ) * q_min_AI
             Qin(ij,    k,l,I_max,1) = (        cmask(ij,k,l,1) ) * q_max_AI &
                                     + ( 1.0_RP-cmask(ij,k,l,1) ) * (-BIG)
             Qin(ip1j,  k,l,I_max,4) = (        cmask(ij,k,l,1) ) * (-BIG)   &
                                     + ( 1.0_RP-cmask(ij,k,l,1) ) * q_max_AI

             Qin(ij,    k,l,I_min,2) = (        cmask(ij,k,l,2) ) * q_min_AIJ &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * BIG
             Qin(ip1jp1,k,l,I_min,5) = (        cmask(ij,k,l,2) ) * BIG       &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * q_min_AIJ
             Qin(ij,    k,l,I_max,2) = (        cmask(ij,k,l,2) ) * q_max_AIJ &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * (-BIG)
             Qin(ip1jp1,k,l,I_max,5) = (        cmask(ij,k,l,2) ) * (-BIG)    &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * q_max_AIJ

             Qin(ij,    k,l,I_min,3) = (        cmask(ij,k,l,3) ) * q_min_AJ &
                                     + ( 1.0_RP-cmask(ij,k,l,3) ) * BIG
             Qin(ijp1,  k,l,I_min,6) = (        cmask(ij,k,l,3) ) * BIG      &
                                     + ( 1.0_RP-cmask(ij,k,l,3) ) * q_min_AJ
             Qin(ij,    k,l,I_max,3) = (        cmask(ij,k,l,3) ) * q_max_AJ &
                                     + ( 1.0_RP-cmask(ij,k,l,3) ) * (-BIG)
             Qin(ijp1,  k,l,I_max,6) = (        cmask(ij,k,l,3) ) * (-BIG)   &
                                     + ( 1.0_RP-cmask(ij,k,l,3) ) * q_max_AJ
          enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then
             !$omp master
             j = gmin-1
             i = gmin-1

             ij     = (j-1)*iall + i
             ijp1   = ij + iall
             ip1jp1 = ij + iall + 1
             ip2jp1 = ij + iall + 2

             q_min_AIJ = min( q(ij,k,l), q(ip1jp1,k,l), q(ip2jp1,k,l), q(ijp1,k,l) )
             q_max_AIJ = max( q(ij,k,l), q(ip1jp1,k,l), q(ip2jp1,k,l), q(ijp1,k,l) )

             Qin(ij,    k,l,I_min,2) = (        cmask(ij,k,l,2) ) * q_min_AIJ &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * BIG
             Qin(ip1jp1,k,l,I_min,5) = (        cmask(ij,k,l,2) ) * BIG       &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * q_min_AIJ
             Qin(ij,    k,l,I_max,2) = (        cmask(ij,k,l,2) ) * q_max_AIJ &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * (-BIG)
             Qin(ip1jp1,k,l,I_max,5) = (        cmask(ij,k,l,2) ) * (-BIG)    &
                                     + ( 1.0_RP-cmask(ij,k,l,2) ) * q_max_AIJ
             !$omp end master
          endif

          !---< (iii) define allowable range of q at next step, eq.(42)&(43) >---
          !$omp do
          do j = gmin, gmax
          do i = gmin, gmax
             ij = (j-1)*iall + i

             qnext_min = min( q(ij,k,l),           &
                              Qin(ij,k,l,I_min,1), &
                              Qin(ij,k,l,I_min,2), &
                              Qin(ij,k,l,I_min,3), &
                              Qin(ij,k,l,I_min,4), &
                              Qin(ij,k,l,I_min,5), &
                              Qin(ij,k,l,I_min,6)  )

             qnext_max = max( q(ij,k,l),           &
                              Qin(ij,k,l,I_max,1), &
                              Qin(ij,k,l,I_max,2), &
                              Qin(ij,k,l,I_max,3), &
                              Qin(ij,k,l,I_max,4), &
                              Qin(ij,k,l,I_max,5), &
                              Qin(ij,k,l,I_max,6)  )

             ch_masked1 = min( ch(ij,k,l,1), 0.0_RP )
             ch_masked2 = min( ch(ij,k,l,2), 0.0_RP )
             ch_masked3 = min( ch(ij,k,l,3), 0.0_RP )
             ch_masked4 = min( ch(ij,k,l,4), 0.0_RP )
             ch_masked5 = min( ch(ij,k,l,5), 0.0_RP )
             ch_masked6 = min( ch(ij,k,l,6), 0.0_RP )

             Cin_sum      = ch_masked1 &
                          + ch_masked2 &
                          + ch_masked3 &
                          + ch_masked4 &
                          + ch_masked5 &
                          + ch_masked6

             Cout_sum     = ch(ij,k,l,1) - ch_masked1 &
                          + ch(ij,k,l,2) - ch_masked2 &
                          + ch(ij,k,l,3) - ch_masked3 &
                          + ch(ij,k,l,4) - ch_masked4 &
                          + ch(ij,k,l,5) - ch_masked5 &
                          + ch(ij,k,l,6) - ch_masked6

             CQin_min_sum = ch_masked1 * Qin(ij,k,l,I_min,1) &
                          + ch_masked2 * Qin(ij,k,l,I_min,2) &
                          + ch_masked3 * Qin(ij,k,l,I_min,3) &
                          + ch_masked4 * Qin(ij,k,l,I_min,4) &
                          + ch_masked5 * Qin(ij,k,l,I_min,5) &
                          + ch_masked6 * Qin(ij,k,l,I_min,6)

             CQin_max_sum = ch_masked1 * Qin(ij,k,l,I_max,1) &
                          + ch_masked2 * Qin(ij,k,l,I_max,2) &
                          + ch_masked3 * Qin(ij,k,l,I_max,3) &
                          + ch_masked4 * Qin(ij,k,l,I_max,4) &
                          + ch_masked5 * Qin(ij,k,l,I_max,5) &
                          + ch_masked6 * Qin(ij,k,l,I_max,6)

             zerosw = 0.5_RP - sign(0.5_RP,abs(Cout_sum)-EPS) ! if Cout_sum = 0, sw = 1

             Qout(ij,k,l,I_min) = ( q(ij,k,l) - CQin_max_sum - qnext_max*(1.0_RP-Cin_sum-Cout_sum+d(ij,k,l)) ) &
                                / ( Cout_sum + zerosw ) * ( 1.0_RP - zerosw )                                         &
                                + q(ij,k,l) * zerosw
             Qout(ij,k,l,I_max) = ( q(ij,k,l) - CQin_min_sum - qnext_min*(1.0_RP-Cin_sum-Cout_sum+d(ij,k,l)) ) &
                                / ( Cout_sum + zerosw ) * ( 1.0_RP - zerosw )                                         &
                                + q(ij,k,l) * zerosw
          enddo
          enddo
          !$omp end do

          !$omp do
          do j = 1, iall
          do i = 1, iall
             if (      i < gmin .OR. i > gmax &
                  .OR. j < gmin .OR. j > gmax ) then
                ij = (j-1)*iall + i

                Qout(ij,k,l,I_min) = q(ij,k,l)
                Qout(ij,k,l,I_min) = q(ij,k,l)
                Qout(ij,k,l,I_max) = q(ij,k,l)
                Qout(ij,k,l,I_max) = q(ij,k,l)
             endif
          enddo
          enddo
          !$omp end do

       enddo ! k loop
       !$omp end parallel
    enddo ! l loop

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

                Qin_pl(ij,k,l,I_min,1) = (        cmask_pl(ij,k,l) ) * q_min_pl &
                                       + ( 1.0_RP-cmask_pl(ij,k,l) ) * BIG
                Qin_pl(ij,k,l,I_min,2) = (        cmask_pl(ij,k,l) ) * BIG      &
                                       + ( 1.0_RP-cmask_pl(ij,k,l) ) * q_min_pl
                Qin_pl(ij,k,l,I_max,1) = (        cmask_pl(ij,k,l) ) * q_max_pl &
                                       + ( 1.0_RP-cmask_pl(ij,k,l) ) * (-BIG)
                Qin_pl(ij,k,l,I_max,2) = (        cmask_pl(ij,k,l) ) * (-BIG)   &
                                       + ( 1.0_RP-cmask_pl(ij,k,l) ) * q_max_pl
             enddo

             qnext_min_pl = q_pl(n,k,l)
             qnext_max_pl = q_pl(n,k,l)
             do v = ADM_gmin_pl, ADM_gmax_pl
                qnext_min_pl = min( qnext_min_pl, Qin_pl(v,k,l,I_min,1) )
                qnext_max_pl = max( qnext_max_pl, Qin_pl(v,k,l,I_max,1) )
             enddo

             Cin_sum_pl      = 0.0_RP
             Cout_sum_pl     = 0.0_RP
             CQin_max_sum_pl = 0.0_RP
             CQin_min_sum_pl = 0.0_RP
             do v = ADM_gmin_pl, ADM_gmax_pl
                ch_masked = cmask_pl(v,k,l) * ch_pl(v,k,l)

                Cin_sum_pl      = Cin_sum_pl      + ch_masked
                Cout_sum_pl     = Cout_sum_pl     - ch_masked + ch_pl (v,k,l)
                CQin_min_sum_pl = CQin_min_sum_pl + ch_masked * Qin_pl(v,k,l,I_min,1)
                CQin_max_sum_pl = CQin_max_sum_pl + ch_masked * Qin_pl(v,k,l,I_max,1)
             enddo

             zerosw = 0.5_RP - sign(0.5_RP,abs(Cout_sum_pl)-EPS) ! if Cout_sum_pl = 0, sw = 1

             Qout_pl(n,k,l,I_min) = ( q_pl(n,k,l) - CQin_max_sum_pl - qnext_max_pl*(1.0_RP-Cin_sum_pl-Cout_sum_pl+d_pl(n,k,l)) ) &
                                  / ( Cout_sum_pl + zerosw ) * ( 1.0_RP - zerosw )                                               &
                                  + q_pl(n,k,l) * zerosw
             Qout_pl(n,k,l,I_max) = ( q_pl(n,k,l) - CQin_min_sum_pl - qnext_min_pl*(1.0_RP-Cin_sum_pl-Cout_sum_pl+d_pl(n,k,l)) ) &
                                  / ( Cout_sum_pl + zerosw ) * ( 1.0_RP - zerosw )                                               &
                                  + q_pl(n,k,l) * zerosw
          enddo
       enddo
    endif

    call COMM_data_transfer( Qout(:,:,:,:), Qout_pl(:,:,:,:) )

    !---- apply inflow/outflow limiter
    do l = 1, ADM_lall
       !$omp parallel do default(none), private(i,j,k,ij,ip1j,ip1jp1,ijp1), &
       !$omp shared(l,gmin,gmax,kall,iall,q_a,cmask,Qin,Qout)
       do k = 1, kall
          do j = gmin-1, gmax
          do i = gmin-1, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1
             ip1jp1 = ij + iall + 1
             ijp1   = ij + iall

             q_a(ij,k,l,1) = (        cmask(ij,k,l,1) ) * min( max( q_a(ij,k,l,1), Qin (ij    ,k,l,I_min,1) ), Qin (ij    ,k,l,I_max,1) ) &
                           + ( 1.0_RP-cmask(ij,k,l,1) ) * min( max( q_a(ij,k,l,1), Qin (ip1j  ,k,l,I_min,4) ), Qin (ip1j  ,k,l,I_max,4) )
             q_a(ij,k,l,1) = (        cmask(ij,k,l,1) ) * max( min( q_a(ij,k,l,1), Qout(ip1j  ,k,l,I_max  ) ), Qout(ip1j  ,k,l,I_min  ) ) &
                           + ( 1.0_RP-cmask(ij,k,l,1) ) * max( min( q_a(ij,k,l,1), Qout(ij    ,k,l,I_max  ) ), Qout(ij    ,k,l,I_min  ) )
             q_a(ip1j,k,l,4) = q_a(ij,k,l,1)

             q_a(ij,k,l,2) = (        cmask(ij,k,l,2) ) * min( max( q_a(ij,k,l,2), Qin (ij    ,k,l,I_min,2) ), Qin (ij    ,k,l,I_max,2) ) &
                           + ( 1.0_RP-cmask(ij,k,l,2) ) * min( max( q_a(ij,k,l,2), Qin (ip1jp1,k,l,I_min,5) ), Qin (ip1jp1,k,l,I_max,5) )
             q_a(ij,k,l,2) = (        cmask(ij,k,l,2) ) * max( min( q_a(ij,k,l,2), Qout(ip1jp1,k,l,I_max  ) ), Qout(ip1jp1,k,l,I_min  ) ) &
                           + ( 1.0_RP-cmask(ij,k,l,2) ) * max( min( q_a(ij,k,l,2), Qout(ij    ,k,l,I_max  ) ), Qout(ij    ,k,l,I_min  ) )
             q_a(ip1jp1,k,l,5) = q_a(ij,k,l,2)

             q_a(ij,k,l,3) = (        cmask(ij,k,l,3) ) * min( max( q_a(ij,k,l,3), Qin (ij    ,k,l,I_min,3) ), Qin (ij    ,k,l,I_max,3) ) &
                           + ( 1.0_RP-cmask(ij,k,l,3) ) * min( max( q_a(ij,k,l,3), Qin (ijp1  ,k,l,I_min,6) ), Qin (ijp1  ,k,l,I_max,6) )
             q_a(ij,k,l,3) = (        cmask(ij,k,l,3) ) * max( min( q_a(ij,k,l,3), Qout(ijp1  ,k,l,I_max  ) ), Qout(ijp1  ,k,l,I_min  ) ) &
                           + ( 1.0_RP-cmask(ij,k,l,3) ) * max( min( q_a(ij,k,l,3), Qout(ij    ,k,l,I_max  ) ), Qout(ij    ,k,l,I_min  ) )
             q_a(ijp1,k,l,6) = q_a(ij,k,l,3)
          enddo
          enddo
       enddo
       !$omp end parallel do
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do v = ADM_gmin_pl, ADM_gmax_pl
          q_a_pl(v,k,l) = (        cmask_pl(v,k,l) ) * min(max(q_a_pl(v,k,l), Qin_pl (v,k,l,I_min,1)), Qin_pl (v,k,l,I_max,1)) &
                        + ( 1.0_RP-cmask_pl(v,k,l) ) * min(max(q_a_pl(v,k,l), Qin_pl (v,k,l,I_min,2)), Qin_pl (v,k,l,I_max,2))
          q_a_pl(v,k,l) = (        cmask_pl(v,k,l) ) * max(min(q_a_pl(v,k,l), Qout_pl(v,k,l,I_max  )), Qout_pl(v,k,l,I_min  )) &
                        + ( 1.0_RP-cmask_pl(v,k,l) ) * max(min(q_a_pl(v,k,l), Qout_pl(n,k,l,I_max  )), Qout_pl(n,k,l,I_min  ))
       enddo
       enddo
       enddo
    endif

    call PROF_rapend  ('____horizontal_adv_limiter',2)

    return
  end subroutine horizontal_limiter_thuburn

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_src_tracer
