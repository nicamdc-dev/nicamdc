  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module is for limitting values of q at cell walls in advection
  !       calculation.
  !       (Reference: Thuburn, 1996) 
  !        These are imported from sub[src_update_tracer/mod_src] and 
  !        sub[OPRT_divergence2/mod_oprt].
  !
  !++ Current COrresponding Author: Y.Niwa
  !
  !++ History:
  !      Version   Date      Comment
  !      -----------------------------------------------------------------------
  !      0.00      08-01-24  Imported from mod_src and mod_oprt
  !                11-09-27  T.Seiki: merge optimization for K by M.Terai and RIST
  !                11-11-28  Y.Yamada: Merge Terai-san timer code
  !                                           into the original code.
  !      -----------------------------------------------------------------------
  !
module mod_trcadv_thuburn
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: src_update_tracer
  public :: advlim_thuburn_v
  public :: OPRT_divergence2_prep
  public :: OPRT_divergence2

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
  ! < for OPRT_divergence2_prep >  ! Y.Niwa add 080130
  real(8), private, allocatable, save :: local_t_var(:,:,:,:,:)
  real(8), private, allocatable, save :: local_t_var_pl(:,:,:,:)

  !-----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------------
  ! Y.Niwa add 080124
  ! This routine is revised version of src_update_tracer
  subroutine src_update_tracer( &
       nqmax,                       & !--- IN    : number of tracers
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
       ADM_gall,       &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_kall,       &
       ADM_lall,       &
       ADM_gall_pl,    &
       ADM_lall_pl,    &
       ADM_kmin,       &
       ADM_kmax,       &
       ADM_prc_me,     &
       ADM_prc_pl,     &
       ADM_AI,ADM_AJ,  &
       ADM_GSLF_PL,    &
       ADM_GMIN_PL,    &
       ADM_gmax_pl,    &
       ADM_gall_1d
    use mod_vmtr, only : &
       VMTR_RGSGAM2,   &
       VMTR_RGSGAM2_pl,&
       VMTR_GZXH,      &
       VMTR_GZXH_pl,   &
       VMTR_GZYH,      &
       VMTR_GZYH_pl,   &
       VMTR_GZZH,      &
       VMTR_GZZH_pl,   &
       VMTR_GSGAMH,    &
       VMTR_GSGAMH_pl, &
       VMTR_RGSH,      &
       VMTR_RGSH_pl,   &
       VMTR_RGAM,      &
       VMTR_RGAM_pl
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac, &
       GRD_XDIR, &
       GRD_ZDIR, &
       GRD_rdgz
    implicit none

    integer, intent(in)    :: nqmax
    real(8), intent(inout) :: rhogq         (ADM_gall,   ADM_kall,ADM_lall,   nqmax)
    real(8), intent(inout) :: rhogq_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,nqmax)
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

    real(8) :: rhog    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rrhog   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rrhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: q       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: q_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: q_h     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: q_h_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: d       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: d_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: flx_v   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_v_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ck      (ADM_gall,   ADM_kall,ADM_lall,   2)
    real(8) :: ck_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(8) :: flx_h   (6,ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_h_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: c       (6,ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: c_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: vx_r    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vx_r_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vy_r    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vy_r_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vz_r    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vz_r_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: hdiv    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: hdiv_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: cp      (ADM_gall,   ADM_kall,ADM_lall,   ADM_AI:ADM_AJ,GRD_XDIR:GRD_ZDIR)
    real(8) :: cp_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_XDIR:GRD_ZDIR)

    real(8), parameter :: b1 = 0.D0
    real(8), parameter :: b2 = 1.D0
    real(8), parameter :: b3 = 1.D0 - (b1+b2)

    integer :: nstart, nend
    integer :: g, k, l, n, nq

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 1st
    !---------------------------------------------------------------------------
    do l = 1, ADM_lall

       do k = 1, ADM_kall
          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog_in(g,k,l)
             d(g,k,l) = b1 * frhog(g,k,l) * dt * rrhog(g,k,l)
          enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
             flx_v(g,k,l) = ( ( ( GRD_afac(k) * VMTR_RGSGAM2(g,k  ,l) * rhogvx_mean(g,k  ,l) &
                                + GRD_bfac(k) * VMTR_RGSGAM2(g,k-1,l) * rhogvx_mean(g,k-1,l) &
                                ) * 0.5D0 * VMTR_GSGAMH(g,k,l) * VMTR_GZXH(g,k,l)            &
                              + ( GRD_afac(k) * VMTR_RGSGAM2(g,k  ,l) * rhogvy_mean(g,k  ,l) &
                                + GRD_bfac(k) * VMTR_RGSGAM2(g,k-1,l) * rhogvy_mean(g,k-1,l) &
                                ) * 0.5D0 * VMTR_GSGAMH(g,k,l) * VMTR_GZYH(g,k,l)            & 
                              + ( GRD_afac(k) * VMTR_RGSGAM2(g,k  ,l) * rhogvz_mean(g,k  ,l) &
                                + GRD_bfac(k) * VMTR_RGSGAM2(g,k-1,l) * rhogvz_mean(g,k-1,l) &
                                ) * 0.5D0 * VMTR_GSGAMH(g,k,l) * VMTR_GZZH(g,k,l)            &
                              ) + rhogw_mean(g,k,l) * VMTR_RGSH(g,k,l) &
                            ) * (0.5D0*dt)
       enddo
       enddo

       do g = 1, ADM_gall
          flx_v(g,ADM_kmin,  l) = 0.D0
          flx_v(g,ADM_kmax+1,l) = 0.D0
       enddo

       !---- Courant numbers at cell boundary
       do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             ck(g,k,l,1) = -flx_v(g,k  ,l) * rrhog(g,k,l) * GRD_rdgz(k)
             ck(g,k,l,2) =  flx_v(g,k+1,l) * rrhog(g,k,l) * GRD_rdgz(k)
          enddo
       enddo

       do g = 1, ADM_gall
          ck(g,ADM_kmin-1,l,1) = 0.D0 
          ck(g,ADM_kmin-1,l,2) = 0.D0
          ck(g,ADM_kmax+1,l,1) = 0.D0
          ck(g,ADM_kmax+1,l,2) = 0.D0
       enddo

    enddo ! l LOOP

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_in_pl(g,k,l)
                d_pl(g,k,l) = b1 * frhog_pl(g,k,l) * dt * rrhog_pl(g,k,l)
             enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
                flx_v_pl(g,k,l) = ( ( ( GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * rhogvx_mean_pl(g,k  ,l) &
                                      + GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * rhogvx_mean_pl(g,k-1,l) &
                                      ) * 0.5D0 * VMTR_GSGAMH_pl(g,k,l) * VMTR_GZXH_pl(g,k,l)            &
                                    + ( GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * rhogvy_mean_pl(g,k  ,l) &
                                      + GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * rhogvy_mean_pl(g,k-1,l) &
                                      ) * 0.5D0 * VMTR_GSGAMH_pl(g,k,l) * VMTR_GZYH_pl(g,k,l)            & 
                                    + ( GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * rhogvz_mean_pl(g,k  ,l) &
                                      + GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * rhogvz_mean_pl(g,k-1,l) &
                                      ) * 0.5D0 * VMTR_GSGAMH_pl(g,k,l) * VMTR_GZZH_pl(g,k,l)            &
                                    ) + rhogw_mean_pl(g,k,l) * VMTR_RGSH_pl(g,k,l)                       &
                                  ) * (0.5d0*dt)
          enddo
          enddo

          do g = 1, ADM_gall_pl
             flx_v_pl(g,ADM_kmin,  l) = 0.D0
             flx_v_pl(g,ADM_kmax+1,l) = 0.D0
          enddo

          !---- Courant numbers at cell boundary
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                ck_pl(g,k,l,1) = -flx_v_pl(g,k  ,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
                ck_pl(g,k,l,2) =  flx_v_pl(g,k+1,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
             enddo
          enddo

          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmin-1,l,1) = 0.D0
             ck_pl(g,ADM_kmin-1,l,2) = 0.D0
             ck_pl(g,ADM_kmax+1,l,1) = 0.D0
             ck_pl(g,ADM_kmax+1,l,2) = 0.D0
          enddo

       enddo
    endif

    !--- basic scheme ( 2nd-order centered difference )
    do nq = 1, nqmax

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall
                q_h(g,k,l) = 0.5D0 * ( GRD_afac(k) * q(g,k,  l) &
                                     + GRD_bfac(k) * q(g,k-1,l) )
             enddo
          enddo

          do g = 1, ADM_gall
             q_h(g,ADM_kmin-1,l) = 0.D0
          enddo
       enddo


       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo

             do k = ADM_kmin, ADM_kmax+1
                do g = 1, ADM_gall_pl
                   q_h_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * q_pl(g,k,  l) &
                                           + GRD_bfac(k) * q_pl(g,k-1,l) )
                enddo
             enddo

             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin-1,l) = 0.D0
             enddo
          enddo
       endif

       ! [mod] 20130613 R.Yoshida
       if (thubern_lim) call advlim_thuburn_v( q_h, q_h_pl, & !--- [INOUT]
                              q,   q_pl,   & !--- [IN]
                              ck,  ck_pl,  & !--- [IN]
                              d,   d_pl    ) !--- [IN]

       !--- update rhogq
       do l = 1, ADM_lall
          do g = 1, ADM_gall
             q_h(g,ADM_kmin  ,l) = 0.D0
             q_h(g,ADM_kmax+1,l) = 0.D0
          enddo

          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) &
                                - ( flx_v(g,k+1,l) * q_h(g,k+1,l) &
                                  - flx_v(g,k,  l) * q_h(g,k,  l) ) * GRD_rdgz(k)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin  ,l) = 0.D0
                q_h_pl(g,ADM_kmax+1,l) = 0.D0
             enddo

             do k = ADM_kmin, ADM_kmax
                do g = 1, ADM_gall_pl
                   rhogq_pl(g,k,l,nq) = rhogq_pl(g,k,l,nq) &
                                      - ( flx_v_pl(g,k+1,l) * q_h_pl(g,k+1,l) &
                                        - flx_v_pl(g,k,  l) * q_h_pl(g,k,  l) ) * GRD_rdgz(k)
                enddo
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

    !--- update rhog
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             rhog(g,k,l) = rhog_in(g,k,l)                                  &
                         - ( flx_v(g,k+1,l) - flx_v(g,k,l) ) * GRD_rdgz(k) &
                         + b1 * frhog(g,k,l) * dt
          enddo
       enddo

       do k = ADM_kmin, ADM_kmax
          rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
          rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
       enddo

       do g = 1, ADM_gall
          rhog(g,ADM_kmin-1,l) = rhog_in(g,ADM_kmin,l)
          rhog(g,ADM_kmax+1,l) = rhog_in(g,ADM_kmax,l)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                rhog_pl(g,k,l) = rhog_in_pl(g,k,l)                                     &
                               - ( flx_v_pl(g,k+1,l) - flx_v_pl(g,k,l) ) * GRD_rdgz(k) &
                               + b1 * frhog_pl(g,k,l) * dt
             enddo
          enddo

          do g = 1, ADM_gall_pl
             rhog_pl(g,ADM_kmin-1,l) = rhog_in_pl(g,ADM_kmin,l)
             rhog_pl(g,ADM_kmax+1,l) = rhog_in_pl(g,ADM_kmax,l)
          enddo
       enddo
    endif


    !---------------------------------------------------------------------------
    ! Horizontal advection by MIURA 2004 scheme
    !---------------------------------------------------------------------------
    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do g = 1, ADM_gall
             vx_r(g,k,l) = rhogvx_mean(g,k,l) * VMTR_RGAM(g,k,l)
             vy_r(g,k,l) = rhogvy_mean(g,k,l) * VMTR_RGAM(g,k,l) 
             vz_r(g,k,l) = rhogvz_mean(g,k,l) * VMTR_RGAM(g,k,l) 
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                vx_r_pl(g,k,l) = rhogvx_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)  
                vy_r_pl(g,k,l) = rhogvy_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)   
                vz_r_pl(g,k,l) = rhogvz_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)   
             enddo
          enddo
       enddo
    endif

    call oprt_divergence2_prep( flx_h,     flx_h_pl,     & !--- [OUT]
                                cp,        cp_pl,        & !--- [OUT]
                                vx_r,      vx_r_pl,      & !--- [IN]
                                vy_r,      vy_r_pl,      & !--- [IN]
                                vz_r,      vz_r_pl,      & !--- [IN]
                                rhog_mean, rhog_mean_pl, & !--- [IN]
                                dt                       ) !--- [IN]

    do l = 1, ADM_lall
    do k = 1, ADM_kall

          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog(g,k,l)
             d(g,k,l) = b2 * frhog(g,k,l) * dt * rrhog(g,k,l)
             !--- Courant number
             c(1,g,k,l) = flx_h(1,g,k,l) * rrhog(g,k,l)
             c(2,g,k,l) = flx_h(2,g,k,l) * rrhog(g,k,l)
             c(3,g,k,l) = flx_h(3,g,k,l) * rrhog(g,k,l)
             c(4,g,k,l) = flx_h(4,g,k,l) * rrhog(g,k,l)
             c(5,g,k,l) = flx_h(5,g,k,l) * rrhog(g,k,l)
             c(6,g,k,l) = flx_h(6,g,k,l) * rrhog(g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall

             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_pl(g,k,l)
             enddo

             do n = ADM_GMIN_PL, ADM_gmax_pl
                c_pl(n,k,l) = flx_h_pl(n,k,l) * rrhog_pl(ADM_GSLF_PL,k,l)
             enddo

             d_pl(ADM_GSLF_PL,k,l) = b2 * frhog_pl(ADM_GSLF_PL,k,l) * dt * rrhog_pl(ADM_GSLF_PL,k,l)
          enddo
       enddo
    endif

    do nq = 1, nqmax

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo
          enddo
       endif


       call oprt_divergence2( hdiv,  hdiv_pl,  & !--- [OUT]
                                  q,     q_pl,     & !--- [IN]
                                  flx_h, flx_h_pl, & !--- [IN]
                                  c,     c_pl,     & !--- [IN]
                                  cp,    cp_pl,    & !--- [IN]
                                  d,     d_pl,     & !--- [IN]
                                  dt               ) !--- [IN]

       !--- update rhogq
       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) - hdiv(g,k,l) * dt
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   rhogq_pl(g,k,l,nq) = rhogq_pl(g,k,l,nq) - hdiv_pl(g,k,l) * dt
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
                                        + flx_h(6,n,k,l) &
                                        ) + b2 * frhog(n,k,l) * dt
          enddo
       enddo

       do k = 1, ADM_kall
       rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
       rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_GSLF_PL

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
             rhog_pl(n,k,l)= rhog_pl(n,k,l) - ( flx_h_pl(ADM_GMIN_PL  ,k,l) &
                                              + flx_h_pl(ADM_GMIN_PL+1,k,l) &
                                              + flx_h_pl(ADM_GMIN_PL+2,k,l) &
                                              + flx_h_pl(ADM_GMIN_PL+3,k,l) &
                                              + flx_h_pl(ADM_GMIN_PL+4,k,l) &
                                              ) + b2 * frhog_pl(n,k,l) * dt
       enddo
       enddo
    endif

    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 2nd
    !---------------------------------------------------------------------------
    do l = 1, ADM_lall

       do k = 1, ADM_kall
          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog(g,k,l)
          enddo
       enddo

       !---- Courant numbers at cell boundary
       do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             ck(g,k,l,1) = -flx_v(g,k  ,l) * rrhog(g,k,l) * GRD_rdgz(k)
             ck(g,k,l,2) =  flx_v(g,k+1,l) * rrhog(g,k,l) * GRD_rdgz(k)
          enddo
       enddo

       do g = 1, ADM_gall
          ck(g,ADM_kmin-1,l,1) = 0.D0 
          ck(g,ADM_kmin-1,l,2) = 0.D0
          ck(g,ADM_kmax+1,l,1) = 0.D0
          ck(g,ADM_kmax+1,l,2) = 0.D0
       enddo

       do k = 1, ADM_kall
          do g = 1, ADM_gall
             d(g,k,l) = b3 * frhog(g,k,l) * dt * rrhog(g,k,l)
          enddo
       enddo

    enddo ! l LOOP

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_pl(g,k,l)
             enddo
          enddo

          !---- Courant numbers at cell boundary
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                ck_pl(g,k,l,1) = -flx_v_pl(g,k  ,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
                ck_pl(g,k,l,2) =  flx_v_pl(g,k+1,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
             enddo
          enddo

          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmin-1,l,1) = 0.D0
             ck_pl(g,ADM_kmin-1,l,2) = 0.D0
             ck_pl(g,ADM_kmax+1,l,1) = 0.D0
             ck_pl(g,ADM_kmax+1,l,2) = 0.D0
          enddo

          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                d_pl(g,k,l) = b3 * frhog_pl(g,k,l) * dt * rrhog_pl(g,k,l)
             enddo
          enddo

       enddo
    endif

    !--- basic scheme ( 2nd-order centered difference )
    do nq = 1, nqmax

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall
                q_h(g,k,l) = 0.5D0 * ( GRD_afac(k) * q(g,k,  l) &
                                     + GRD_bfac(k) * q(g,k-1,l) )
             enddo
          enddo

          do g = 1, ADM_gall
             q_h(g,ADM_kmin-1,l) = 0.D0
          enddo
       enddo


       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo

             do k = ADM_kmin, ADM_kmax+1
                do g = 1, ADM_gall_pl
                   q_h_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * q_pl(g,k,  l) &
                                           + GRD_bfac(k) * q_pl(g,k-1,l) )
                enddo
             enddo

             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin-1,l) = 0.D0
             enddo
          enddo
       endif

       ! [mod] 20130613 R.Yoshida
       if (thubern_lim) call advlim_thuburn_v( q_h, q_h_pl, & !--- [INOUT]
                              q,   q_pl,   & !--- [IN]
                              ck,  ck_pl,  & !--- [IN]
                              d,   d_pl    ) !--- [IN]

       !--- update rhogq
       do l = 1, ADM_lall
          do g = 1, ADM_gall
             q_h(g,ADM_kmin  ,l) = 0.D0
             q_h(g,ADM_kmax+1,l) = 0.D0
          enddo

          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) &
                                - ( flx_v(g,k+1,l) * q_h(g,k+1,l) &
                                  - flx_v(g,k,  l) * q_h(g,k,  l) ) * GRD_rdgz(k)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin  ,l) = 0.D0
                q_h_pl(g,ADM_kmax+1,l) = 0.D0
             enddo

             do k = ADM_kmin, ADM_kmax
                do g = 1, ADM_gall_pl
                   rhogq_pl(g,k,l,nq) = rhogq_pl(g,k,l,nq) &
                                      - ( flx_v_pl(g,k+1,l) * q_h_pl(g,k+1,l) &
                                        - flx_v_pl(g,k,  l) * q_h_pl(g,k,  l) ) * GRD_rdgz(k)
                enddo
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

    return
  end subroutine src_update_tracer

  !-----------------------------------------------------------------------------
  subroutine advlim_thuburn_v( &
       q_h, q_h_pl, &
       q,   q_pl,   &
       ck,  ck_pl,  &
       d,   d_pl    )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
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
    real(8), intent(in)    :: ck    (ADM_gall,   ADM_kall,ADM_lall   ,2)
    real(8), intent(in)    :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(8), intent(in)    :: d     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: qin_min    (ADM_gall,   ADM_kall,2)
    real(8) :: qin_min_pl (ADM_gall_pl,ADM_kall,2)
    real(8) :: qin_max    (ADM_gall,   ADM_kall,2)
    real(8) :: qin_max_pl (ADM_gall_pl,ADM_kall,2)

    real(8) :: qout_min   (ADM_gall,   ADM_kall)
    real(8) :: qout_min_pl(ADM_gall_pl,ADM_kall)
    real(8) :: qout_max   (ADM_gall,   ADM_kall)
    real(8) :: qout_max_pl(ADM_gall_pl,ADM_kall)

    real(8) :: q_mp1_min
    real(8) :: q_mp1_min_pl
    real(8) :: q_mp1_max
    real(8) :: q_mp1_max_pl

    real(8) :: c_in
    real(8) :: c_in_pl
    real(8) :: c_out
    real(8) :: c_out_pl

    real(8) :: c_qin_min
    real(8) :: c_qin_min_pl
    real(8) :: c_qin_max
    real(8) :: c_qin_max_pl

    real(8) :: mask, mask1, mask2
    real(8) :: tmp

    integer :: n, k, l
    !---------------------------------------------------------------------------

    !--- inflow & outflow limiter
    do l = 1, ADM_lall

       do k = 1, ADM_kall
          do n = 1, ADM_gall
             qin_min(n,k,1) = CNST_MAX_REAL
             qin_max(n,k,1) =-CNST_MAX_REAL
          enddo
          do n = 1, ADM_gall
             qin_min(n,k,2) = CNST_MAX_REAL
             qin_max(n,k,2) =-CNST_MAX_REAL
          enddo
       enddo

       do k = ADM_kmin, ADM_kmax+1
          do n = 1,ADM_gall
             if ( ck(n,k,l,1) > 0.D0 ) then
                qin_min(n,k-1,2) = min( q(n,k,l), q(n,k-1,l) )
                qin_max(n,k-1,2) = max( q(n,k,l), q(n,k-1,l) )
             else
                qin_min(n,k,  1) = min( q(n,k,l), q(n,k-1,l) )
                qin_max(n,k,  1) = max( q(n,k,l), q(n,k-1,l) )
             endif
          enddo
       enddo

       do k = 1, ADM_kall
          do n = 1, ADM_gall

             q_mp1_min = min( qin_min(n,k,1), qin_min(n,k,2) )
             if( q_mp1_min ==  CNST_MAX_REAL ) q_mp1_min = q(n,k,l)
             q_mp1_min = max( 0.D0, q_mp1_min )

             q_mp1_max = max( qin_max(n,k,1), qin_max(n,k,2) )
             if( q_mp1_max == -CNST_MAX_REAL ) q_mp1_max = q(n,k,l)

             mask1 = 0.5D0 - sign( 0.5D0, ck(n,k,l,1) )
             mask2 = 0.5D0 - sign( 0.5D0, ck(n,k,l,2) )

             c_in  = (        mask1 ) * ck(n,k,l,1) &
                   + (        mask2 ) * ck(n,k,l,2)
             c_out = ( 1.D0 - mask1 ) * ck(n,k,l,1) &
                   + ( 1.D0 - mask2 ) * ck(n,k,l,2)

             c_qin_max = ( mask1 ) * ( ck(n,k,l,1) * qin_max(n,k,1) ) &
                       + ( mask2 ) * ( ck(n,k,l,2) * qin_max(n,k,2) )
             c_qin_min = ( mask1 ) * ( ck(n,k,l,1) * qin_min(n,k,1) ) &
                       + ( mask2 ) * ( ck(n,k,l,2) * qin_min(n,k,2) )

             if ( abs(c_out) <= CNST_EPS_ZERO ) then
                qout_min(n,k) = q(n,k,l)
                qout_max(n,k) = q(n,k,l)
             else
                qout_min(n,k) = ( q(n,k,l)                                  &
                                 - c_qin_max                                &
                                 - q_mp1_max * ( 1.D0-c_in-c_out+d(n,k,l) ) &
                                 ) / c_out
                qout_max(n,k) = ( q(n,k,l)                                  &
                                 - c_qin_min                                &
                                 - q_mp1_min * ( 1.D0-c_in-c_out+d(n,k,l) ) &
                                 ) / c_out
             endif

          enddo
       enddo

       !--- apply limiter
       do k = ADM_kmin, ADM_kmax+1
          do n = 1, ADM_gall
             mask = 0.5D0 - sign( 0.5D0, ck(n,k,l,1) )

             tmp = (        mask ) * min( max( q_h(n,k,l), qin_min(n,k  ,1) ), qin_max(n,k  ,1) ) &
                 + ( 1.D0 - mask ) * min( max( q_h(n,k,l), qin_min(n,k-1,2) ), qin_max(n,k-1,2) )

             q_h(n,k,l) = (        mask ) * max( min( tmp, qout_max(n,k-1) ), qout_min(n,k-1) ) &
                        + ( 1.D0 - mask ) * max( min( tmp, qout_max(n,k  ) ), qout_min(n,k  ) )
          enddo
       enddo

    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          do k = 1, ADM_kall
             do n = 1, ADM_gall_pl
                qin_min_pl(n,k,1) =  CNST_MAX_REAL
                qin_max_pl(n,k,1) = -CNST_MAX_REAL
             enddo
             do n = 1, ADM_gall_pl
                qin_min_pl(n,k,2) =  CNST_MAX_REAL
                qin_max_pl(n,k,2) = -CNST_MAX_REAL
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
             do n = 1, ADM_gall_pl
                if ( ck_pl(n,k,l,1) > 0.D0 ) then
                   qin_min_pl(n,k-1,2) = min( q_pl(n,k,l), q_pl(n,k-1,l) )
                   qin_max_pl(n,k-1,2) = max( q_pl(n,k,l), q_pl(n,k-1,l) )
                else
                   qin_min_pl(n,k,  1) = min( q_pl(n,k,l), q_pl(n,k-1,l) )
                   qin_max_pl(n,k,  1) = max( q_pl(n,k,l), q_pl(n,k-1,l) )
                endif
             enddo
          enddo

          do k = 1, ADM_kall
             do n = 1, ADM_gall_pl

                q_mp1_min_pl = min( qin_min_pl(n,k,1), qin_min_pl(n,k,2) )
                if( q_mp1_min_pl ==  CNST_MAX_REAL ) q_mp1_min_pl = q_pl(n,k,l)
                q_mp1_min_pl = max( 0.D0, q_mp1_min_pl)

                q_mp1_max_pl = max( qin_max_pl(n,k,1), qin_max_pl(n,k,2) )
                if( q_mp1_max_pl == -CNST_MAX_REAL ) q_mp1_max_pl = q_pl(n,k,l)

                mask1 = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,1) )
                mask2 = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,2) )

                c_in_pl  = (        mask1 ) * ck_pl(n,k,l,1) &
                         + (        mask2 ) * ck_pl(n,k,l,2)
                c_out_pl = ( 1.D0 - mask1 ) * ck_pl(n,k,l,1) &
                         + ( 1.D0 - mask2 ) * ck_pl(n,k,l,2)

                c_qin_max_pl = ( mask1 ) * ( ck_pl(n,k,l,1) * qin_max_pl(n,k,1) ) &
                             + ( mask2 ) * ( ck_pl(n,k,l,2) * qin_max_pl(n,k,2) )
                c_qin_min_pl = ( mask1 ) * ( ck_pl(n,k,l,1) * qin_min_pl(n,k,1) ) &
                             + ( mask2 ) * ( ck_pl(n,k,l,2) * qin_min_pl(n,k,2) )

                if ( abs(c_out_pl) < CNST_EPS_ZERO ) then
                   qout_min_pl(n,k) = q_pl(n,k,l)
                   qout_max_pl(n,k) = q_pl(n,k,l)
                else
                   qout_min_pl(n,k) = ( q_pl(n,k,l)                                          &
                                      - c_qin_max_pl                                         &
                                      - q_mp1_max_pl * ( 1.D0-c_in_pl-c_out_pl+d_pl(n,k,l) ) &
                                      ) / c_out_pl
                   qout_max_pl(n,k) = ( q_pl(n,k,l)                                          &
                                      - c_qin_min_pl                                         &
                                      - q_mp1_min_pl * ( 1.D0-c_in_pl-c_out_pl+d_pl(n,k,l) ) &
                                      ) / c_out_pl
                endif

             enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
             do n = 1, ADM_gall_pl
                mask = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,1) )

                tmp = (        mask ) * min( max( q_h_pl(n,k,l), qin_min_pl(n,k  ,1) ), qin_max_pl(n,k  ,1) ) &
                    + ( 1.D0 - mask ) * min( max( q_h_pl(n,k,l), qin_min_pl(n,k-1,2) ), qin_max_pl(n,k-1,2) )

                q_h_pl(n,k,l) = (        mask ) * max( min( tmp, qout_max_pl(n,k-1) ), qout_min_pl(n,k-1) ) &
                              + ( 1.D0 - mask ) * max( min( tmp, qout_max_pl(n,k  ) ), qout_min_pl(n,k  ) )
             enddo
          enddo

       enddo
    endif

    return
  end subroutine advlim_thuburn_v

  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence2_prep( &
       flx_h,   flx_h_pl,  &
       cp,      cp_pl,     &
       rhovx,  rhovx_pl, &
       rhovy,  rhovy_pl, &
       rhovz,  rhovz_pl, &
       rho,    rho_pl,   &
       dt                  &
       )
    !
    !--- Miura(2004)'s scheme : path1
    !
    use mod_adm, only :   &
         !--- public parameters
         ADM_W,           &
         ADM_prc_pl,      &
         ADM_TI,          &
         ADM_TJ,          &
         ADM_AI,          &
         ADM_AIJ,         &
         ADM_AJ,          &
         ADM_KNONE,       &
         ADM_lall_pl,     &
         ADM_gmin_pl,     &
         ADM_gslf_pl,     &
         ADM_gall_pl,     &
         ADM_gmax_pl,     &
         !--- public variables
         ADM_prc_me,      &
         ADM_prc_tab,     &
         ADM_rgn_vnum,    &
         ADM_gall_1d,     &
         ADM_kall,        &
         ADM_lall,        &
         ADM_kmin,        &
         ADM_kmax,        &
         ADM_gmin,        &
         ADM_gmax,        &
         ADM_gall,        &
         ADM_ImoJmo_nmax, &  ! Y.Niwa add 080130
         ADM_ImoJmo,      &  ! Y.Niwa add 080130
         ADM_GIoJo,       &  ! Y.Niwa add 080130
         ADM_VMISS           ! Y.Niwa add 080130
    use mod_gmtr, only :  &
         !--- public parameters
         GMTR_T_W1,       &
         GMTR_T_W2,       &
         GMTR_T_W3,       &
         GMTR_T_rarea,    &
         GMTR_A_hnx,      &
         GMTR_A_hny,      &
         GMTR_A_hnz,      &
         GMTR_A_tnx,      &
         GMTR_A_tny,      &
         GMTR_A_tnz,      &
         GMTR_A_tn2x,     &
         GMTR_A_tn2y,     &
         GMTR_A_tn2z,     &
         GMTR_P_rarea,    &
         !--- public variables
         GMTR_T_var,      &
         GMTR_T_var_pl,   &
         GMTR_P_var,      &
         GMTR_P_var_pl,   &
         GMTR_A_var,      &
         GMTR_A_var_pl
    use mod_grd, only :    &
         GRD_x, GRD_x_pl,  &
         GRD_xt, GRD_xt_pl,&
         GRD_XDIR,         &  
         GRD_YDIR,         &
         GRD_ZDIR
    use mod_comm, only : &
         COMM_data_transfer
    use mod_cnst, only : &
         CNST_MAX_REAL,  &
         CNST_EPS_ZERO, &
         CNST_UNDEF
    implicit none

    real(8), intent(out) :: flx_h(6,ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: flx_h_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: cp(ADM_gall,ADM_kall,ADM_lall,ADM_AI:ADM_AJ,GRD_XDIR:GRD_ZDIR)
    real(8), intent(out) :: cp_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,GRD_XDIR:GRD_ZDIR)
    real(8), intent(in)  :: rhovx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: rhovx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in)  :: rhovy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: rhovy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in)  :: rhovz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: rhovz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in)  :: rho(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: rho_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in)  :: dt

    real(8) :: rhovxt(ADM_gall,ADM_kall,ADM_lall,ADM_TI:ADM_TJ)
    real(8) :: rhovxt_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rhovyt(ADM_gall,ADM_kall,ADM_lall,ADM_TI:ADM_TJ)
    real(8) :: rhovyt_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rhovzt(ADM_gall,ADM_kall,ADM_lall,ADM_TI:ADM_TJ)
    real(8) :: rhovzt_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rhot(ADM_gall,ADM_kall,ADM_lall,ADM_TI:ADM_TJ)
    real(8) :: rhot_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    real(8)  :: ccc
    real(8)  :: ccc_pl

    real(8)  :: rpx(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rpx_pl(ADM_gall_pl)
    real(8)  :: rpy(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rpy_pl(ADM_gall_pl)
    real(8)  :: rpz(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rpz_pl(ADM_gall_pl)
    !
    real(8)  :: rhovxa(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rhovxa_pl(ADM_GALL_PL)
    real(8)  :: rhovya(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rhovya_pl(ADM_GALL_PL)
    real(8)  :: rhovza(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rhovza_pl(ADM_GALL_PL)
    real(8)  :: rhoa(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rhoa_pl(ADM_GALL_PL)
    !
    integer :: l,n,k, k0, ij
    integer :: rgnid
    !
    integer :: np1(ADM_gall_pl)
    integer :: nm1(ADM_gall_pl)
    !
    integer :: nstart,nend
    logical, save :: first = .true.  ! Y.Niwa add 080130
    integer :: t  ! Y.Niwa add 080130
    !
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    if (first) then
       first = .false.

       allocate( local_t_var   (ADM_gall   ,k0,ADM_lall   ,ADM_TI:ADM_TJ,GMTR_T_W1:GMTR_T_W3) )
       allocate( local_t_var_pl(ADM_gall_pl,k0,ADM_lall_pl,              GMTR_T_W1:GMTR_T_W3) )
       local_t_var = CNST_UNDEF

       do l = 1, ADM_lall
       do t = ADM_TI, ADM_TJ
       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          local_t_var(ij,k0,l,t,GMTR_T_W1) = GMTR_t_var(ij,k0,l,t,GMTR_T_W1)
          local_t_var(ij,k0,l,t,GMTR_T_W2) = GMTR_t_var(ij,k0,l,t,GMTR_T_W2)
          local_t_var(ij,k0,l,t,GMTR_T_W3) = GMTR_t_var(ij,k0,l,t,GMTR_T_W3)
       enddo
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do n = ADM_gmin_pl, ADM_gmax_pl
             local_t_var_pl(n,k0,l,GMTR_T_W1) = GMTR_t_var_pl(n,k0,l,GMTR_T_W1)
             local_t_var_pl(n,k0,l,GMTR_T_W2) = GMTR_t_var_pl(n,k0,l,GMTR_T_W2)
             local_t_var_pl(n,k0,l,GMTR_T_W3) = GMTR_t_var_pl(n,k0,l,GMTR_T_W3)
          enddo
          enddo
       endif
    endif

    flx_h   (:,:,:,:)   = 0.D0
    flx_h_pl(:,:,:)     = 0.D0
    cp      (:,:,:,:,:) = 0.D0
    cp_pl   (:,:,:,:)   = 0.D0

    nm1(ADM_gmin_pl) = ADM_gmax_pl
    do n = ADM_gmin_pl+1, ADM_gmax_pl
       nm1(n) = n - 1
    enddo

    do n = ADM_gmin_pl, ADM_gmax_pl-1
       np1(n) = n + 1
    enddo
    np1(ADM_gmax_pl) = ADM_gmin_pl

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do k = 1, ADM_kall

          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart, nend
             rhovxt(n,k,l,ADM_TI) = GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W1) * rhovx(n              ,k,l) &
                                  + GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W2) * rhovx(n+1            ,k,l) &
                                  + GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W3) * rhovx(n+1+ADM_gall_1d,k,l)
             rhovxt(n,k,l,ADM_TJ) = GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W1) * rhovx(n              ,k,l) &
                                  + GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W2) * rhovx(n+1+ADM_gall_1d,k,l) &
                                  + GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W3) * rhovx(n+ADM_gall_1d,  k,l)
             rhovyt(n,k,l,ADM_TI) =GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W1) *rhovy(n,k,l) &
                  +GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W2) *rhovy(n+1,k,l) &
                  +GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W3) *rhovy(n+1+ADM_gall_1d,k,l)
             rhovyt(n,k,l,ADM_TJ) =GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W1) *rhovy(n,k,l) &
                  +GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W2) *rhovy(n+1+ADM_gall_1d,k,l) &
                  +GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W3) *rhovy(n+ADM_gall_1d,k,l)
             rhovzt(n,k,l,ADM_TI) =GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W1) *rhovz(n,k,l) &
                  +GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W2) *rhovz(n+1,k,l) &
                  +GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W3) *rhovz(n+1+ADM_gall_1d,k,l)
             rhovzt(n,k,l,ADM_TJ) =GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W1) *rhovz(n,k,l) &
                  +GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W2) *rhovz(n+1+ADM_gall_1d,k,l) &
                  +GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W3) *rhovz(n+ADM_gall_1d,k,l)
             rhot(n,k,l,ADM_TI) =local_T_var(n,k0,l,ADM_TI,GMTR_T_W1) *rho(n,k,l) &
                  +local_T_var(n,k0,l,ADM_TI,GMTR_T_W2) *rho(n+1,k,l) &
                  +local_T_var(n,k0,l,ADM_TI,GMTR_T_W3) *rho(n+1+ADM_gall_1d,k,l)
!                 =GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W1)&  ! A.Noda 110728
!                 *rho(n,k,l)                         &
!                 +GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W2)&
!                 *rho(n+1,k,l)                       &
!                 +GMTR_T_var(n,k0,l,ADM_TI,GMTR_T_W3)&
!                 *rho(n+1+ADM_gall_1d,k,l)
             rhot(n,k,l,ADM_TJ)                            &
                 =local_T_var(n,k0,l,ADM_TJ,GMTR_T_W1)&
                  *rho(n,k,l)                         &
                  +local_T_var(n,k0,l,ADM_TJ,GMTR_T_W2)&
                  *rho(n+1+ADM_gall_1d,k,l)                     &
                  +local_T_var(n,k0,l,ADM_TJ,GMTR_T_W3)&
                  *rho(n+ADM_gall_1d,k,l)
!!$                  =GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W1)&   ! Y.Niwa 080130
!!$                  *rho(n,k,l)                         &
!!$                  +GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W2)&
!!$                  *rho(n+1+ADM_gall_1d,k,l)                     &
!!$                  +GMTR_T_var(n,k0,l,ADM_TJ,GMTR_T_W3)&
!!$                  *rho(n+ADM_gall_1d,k,l)
          enddo
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             rhovxt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =rhovxt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             rhovyt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =rhovyt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             rhovzt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =rhovzt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             rhot(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =rhot(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
          !
       enddo
       !
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
       !
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n=ADM_gmin_pl,ADM_gmax_pl
                rhovxt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,k0,l,GMTR_T_W1)&
                     *rhovx_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,k0,l,GMTR_T_W2)&
                     *rhovx_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,k0,l,GMTR_T_W3)&
                     *rhovx_pl(np1(n),k,l)
                rhovyt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,k0,l,GMTR_T_W1)&
                     *rhovy_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,k0,l,GMTR_T_W2)&
                     *rhovy_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,k0,l,GMTR_T_W3)&
                     *rhovy_pl(np1(n),k,l)
                rhovzt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,k0,l,GMTR_T_W1)&
                     *rhovz_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,k0,l,GMTR_T_W2)&
                     *rhovz_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,k0,l,GMTR_T_W3)&
                     *rhovz_pl(np1(n),k,l)
                rhot_pl(n,k,l)                            &
                     =local_T_var_pl(n,k0,l,GMTR_T_W1)&
                     *rho_pl(ADM_gslf_pl,k,l)              &
                     +local_T_var_pl(n,k0,l,GMTR_T_W2)&
                     *rho_pl(n,k,l)                        &
                     +local_T_var_pl(n,k0,l,GMTR_T_W3)&
                     *rho_pl(np1(n),k,l)
!!$                     =GMTR_T_var_pl(n,k0,l,GMTR_T_W1)&  ! Y.Niwa 080130
!!$                     *rho_pl(ADM_gslf_pl,k,l)              &
!!$                     +GMTR_T_var_pl(n,k0,l,GMTR_T_W2)&
!!$                     *rho_pl(n,k,l)                        &
!!$                     +GMTR_T_var_pl(n,k0,l,GMTR_T_W3)&
!!$                     *rho_pl(np1(n),k,l)
             enddo
          enddo
       enddo
    end if

    !--- calculation of courant number
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)

       do k = 1, ADM_kall
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax,  ADM_gmax  )
          do n = nstart,nend
             rpx(n,ADM_AI) = 0.5D0 * ( GRD_xt(n,k0,l,ADM_TI,GRD_XDIR) + GRD_xt(n-ADM_gall_1d,k0,l,ADM_TJ,GRD_XDIR) )
             rpy(n,ADM_AI) = 0.5D0 * ( GRD_xt(n,k0,l,ADM_TI,GRD_YDIR) + GRD_xt(n-ADM_gall_1d,k0,l,ADM_TJ,GRD_YDIR) )
             rpz(n,ADM_AI) = 0.5D0 * ( GRD_xt(n,k0,l,ADM_TI,GRD_ZDIR) + GRD_xt(n-ADM_gall_1d,k0,l,ADM_TJ,GRD_ZDIR) )

             rhovxa(n,ADM_AI) = 0.5D0 * ( rhovxt(n-ADM_gall_1d,k,l,ADM_TJ) + rhovxt(n,k,l,ADM_TI) )
             rhovya(n,ADM_AI) = 0.5D0 * ( rhovyt(n-ADM_gall_1d,k,l,ADM_TJ) + rhovyt(n,k,l,ADM_TI) )
             rhovza(n,ADM_AI) = 0.5D0 * ( rhovzt(n-ADM_gall_1d,k,l,ADM_TJ) + rhovzt(n,k,l,ADM_TI) )
             rhoa  (n,ADM_AI) = 0.5D0 * ( rhot  (n-ADM_gall_1d,k,l,ADM_TJ) + rhot  (n,k,l,ADM_TI) )

             cp(n,k,l,ADM_AI,GRD_XDIR) = rpx(n,ADM_AI) - rhovxa(n,ADM_AI) / rhoa(n,ADM_AI) * dt * 0.5D0
             cp(n,k,l,ADM_AI,GRD_YDIR) = rpy(n,ADM_AI) - rhovya(n,ADM_AI) / rhoa(n,ADM_AI) * dt * 0.5D0
             cp(n,k,l,ADM_AI,GRD_ZDIR) = rpz(n,ADM_AI) - rhovza(n,ADM_AI) / rhoa(n,ADM_AI) * dt * 0.5D0

             ccc = ( rhovxa(n,ADM_AI) * GMTR_A_var(n,k0,l,ADM_AI,GMTR_A_HNX) &
                   + rhovya(n,ADM_AI) * GMTR_A_var(n,k0,l,ADM_AI,GMTR_A_HNY) &
                   + rhovza(n,ADM_AI) * GMTR_A_var(n,k0,l,ADM_AI,GMTR_A_HNZ) )

             flx_h(1,n  ,k,l) =  ccc * dt * GMTR_P_var(n  ,k0,l,GMTR_P_RAREA)
             flx_h(4,n+1,k,l) = -ccc * dt * GMTR_P_var(n+1,k0,l,GMTR_P_RAREA)
          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax,  ADM_gmax  )
          do n = nstart,nend
             rpx(n,ADM_AIJ) &
                  = 0.5D0*(GRD_xt(n,k0,l,ADM_TI,GRD_XDIR)+GRD_xt(n,k0,l,ADM_TJ,GRD_XDIR))
             rpy(n,ADM_AIJ) &
                  = 0.5D0*(GRD_xt(n,k0,l,ADM_TI,GRD_YDIR)+GRD_xt(n,k0,l,ADM_TJ,GRD_YDIR))
             rpz(n,ADM_AIJ) &
                  = 0.5D0*(GRD_xt(n,k0,l,ADM_TI,GRD_ZDIR)+GRD_xt(n,k0,l,ADM_TJ,GRD_ZDIR))
             rhovxa(n,ADM_AIJ) = (rhovxt(n,k,l,ADM_TI)+rhovxt(n,k,l,ADM_TJ))*0.5D0
             rhovya(n,ADM_AIJ) = (rhovyt(n,k,l,ADM_TI)+rhovyt(n,k,l,ADM_TJ))*0.5D0
             rhovza(n,ADM_AIJ) = (rhovzt(n,k,l,ADM_TI)+rhovzt(n,k,l,ADM_TJ))*0.5D0
             rhoa(n,ADM_AIJ) = (rhot(n,k,l,ADM_TI)+rhot(n,k,l,ADM_TJ))*0.5D0
             !
             cp(n,k,l,ADM_AIJ,GRD_XDIR) = rpx(n,ADM_AIJ) - rhovxa(n,ADM_AIJ)/rhoa(n,ADM_AIJ)*dt*0.5D0
             cp(n,k,l,ADM_AIJ,GRD_YDIR) = rpy(n,ADM_AIJ) - rhovya(n,ADM_AIJ)/rhoa(n,ADM_AIJ)*dt*0.5D0
             cp(n,k,l,ADM_AIJ,GRD_ZDIR) = rpz(n,ADM_AIJ) - rhovza(n,ADM_AIJ)/rhoa(n,ADM_AIJ)*dt*0.5D0
             !
             ccc = &
                  (rhovxa(n,ADM_AIJ)*GMTR_A_var(n,k0,l,ADM_AIJ,GMTR_A_HNX)&
                  +rhovya(n,ADM_AIJ)*GMTR_A_var(n,k0,l,ADM_AIJ,GMTR_A_HNY)&
                  +rhovza(n,ADM_AIJ)*GMTR_A_var(n,k0,l,ADM_AIJ,GMTR_A_HNZ))
             flx_h(2,n              ,k,l) = ccc*dt*GMTR_P_var(n,k0,l,GMTR_P_RAREA)
             flx_h(5,n+1+ADM_gall_1d,k,l) =-ccc*dt*GMTR_P_var(n+1+ADM_gall_1d,k0,l,GMTR_P_RAREA)
             !
          enddo

          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             rpx(n,ADM_AJ) &
                  = 0.5D0*(GRD_xt(n,k0,l,ADM_TJ,GRD_XDIR)+GRD_xt(n-1,k0,l,ADM_TI,GRD_XDIR))
             rpy(n,ADM_AJ) &
                  = 0.5D0*(GRD_xt(n,k0,l,ADM_TJ,GRD_YDIR)+GRD_xt(n-1,k0,l,ADM_TI,GRD_YDIR))
             rpz(n,ADM_AJ) &
                  = 0.5D0*(GRD_xt(n,k0,l,ADM_TJ,GRD_ZDIR)+GRD_xt(n-1,k0,l,ADM_TI,GRD_ZDIR))
             rhovxa(n,ADM_AJ)=  (rhovxt(n,k,l,ADM_TJ)+rhovxt(n-1,k,l,ADM_TI))*0.5D0
             rhovya(n,ADM_AJ)=  (rhovyt(n,k,l,ADM_TJ)+rhovyt(n-1,k,l,ADM_TI))*0.5D0
             rhovza(n,ADM_AJ)=  (rhovzt(n,k,l,ADM_TJ)+rhovzt(n-1,k,l,ADM_TI))*0.5D0
             rhoa(n,ADM_AJ)=  (rhot(n,k,l,ADM_TJ)+rhot(n-1,k,l,ADM_TI))*0.5D0
             !
             cp(n,k,l,ADM_AJ,GRD_XDIR) = rpx(n,ADM_AJ) - rhovxa(n,ADM_AJ)/rhoa(n,ADM_AJ)*dt*0.5D0
             cp(n,k,l,ADM_AJ,GRD_YDIR) = rpy(n,ADM_AJ) - rhovya(n,ADM_AJ)/rhoa(n,ADM_AJ)*dt*0.5D0
             cp(n,k,l,ADM_AJ,GRD_ZDIR) = rpz(n,ADM_AJ) - rhovza(n,ADM_AJ)/rhoa(n,ADM_AJ)*dt*0.5D0
             !
             ccc = &
                  (rhovxa(n,ADM_AJ)*GMTR_A_var(n,k0,l,ADM_AJ,GMTR_A_HNX)&
                  +rhovya(n,ADM_AJ)*GMTR_A_var(n,k0,l,ADM_AJ,GMTR_A_HNY)&
                  +rhovza(n,ADM_AJ)*GMTR_A_var(n,k0,l,ADM_AJ,GMTR_A_HNZ))
             flx_h(3,n            ,k,l) = ccc*dt*GMTR_P_var(n,k0,l,GMTR_P_RAREA)
             flx_h(6,n+ADM_gall_1d,k,l) =-ccc*dt*GMTR_P_var(n+ADM_gall_1d,k0,l,GMTR_P_RAREA)
          enddo
          !
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             flx_h(6,suf(ADM_gmin,ADM_gmin),k,l) = 0.0D0
          end if
       enddo
    enddo

    if(ADM_prc_me==ADM_prc_pl) then
       !
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n=ADM_gmin_pl,ADM_gmax_pl
                rpx_pl(n) &
                     = 0.5D0*(GRD_xt_pl(nm1(n),k0,l,GRD_XDIR)+GRD_xt_pl(n,k0,l,GRD_XDIR))
                rpy_pl(n) &
                     = 0.5D0*(GRD_xt_pl(nm1(n),k0,l,GRD_YDIR)+GRD_xt_pl(n,k0,l,GRD_YDIR))
                rpz_pl(n) &
                     = 0.5D0*(GRD_xt_pl(nm1(n),k0,l,GRD_ZDIR)+GRD_xt_pl(n,k0,l,GRD_ZDIR))
                !
                rhovxa_pl(n)=  (rhovxt_pl(nm1(n),k,l)+rhovxt_pl(n,k,l))*0.5D0
                rhovya_pl(n)=  (rhovyt_pl(nm1(n),k,l)+rhovyt_pl(n,k,l))*0.5D0
                rhovza_pl(n)=  (rhovzt_pl(nm1(n),k,l)+rhovzt_pl(n,k,l))*0.5D0
                rhoa_pl(n)=  (rhot_pl(nm1(n),k,l)+rhot_pl(n,k,l))*0.5D0
                !
                cp_pl(n,k,l,GRD_XDIR) = rpx_pl(n) - rhovxa_pl(n)/rhoa_pl(n)*dt*0.5D0
                cp_pl(n,k,l,GRD_YDIR) = rpy_pl(n) - rhovya_pl(n)/rhoa_pl(n)*dt*0.5D0
                cp_pl(n,k,l,GRD_ZDIR) = rpz_pl(n) - rhovza_pl(n)/rhoa_pl(n)*dt*0.5D0
                !
                ccc_pl = &
                     (rhovxa_pl(n)*GMTR_A_var_pl(n,k0,l,GMTR_A_HNX)&
                     +rhovya_pl(n)*GMTR_A_var_pl(n,k0,l,GMTR_A_HNY)&
                     +rhovza_pl(n)*GMTR_A_var_pl(n,k0,l,GMTR_A_HNZ))
                flx_h_pl(n,k,l) = ccc_pl*dt*GMTR_P_var_pl(ADM_gslf_pl,k0,l,GMTR_P_RAREA)
             enddo
          enddo
       enddo
    end if

    return
  end subroutine OPRT_divergence2_prep

  !-----------------------------------------------------------------------------
  ! Miura(2004)'s scheme with Thuburn(1996) limiter : path2 
  subroutine OPRT_divergence2( &
       scl,   scl_pl,   &
       s,     s_pl,     &
       flx_h, flx_h_pl, &
       c,     c_pl,     &
       cp,    cp_pl,    &
       d,     d_pl,     &
       dt,              &
       limiter,         &
       mfact            )
    use mod_adm, only: &
       ADM_prc_pl,      &
       ADM_AI,          &
       ADM_AIJ,         &
       ADM_AJ,          &
       ADM_KNONE,       &
       ADM_lall_pl,     &
       ADM_gmin_pl,     &
       ADM_gslf_pl,     &
       ADM_gall_pl,     &
       ADM_gmax_pl,     &
       ADM_prc_me,      &
       ADM_prc_tab,     &
       ADM_gall_1d,     &
       ADM_kall,        &
       ADM_lall,        &
       ADM_kmin,        &
       ADM_kmax,        &
       ADM_gmin,        &
       ADM_gmax,        &
       ADM_gall,        &
       ADM_W,           &
       ADM_rgn_vnum
    use mod_grd, only: &
       GRD_XDIR, &  
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_x,    &
       GRD_x_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_cnst, only: &
       CNST_EPS_ZERO, & 
       CNST_MAX_REAL
    use mod_oprt, only: &
       OPRT_gradient
    implicit none

    real(8), intent(out) :: scl     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: scl_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: s       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: s_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: flx_h   (6,ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: flx_h_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: c       (6,ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: c_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: cp      (ADM_gall,ADM_kall,ADM_lall,ADM_AI:ADM_AJ,GRD_XDIR:GRD_ZDIR)
    real(8), intent(in)  :: cp_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_XDIR:GRD_ZDIR)
    real(8), intent(in)  :: d       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: d_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: dt

    character(len=*), intent(in), optional :: limiter
    real(8),          intent(in), optional :: mfact

    real(8)  :: sa(ADM_AI:ADM_AJ,ADM_gall,ADM_kall,ADM_lall)
    real(8)  :: sa_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8)  :: sa_p,sa_m

    real(8)  :: sw(ADM_AI:ADM_AJ,ADM_gall)

    real(8)  :: s_in_min(6,ADM_gall,ADM_kall,ADM_lall)
    real(8)  :: s_in_min_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(8)  :: s_in_max(6,ADM_gall,ADM_kall,ADM_lall)
    real(8)  :: s_in_max_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    !
    real(8)  :: s_m1_k_min(ADM_gall)
    real(8)  :: s_m1_k_min_pl
    real(8)  :: s_m1_k_max(ADM_gall)
    real(8)  :: s_m1_k_max_pl
    real(8)  :: c_in_sum(ADM_gall)
    real(8)  :: c_in_sum_pl
    real(8)  :: c_out_sum(ADM_gall)
    real(8)  :: c_out_sum_pl
    real(8)  :: c_qin_sum_max(ADM_gall)
    real(8)  :: c_qin_sum_max_pl
    real(8)  :: c_qin_sum_min(ADM_gall)
    real(8)  :: c_qin_sum_min_pl
    !
    real(8)  :: s_m1_k_min_n
    real(8)  :: s_m1_k_max_n
    real(8)  :: c_in_sum_n
    real(8)  :: c_out_sum_n
    real(8)  :: c_qin_sum_max_n
    real(8)  :: c_qin_sum_min_n
    !
    integer :: np1(ADM_gall_pl)
    integer :: nm1(ADM_gall_pl)
    !
    real(8) :: SMALL_ZERO = 0.0d0
    !
    integer, parameter :: dsx=1
    integer, parameter :: dsy=2
    integer, parameter :: dsz=3
    integer, parameter :: s_out_k_min=4
    integer, parameter :: s_out_k_max=5
    real(8)  :: wrk(ADM_gall,ADM_kall,ADM_lall,dsx:s_out_k_max)
    real(8)  :: wrk_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,dsx:s_out_k_max)
    !
!    logical :: NON_NEG=.false.
    !
    integer :: l,n,k,m
    integer :: rgnid
    real(8) :: fact
    !
    integer :: nstart,nend
    integer :: nstart2,nstart3

    real(8) :: AI_min, AIJ_min, AJ_min
    real(8) :: AI_max, AIJ_max, AJ_max

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1, ip2jp1
    integer :: im1j, ijm1

    integer :: ierr

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    if ( present(mfact) ) then
       fact = mfact / dt
    else
       fact = 1.D0 / dt
    endif

    nm1(ADM_gmin_pl) = ADM_gmax_pl
    do n = ADM_gmin_pl+1, ADM_gmax_pl
       nm1(n) = n-1
    enddo

    do n = ADM_gmin_pl, ADM_gmax_pl-1
       np1(n) = n+1
    enddo
    np1(ADM_gmax_pl) = ADM_gmin_pl

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_gall
       scl(n,k,l) = 0.D0
    enddo
    enddo
    enddo

    do l = 1, ADM_lall_pl
    do k = 1, ADM_kall
    do n = 1, ADM_gall_pl
       scl_pl(n,k,l) = 0.D0
    enddo
    enddo
    enddo

    ! [Fix] 08/04/28 T.Mitsui, to avoid undefined reference in outflow limiter
    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do n = 1, ADM_gall
             wrk(n,k,l,s_out_k_min) = s(n,k,l)
          enddo
          do n = 1, ADM_gall
             wrk(n,k,l,s_out_k_max) = s(n,k,l)
          enddo
       enddo
    enddo

    do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          do n = 1,ADM_gall_pl
             wrk_pl(n,k,l,s_out_k_min) = s_pl(n,k,l)
          enddo
          do n = 1,ADM_gall_pl
             wrk_pl(n,k,l,s_out_k_max) = s_pl(n,k,l)
          enddo
       enddo
    enddo

    if(present(limiter)) then
       if(trim(limiter)=='NON_LIM') go to 1000
    endif

    !--- calculation of inflow limiter
!    do l = 1, ADM_lall
!    do k = 1, ADM_kall
!    do n = 1, ADM_gall
!    do m = 1, 6
!       s_in_min(m,n,k,l) = CNST_MAX_REAL
!    enddo
!    enddo
!    enddo
!    enddo
!
!    do l = 1, ADM_lall
!    do k = 1, ADM_kall
!    do n = 1, ADM_gall
!    do m = 1, 6
!       s_in_max(m,n,k,l) = -CNST_MAX_REAL
!    enddo
!    enddo
!    enddo
!    enddo

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do k = 1, ADM_kall

          do n = 1, ADM_gall
             do m = ADM_AI, ADM_AJ
                sw(m,n) = 0.D0
             enddo
          enddo

          do n = 1, ADM_gall
             if ( c(1,n,k,l) <= 0.D0 ) sw(ADM_AI, n) = 1.D0
             if ( c(2,n,k,l) <= 0.D0 ) sw(ADM_AIJ,n) = 1.D0
             if ( c(3,n,k,l) <= 0.D0 ) sw(ADM_AJ, n) = 1.D0
          enddo

          nstart = suf(ADM_gmin,ADM_gmin)
          nend   = suf(ADM_gmax,ADM_gmax)

          do n = nstart, nend
             ij     = n
             ip1j   = n+1
             ijp1   = n  +ADM_gall_1d
             ip1jp1 = n+1+ADM_gall_1d
             im1j   = n-1
             ijm1   = n  -ADM_gall_1d

             AI_min  = min( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
             AI_max  = max( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
             AIJ_min = min( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
             AIJ_max = max( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
             AJ_min  = min( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
             AJ_max  = max( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )

!             if ( c(1,n,k,l) <= 0.D0 ) then
!                s_in_min(1,ij,  k,l) = AI_min * 1.D0 + 0.D0
!                s_in_max(1,ij,  k,l) = AI_max * 1.D0 + 0.D0
!                s_in_min(4,ip1j,k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(4,ip1j,k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!             else
!                s_in_min(1,ij,  k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(1,ij,  k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_min(4,ip1j,k,l) = AI_min * 1.D0 + 0.D0
!                s_in_max(4,ip1j,k,l) = AI_max * 1.D0 + 0.D0
!             endif
!             if ( c(2,n,k,l) <= 0.D0 ) then
!                s_in_min(2,ij,    k,l) = AIJ_min * 1.D0 + 0.D0
!                s_in_max(2,ij,    k,l) = AIJ_max * 1.D0 + 0.D0
!                s_in_min(5,ip1jp1,k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(5,ip1jp1,k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!             else
!                s_in_min(2,ij,    k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(2,ij,    k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_min(5,ip1jp1,k,l) = AIJ_min * 1.D0 + 0.D0
!                s_in_max(5,ip1jp1,k,l) = AIJ_max * 1.D0 + 0.D0
!             endif
!             if ( c(3,n,k,l) <= 0.D0 ) then
!                s_in_min(3,ij,  k,l) = AJ_min * 1.D0 + 0.D0
!                s_in_max(3,ij,  k,l) = AJ_max * 1.D0 + 0.D0
!                s_in_min(6,ijp1,k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(6,ijp1,k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!             else
!                s_in_min(3,ij,  k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(3,ij,  k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_min(6,ijp1,k,l) = AJ_min * 1.D0 + 0.D0
!                s_in_max(6,ijp1,k,l) = AJ_max * 1.D0 + 0.D0
!             endif

             s_in_min(1,ij,    k,l) = (        sw(ADM_AI, n) ) * AI_min        &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * CNST_MAX_REAL
             s_in_min(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_min       &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * CNST_MAX_REAL
             s_in_min(3,ij,    k,l) = (        sw(ADM_AJ, n) ) * AJ_min        &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * CNST_MAX_REAL
             s_in_min(4,ip1j,  k,l) = (        sw(ADM_AI, n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * AI_min
             s_in_min(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_min
             s_in_min(6,ijp1,  k,l) = (        sw(ADM_AJ, n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * AJ_min

             s_in_max(1,ij,    k,l) = (        sw(ADM_AI, n) ) * AI_max         &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * (-CNST_MAX_REAL)
             s_in_max(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_max        &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL)
             s_in_max(3,ij,    k,l) = (        sw(ADM_AJ, n) ) * AJ_max         &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * (-CNST_MAX_REAL)
             s_in_max(4,ip1j,  k,l) = (        sw(ADM_AI, n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * AI_max
             s_in_max(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_max
             s_in_max(6,ijp1,  k,l) = (        sw(ADM_AJ, n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * AJ_max
          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmin-1,ADM_gmin  )

          do n = nstart, nend
             ij     = n
             ip1j   = n+1
             ip1jp1 = n+1+ADM_gall_1d
             ijm1   = n  -ADM_gall_1d

             AI_min  = min( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
             AI_max  = max( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )

             s_in_min(1,ij,    k,l) = (        sw(ADM_AI, n) ) * AI_min        &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * CNST_MAX_REAL
             s_in_min(4,ip1j,  k,l) = (        sw(ADM_AI, n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * AI_min

             s_in_max(1,ij,    k,l) = (        sw(ADM_AI, n) ) * AI_max         &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * (-CNST_MAX_REAL)
             s_in_max(4,ip1j,  k,l) = (        sw(ADM_AI, n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * AI_max

!             if ( c(1,n,k,l) <= 0.D0 ) then
!                s_in_min(1,ij,  k,l) = min( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
!                s_in_max(1,ij,  k,l) = max( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
!             else
!                s_in_min(4,ip1j,k,l) = min( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
!                s_in_max(4,ip1j,k,l) = max( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
!             endif
          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmin-1,ADM_gmin  )

          do n = nstart, nend
             ij     = n
             ip1j   = n+1
             ip1jp1 = n+1+ADM_gall_1d
             ijp1   = n  +ADM_gall_1d

             AIJ_min = min( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
             AIJ_max = max( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )

             s_in_min(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_min       &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * CNST_MAX_REAL
             s_in_min(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_min
             s_in_max(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_max        &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL)
             s_in_max(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_max

!             if ( c(2,n,k,l) <= 0.D0 ) then
!                s_in_min(2,ij,    k,l) = min( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
!                s_in_max(2,ij,    k,l) = max( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
!             else
!                s_in_min(5,ip1jp1,k,l) = min( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
!                s_in_max(5,ip1jp1,k,l) = max( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
!             endif
          enddo

          nstart = suf(ADM_gmin,  ADM_gmin-1)
          nend   = suf(ADM_gmin-1,ADM_gmin  )

          do n = nstart, nend
             ij     = n
             ip1jp1 = n+1+ADM_gall_1d
             ijp1   = n  +ADM_gall_1d
             im1j   = n-1

             AJ_min  = min( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
             AJ_max  = max( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )

             s_in_min(3,ij,    k,l) = (        sw(ADM_AJ, n) ) * AJ_min        &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * CNST_MAX_REAL
             s_in_min(6,ijp1,  k,l) = (        sw(ADM_AJ, n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * AJ_min
             s_in_max(3,ij,    k,l) = (        sw(ADM_AJ, n) ) * AJ_max         &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * (-CNST_MAX_REAL)
             s_in_max(6,ijp1,  k,l) = (        sw(ADM_AJ, n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * AJ_max

!             if ( c(3,n,k,l) <= 0.D0 ) then
!                s_in_min(3,ij,  k,l) = min( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
!                s_in_max(3,ij,  k,l) = max( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
!             else
!                s_in_min(6,ijp1,k,l) = min( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
!                s_in_max(6,ijp1,k,l) = max( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
!             endif
          enddo

          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
             n = suf(ADM_gmin-1,ADM_gmin-1)
             ij     = n
             ip1jp1 = n+1+ADM_gall_1d
             ip2jp1 = n+2+ADM_gall_1d
             ijp1   = n  +ADM_gall_1d

             AIJ_min = min( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
             AIJ_max = max( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )

             s_in_min(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_min       &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * CNST_MAX_REAL
             s_in_min(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_min
             s_in_max(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_max        &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL)
             s_in_max(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_max

!             if ( c(2,n,k,l) <= 0.D0 ) then
!                s_in_min(2,ij,    k,l) = min( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
!                s_in_max(2,ij,    k,l) = max( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
!             else
!                s_in_min(5,ip1jp1,k,l) = min( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
!                s_in_max(5,ip1jp1,k,l) = max( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
!             endif
          endif

       enddo
    enddo

    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             s_in_min_pl(:,k,l,:) = CNST_MAX_REAL
             s_in_max_pl(:,k,l,:) =-CNST_MAX_REAL
             do n=ADM_gmin_pl,ADM_gmax_pl
                if(c_pl(n,k,l)<=0.0D0) then
                   s_in_min_pl(n,k,l,1) = min(s_pl(ADM_gslf_pl,k,l),s_pl(n,k,l),s_pl(nm1(n),k,l),s_pl(np1(n),k,l))
                   s_in_max_pl(n,k,l,1) = max(s_pl(ADM_gslf_pl,k,l),s_pl(n,k,l),s_pl(nm1(n),k,l),s_pl(np1(n),k,l))
                else
                   s_in_min_pl(n,k,l,2) = min(s_pl(ADM_gslf_pl,k,l),s_pl(n,k,l),s_pl(nm1(n),k,l),s_pl(np1(n),k,l))
                   s_in_max_pl(n,k,l,2) = max(s_pl(ADM_gslf_pl,k,l),s_pl(n,k,l),s_pl(nm1(n),k,l),s_pl(np1(n),k,l))
                end if
             enddo

          enddo
          !
       enddo
       !
    end if
    !
    !--- calcluation outflow limiter
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k = 1, ADM_kall
       nstart = suf(ADM_gmin, ADM_gmin )
       nend  = suf(ADM_gmax ,ADM_gmax )
       do n = nstart,nend
          s_m1_k_min_n = min(s_in_min(1,n,k,l),s_in_min(2,n,k,l),s_in_min(3,n,k,l),&
                             s_in_min(4,n,k,l),s_in_min(5,n,k,l),s_in_min(6,n,k,l))
          if(s_m1_k_min_n==CNST_MAX_REAL) s_m1_k_min_n = s(n,k,l)
          s_m1_k_max_n = max(s_in_max(1,n,k,l),s_in_max(2,n,k,l),s_in_max(3,n,k,l),&
                             s_in_max(4,n,k,l),s_in_max(5,n,k,l),s_in_max(6,n,k,l))
          if(s_m1_k_max_n==-CNST_MAX_REAL) s_m1_k_max_n = s(n,k,l)

          c_in_sum_n &
            = (0.5D0-sign(0.5D0,c(1,n,k,l)))*c(1,n,k,l)&
            + (0.5D0-sign(0.5D0,c(2,n,k,l)))*c(2,n,k,l)&
            + (0.5D0-sign(0.5D0,c(3,n,k,l)))*c(3,n,k,l)&
            + (0.5D0-sign(0.5D0,c(4,n,k,l)))*c(4,n,k,l)&
            + (0.5D0-sign(0.5D0,c(5,n,k,l)))*c(5,n,k,l)&
            + (0.5D0-sign(0.5D0,c(6,n,k,l)))*c(6,n,k,l)

          c_out_sum_n &
            = (0.5D0+sign(0.5D0,c(1,n,k,l)))*c(1,n,k,l)&
            + (0.5D0+sign(0.5D0,c(2,n,k,l)))*c(2,n,k,l)&
            + (0.5D0+sign(0.5D0,c(3,n,k,l)))*c(3,n,k,l)&
            + (0.5D0+sign(0.5D0,c(4,n,k,l)))*c(4,n,k,l)&
            + (0.5D0+sign(0.5D0,c(5,n,k,l)))*c(5,n,k,l)&
            + (0.5D0+sign(0.5D0,c(6,n,k,l)))*c(6,n,k,l)
          
          c_qin_sum_max_n &
            = (0.5D0-sign(0.5D0,c(1,n,k,l)))*(c(1,n,k,l)*s_in_max(1,n,k,l))&
            + (0.5D0-sign(0.5D0,c(2,n,k,l)))*(c(2,n,k,l)*s_in_max(2,n,k,l))&
            + (0.5D0-sign(0.5D0,c(3,n,k,l)))*(c(3,n,k,l)*s_in_max(3,n,k,l))&
            + (0.5D0-sign(0.5D0,c(4,n,k,l)))*(c(4,n,k,l)*s_in_max(4,n,k,l))&
            + (0.5D0-sign(0.5D0,c(5,n,k,l)))*(c(5,n,k,l)*s_in_max(5,n,k,l))&
            + (0.5D0-sign(0.5D0,c(6,n,k,l)))*(c(6,n,k,l)*s_in_max(6,n,k,l))

          c_qin_sum_min_n &
            = (0.5D0-sign(0.5D0,c(1,n,k,l)))*(c(1,n,k,l)*s_in_min(1,n,k,l))&
            + (0.5D0-sign(0.5D0,c(2,n,k,l)))*(c(2,n,k,l)*s_in_min(2,n,k,l))&
            + (0.5D0-sign(0.5D0,c(3,n,k,l)))*(c(3,n,k,l)*s_in_min(3,n,k,l))&
            + (0.5D0-sign(0.5D0,c(4,n,k,l)))*(c(4,n,k,l)*s_in_min(4,n,k,l))&
            + (0.5D0-sign(0.5D0,c(5,n,k,l)))*(c(5,n,k,l)*s_in_min(5,n,k,l))&
            + (0.5D0-sign(0.5D0,c(6,n,k,l)))*(c(6,n,k,l)*s_in_min(6,n,k,l))

          if(abs(c_out_sum_n)<CNST_EPS_ZERO) then
            wrk(n,k,l,s_out_k_min) = s(n,k,l)
            wrk(n,k,l,s_out_k_max) = s(n,k,l)
          else
            wrk(n,k,l,s_out_k_min) = ( &
              s(n,k,l)-c_qin_sum_max_n&
              -s_m1_k_max_n*(1.0D0-c_in_sum_n-c_out_sum_n+d(n,k,l)) &
              )/c_out_sum_n
            wrk(n,k,l,s_out_k_max) = ( &
              s(n,k,l)-c_qin_sum_min_n&
              -s_m1_k_min_n*(1.0D0-c_in_sum_n-c_out_sum_n+d(n,k,l)) &
              )/c_out_sum_n
          end if
       enddo !N
    enddo !K
    enddo !L
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             s_m1_k_min_pl &
                  = min(s_in_min_pl(ADM_gmin_pl,k,l,1),s_in_min_pl(ADM_gmin_pl+1,k,l,1),&
                  s_in_min_pl(ADM_gmin_pl+2,k,l,1),s_in_min_pl(ADM_gmin_pl+3,k,l,1),s_in_min_pl(ADM_gmin_pl+4,k,l,1))
             if(s_m1_k_min_pl== CNST_MAX_REAL) s_m1_k_min_pl=s_pl(ADM_gslf_pl,k,l)
             !          s_m1_k_min_pl = max(SMALL_ZERO,s_m1_k_min_pl)
             s_m1_k_max_pl &
                  = max(s_in_min_pl(ADM_gmin_pl,k,l,1),s_in_min_pl(ADM_gmin_pl+1,k,l,1),&
                  s_in_min_pl(ADM_gmin_pl+2,k,l,1),s_in_min_pl(ADM_gmin_pl+3,k,l,1),s_in_min_pl(ADM_gmin_pl+4,k,l,1))
             if(s_m1_k_max_pl==-CNST_MAX_REAL) s_m1_k_max_pl=s_pl(ADM_gslf_pl,k,l)
             !
             c_in_sum_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl  ,k,l)))*c_pl(ADM_gmin_pl  ,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+1,k,l)))*c_pl(ADM_gmin_pl+1,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+2,k,l)))*c_pl(ADM_gmin_pl+2,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+3,k,l)))*c_pl(ADM_gmin_pl+3,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+4,k,l)))*c_pl(ADM_gmin_pl+4,k,l)
             c_out_sum_pl &
                  = (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl  ,k,l)))*c_pl(ADM_gmin_pl  ,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl+1,k,l)))*c_pl(ADM_gmin_pl+1,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl+2,k,l)))*c_pl(ADM_gmin_pl+2,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl+3,k,l)))*c_pl(ADM_gmin_pl+3,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl+4,k,l)))*c_pl(ADM_gmin_pl+4,k,l)
             c_qin_sum_max_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl  ,k,l)))*(c_pl(ADM_gmin_pl  ,k,l)*s_in_max_pl(ADM_gmin_pl  ,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+1,k,l)))*(c_pl(ADM_gmin_pl+1,k,l)*s_in_max_pl(ADM_gmin_pl+1,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+2,k,l)))*(c_pl(ADM_gmin_pl+2,k,l)*s_in_max_pl(ADM_gmin_pl+2,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+3,k,l)))*(c_pl(ADM_gmin_pl+3,k,l)*s_in_max_pl(ADM_gmin_pl+3,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+4,k,l)))*(c_pl(ADM_gmin_pl+4,k,l)*s_in_max_pl(ADM_gmin_pl+4,k,l,1))
             c_qin_sum_min_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl  ,k,l)))*(c_pl(ADM_gmin_pl  ,k,l)*s_in_min_pl(ADM_gmin_pl  ,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+1,k,l)))*(c_pl(ADM_gmin_pl+1,k,l)*s_in_min_pl(ADM_gmin_pl+1,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+2,k,l)))*(c_pl(ADM_gmin_pl+2,k,l)*s_in_min_pl(ADM_gmin_pl+2,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+3,k,l)))*(c_pl(ADM_gmin_pl+3,k,l)*s_in_min_pl(ADM_gmin_pl+3,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+4,k,l)))*(c_pl(ADM_gmin_pl+4,k,l)*s_in_min_pl(ADM_gmin_pl+4,k,l,1))
             !
             if(abs(c_out_sum_pl)<CNST_EPS_ZERO) then
                wrk_pl(ADM_gslf_pl,k,l,s_out_k_min) = s_pl(ADM_gslf_pl,k,l)
                wrk_pl(ADM_gslf_pl,k,l,s_out_k_max) = s_pl(ADM_gslf_pl,k,l)
             else
                wrk_pl(ADM_gslf_pl,k,l,s_out_k_min) = ( &
                     s_pl(ADM_gslf_pl,k,l)-c_qin_sum_max_pl&
                     -s_m1_k_max_pl*(1.0D0-c_in_sum_pl-c_out_sum_pl+d_pl(ADM_gslf_pl,k,l)) &
                     )/c_out_sum_pl
                wrk_pl(ADM_gslf_pl,k,l,s_out_k_max) = ( &
                     s_pl(ADM_gslf_pl,k,l)-c_qin_sum_min_pl&
                     -s_m1_k_min_pl*(1.0D0-c_in_sum_pl-c_out_sum_pl+d_pl(ADM_gslf_pl,k,l)) &
                     )/c_out_sum_pl
             end if
          enddo
       enddo
    endif
    !
    !
1000 continue
    !
    !--- H.Tomita 090414 
    do l=1,ADM_lall
      do k=1,ADM_kall
        do n=1,ADM_gall
          wrk(n,k,l,dsx:dsz)=0.0D0
        enddo
      enddo
    enddo

    do l=1,ADM_lall_pl
      do k=1,ADM_kall
        do n=1,ADM_gall_pl
          wrk_pl(n,k,l,dsx:dsz)=0.0D0
        enddo
      enddo
    enddo

    call OPRT_gradient(                    &
         wrk(:,:,:,dsx), wrk_pl(:,:,:,dsx),&
         wrk(:,:,:,dsy), wrk_pl(:,:,:,dsy),&
         wrk(:,:,:,dsz), wrk_pl(:,:,:,dsz),&
         s, s_pl )

    call COMM_data_transfer(wrk,wrk_pl)

    !--- basic scheme
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       do k = 1, ADM_kall
        nstart = suf(ADM_gmin-1,ADM_gmin-1 )
        nend  = suf(ADM_gmax, ADM_gmax )
        do n = nstart, nend
          sa_p = s(n,k,l)  &
            +wrk(n,k,l,dsx)*(cp(n,k,l,ADM_AI,GRD_XDIR)-GRD_x(n,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n,k,l,dsy)*(cp(n,k,l,ADM_AI,GRD_YDIR)-GRD_x(n,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n,k,l,dsz)*(cp(n,k,l,ADM_AI,GRD_ZDIR)-GRD_x(n,ADM_KNONE,l,GRD_ZDIR))
          sa_m = s(n+1,k,l) &
            +wrk(n+1,k,l,dsx)*(cp(n,k,l,ADM_AI,GRD_XDIR)-GRD_x(n+1,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n+1,k,l,dsy)*(cp(n,k,l,ADM_AI,GRD_YDIR)-GRD_x(n+1,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n+1,k,l,dsz)*(cp(n,k,l,ADM_AI,GRD_ZDIR)-GRD_x(n+1,ADM_KNONE,l,GRD_ZDIR))
          sa(ADM_AI,n,k,l) &
            =(0.5D0+sign(0.5D0,c(1,n,k,l)))*sa_p+(0.5D0-sign(0.5D0,c(1,n,k,l)))*sa_m

          sa_p = s(n,k,l) &
            +wrk(n,k,l,dsx)*(cp(n,k,l,ADM_AIJ,GRD_XDIR)-GRD_x(n,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n,k,l,dsy)*(cp(n,k,l,ADM_AIJ,GRD_YDIR)-GRD_x(n,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n,k,l,dsz)*(cp(n,k,l,ADM_AIJ,GRD_ZDIR)-GRD_x(n,ADM_KNONE,l,GRD_ZDIR))
          sa_m = s(n+1+ADM_gall_1d,k,l) &
            +wrk(n+1+ADM_gall_1d,k,l,dsx)*(cp(n,k,l,ADM_AIJ,GRD_XDIR)- &
            GRD_x(n+1+ADM_gall_1d,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n+1+ADM_gall_1d,k,l,dsy)*(cp(n,k,l,ADM_AIJ,GRD_YDIR)- &
            GRD_x(n+1+ADM_gall_1d,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n+1+ADM_gall_1d,k,l,dsz)*(cp(n,k,l,ADM_AIJ,GRD_ZDIR)- &
            GRD_x(n+1+ADM_gall_1d,ADM_KNONE,l,GRD_ZDIR))
          sa(ADM_AIJ,n,k,l) &
            =(0.5D0+sign(0.5D0,c(2,n,k,l)))*sa_p+(0.5D0-sign(0.5D0,c(2,n,k,l)))*sa_m

          sa_p = s(n,k,l) &
            +wrk(n,k,l,dsx)*(cp(n,k,l,ADM_AJ,GRD_XDIR)-GRD_x(n,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n,k,l,dsy)*(cp(n,k,l,ADM_AJ,GRD_YDIR)-GRD_x(n,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n,k,l,dsz)*(cp(n,k,l,ADM_AJ,GRD_ZDIR)-GRD_x(n,ADM_KNONE,l,GRD_ZDIR))
          sa_m = s(n+ADM_gall_1d,k,l) &
            +wrk(n+ADM_gall_1d,k,l,dsx)*(cp(n,k,l,ADM_AJ,GRD_XDIR)- &
            GRD_x(n+ADM_gall_1d,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n+ADM_gall_1d,k,l,dsy)*(cp(n,k,l,ADM_AJ,GRD_YDIR)- GRD_x(n+ADM_gall_1d,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n+ADM_gall_1d,k,l,dsz)*(cp(n,k,l,ADM_AJ,GRD_ZDIR)- GRD_x(n+ADM_gall_1d,ADM_KNONE,l,GRD_ZDIR))
          sa(ADM_AJ,n,k,l) &
            =(0.5D0+sign(0.5D0,c(3,n,k,l)))*sa_p+(0.5D0-sign(0.5D0,c(3,n,k,l)))*sa_m
        enddo !N
      enddo !K
    enddo !L
    !
    if(ADM_prc_me==ADM_prc_pl) then
       !
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n=ADM_gmin_pl,ADM_gmax_pl
                sa_p = s_pl(ADM_gslf_pl,k,l) &
                     +wrk_pl(ADM_gslf_pl,k,l,dsx)*(cp_pl(n,k,l,GRD_XDIR)-GRD_x_pl(ADM_gslf_pl,ADM_KNONE,l,GRD_XDIR)) &
                     +wrk_pl(ADM_gslf_pl,k,l,dsy)*(cp_pl(n,k,l,GRD_YDIR)-GRD_x_pl(ADM_gslf_pl,ADM_KNONE,l,GRD_YDIR)) &
                     +wrk_pl(ADM_gslf_pl,k,l,dsz)*(cp_pl(n,k,l,GRD_ZDIR)-GRD_x_pl(ADM_gslf_pl,ADM_KNONE,l,GRD_ZDIR))
                sa_m = s_pl(n,k,l) &
                     +wrk_pl(n,k,l,dsx)*(cp_pl(n,k,l,GRD_XDIR)-GRD_x_pl(n,ADM_KNONE,l,GRD_XDIR)) &
                     +wrk_pl(n,k,l,dsy)*(cp_pl(n,k,l,GRD_YDIR)-GRD_x_pl(n,ADM_KNONE,l,GRD_YDIR)) &
                     +wrk_pl(n,k,l,dsz)*(cp_pl(n,k,l,GRD_ZDIR)-GRD_x_pl(n,ADM_KNONE,l,GRD_ZDIR))
                sa_pl(n,k,l) &
                     =(0.5D0+sign(0.5D0,c_pl(n,k,l)))*sa_p+(0.5D0-sign(0.5D0,c_pl(n,k,l)))*sa_m
             enddo
          enddo
       enddo
       !
    end if
    !
    if(present(limiter)) then
       if(trim(limiter)=='NON_LIM') goto 2000
    end if
    !
    !---- apply inflow limiter
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k = 1, ADM_kall
        nstart = suf(ADM_gmin-1,ADM_gmin-1 )
        nend = suf(ADM_gmax ,ADM_gmax )
        do n=nstart,nend
          sa(ADM_AI,n,k,l) &
            =(0.5D0-sign(0.5D0,c(1,n,k,l)))&
            *min(max(sa(ADM_AI,n,k,l),s_in_min(1,n,k,l)),s_in_max(1,n,k,l))&
            +(0.5D0+sign(0.5D0,c(1,n,k,l)))&
            *min(max(sa(ADM_AI,n,k,l),s_in_min(4,n+1,k,l)),s_in_max(4,n+1,k,l))

          sa(ADM_AIJ,n,k,l) &
            =(0.5D0-sign(0.5D0,c(2,n,k,l)))&
            *min(max(sa(ADM_AIJ,n,k,l),s_in_min(2,n,k,l)),s_in_max(2,n,k,l))&
            +(0.5D0+sign(0.5D0,c(2,n,k,l)))&
            *min(max(sa(ADM_AIJ,n,k,l),s_in_min(5,n+1+ADM_gall_1d,k,l)),s_in_max(5,n+1+ADM_gall_1d,k,l))

          sa(ADM_AJ,n,k,l) &
            =(0.5D0-sign(0.5D0,c(3,n,k,l)))&
            *min(max(sa(ADM_AJ,n,k,l),s_in_min(3,n,k,l)),s_in_max(3,n,k,l))&
            +(0.5D0+sign(0.5D0,c(3,n,k,l)))&
            *min(max(sa(ADM_AJ,n,k,l),s_in_min(6,n+ADM_gall_1d,k,l)),s_in_max(6,n+ADM_gall_1d,k,l))
        enddo
      enddo
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n=ADM_gmin_pl,ADM_gmax_pl
                sa_pl(n,k,l) &
                  =(0.5D0-sign(0.5D0,c_pl(n,k,l)))&
                  *min(max(sa_pl(n,k,l),s_in_min_pl(n,k,l,1)),s_in_max_pl(n,k,l,1))&
                  +(0.5D0+sign(0.5D0,c_pl(n,k,l)))&
                  *min(max(sa_pl(n,k,l),s_in_min_pl(n,k,l,2)),s_in_max_pl(n,k,l,2))
             enddo
          enddo
       enddo
    end if
    !
    !---- apply outflow limitter
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k = 1, ADM_kall
        nstart = suf(ADM_gmin-1,ADM_gmin-1 )
        nend = suf(ADM_gmax ,ADM_gmax )
        do n = nstart,nend
          sa(ADM_AI,n,k,l) &
            =(0.5D0-sign(0.5D0,c(1,n,k,l)))&
            *max(min(sa(ADM_AI,n,k,l),wrk(n+1,k,l,s_out_k_max)),wrk(n+1,k,l,s_out_k_min))&
            +(0.5D0+sign(0.5D0,c(1,n,k,l)))&
            *max(min(sa(ADM_AI,n,k,l),wrk(n,k,l,s_out_k_max)),wrk(n,k,l,s_out_k_min))

          sa(ADM_AIJ,n,k,l) &
            =(0.5D0-sign(0.5D0,c(2,n,k,l)))&
            *max(min(sa(ADM_AIJ,n,k,l),wrk(n+1+ADM_gall_1d,k,l,s_out_k_max)),wrk(n+1+ADM_gall_1d,k,l,s_out_k_min))&
            +(0.5D0+sign(0.5D0,c(2,n,k,l)))&
            *max(min(sa(ADM_AIJ,n,k,l),wrk(n,k,l,s_out_k_max)),wrk(n,k,l,s_out_k_min))

          sa(ADM_AJ,n,k,l) &
            =(0.5D0-sign(0.5D0,c(3,n,k,l)))&
            *max(min(sa(ADM_AJ,n,k,l),wrk(n+ADM_gall_1d,k,l,s_out_k_max)),wrk(n+ADM_gall_1d,k,l,s_out_k_min))&
            +(0.5D0+sign(0.5D0,c(3,n,k,l)))&
            *max(min(sa(ADM_AJ,n,k,l),wrk(n,k,l,s_out_k_max)),wrk(n,k,l,s_out_k_min))
        enddo
       enddo
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n = ADM_gmin_pl,ADM_gmax_pl
                sa_pl(n,k,l)&
                     =(0.5D0-sign(0.5D0,c_pl(n,k,l)))&
                     *max(min(sa_pl(n,k,l),wrk_pl(n,k,l,s_out_k_max)),wrk_pl(n,k,l,s_out_k_min))&
                     +(0.5D0+sign(0.5D0,c_pl(n,k,l)))&
                     *max(min(sa_pl(n,k,l),wrk_pl(ADM_gslf_pl,k,l,s_out_k_max)),wrk_pl(ADM_gslf_pl,k,l,s_out_k_min))
             enddo
          enddo
       enddo
    end if
    !
    !
2000 continue
    !
    !--- update
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       nstart = suf(ADM_gmin  ,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )
       do k = 1, ADM_kall
          do n = nstart,nend
             scl(n,k,l) = &
                  ( flx_h(1,n,k,l)*sa(ADM_AI,n,k,l)   &
                  + flx_h(2,n,k,l)*sa(ADM_AIJ,n,k,l)  &
                  + flx_h(3,n,k,l)*sa(ADM_AJ,n,k,l)   &
                  + flx_h(4,n,k,l)*sa(ADM_AI,n-1,k,l) &
                  + flx_h(5,n,k,l)*sa(ADM_AIJ,n-1-ADM_gall_1d,k,l) &
                  + flx_h(6,n,k,l)*sa(ADM_AJ,n-ADM_gall_1d,k,l)    &
                  ) * fact
          enddo
       enddo
    enddo

    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             scl_pl(ADM_gslf_pl,k,l)=  &
                  ( flx_h_pl(ADM_gmin_pl  ,k,l)*sa_pl(ADM_gmin_pl  ,k,l) &
                  + flx_h_pl(ADM_gmin_pl+1,k,l)*sa_pl(ADM_gmin_pl+1,k,l) &
                  + flx_h_pl(ADM_gmin_pl+2,k,l)*sa_pl(ADM_gmin_pl+2,k,l) &
                  + flx_h_pl(ADM_gmin_pl+3,k,l)*sa_pl(ADM_gmin_pl+3,k,l) &
                  + flx_h_pl(ADM_gmin_pl+4,k,l)*sa_pl(ADM_gmin_pl+4,k,l) &
                  ) * fact
          enddo
       enddo
    end if

    return
  end subroutine OPRT_divergence2

end module mod_trcadv_thuburn
!-------------------------------------------------------------------------------------
