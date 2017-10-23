!-------------------------------------------------------------------------------
!> Module vertical implicit scheme
!!
!! @par Description
!!          This module is for the caluculation of vertical implicit scheme
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_vi
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof

  use mod_adm, only: &
     ADM_lall,    &
     ADM_lall_pl, &
     ADM_gall,    &
     ADM_gall_pl, &
     ADM_kall
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: vi_setup
  public :: vi_small_step

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: vi_main
  private :: vi_rhow_update_matrix
  private :: vi_rhow_solver

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
#ifdef _FIXEDINDEX_
  real(RP), public               :: Mc   (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), private              :: Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public               :: Ml   (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), private              :: Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public               :: Mu   (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), private              :: Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
#else
  real(RP), public,  allocatable :: Mc   (:,:,:)
  real(RP), private, allocatable :: Mc_pl(:,:,:)
  real(RP), public,  allocatable :: Ml   (:,:,:)
  real(RP), private, allocatable :: Ml_pl(:,:,:)
  real(RP), public,  allocatable :: Mu   (:,:,:)
  real(RP), private, allocatable :: Mu_pl(:,:,:)
#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine vi_setup
    implicit none
    !---------------------------------------------------------------------------

#ifndef _FIXEDINDEX_
    allocate( Mc   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Mu   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Ml   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
#endif

    return
  end subroutine vi_setup

  !-----------------------------------------------------------------------------
  !> Driver of the small step integration
  subroutine vi_small_step( &
       PROG,       PROG_pl,       &
       vx,         vx_pl,         &
       vy,         vy_pl,         &
       vz,         vz_pl,         &
       eth,        eth_pl,        &
       rhog_prim,  rhog_prim_pl,  &
       preg_prim,  preg_prim_pl,  &
       g_TEND0,    g_TEND0_pl,    &
       PROG_split, PROG_split_pl, &
       PROG_mean,  PROG_mean_pl,  &
       num_of_itr,                &
       dt                         )
    use mod_const, only: &
       CONST_GRAV, &
       CONST_Rdry, &
       CONST_CVdry
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_kmax,    &
       ADM_kmin
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR, &
       GRD_afact,        &
       GRD_bfact
    use mod_oprt, only: &
       OPRT_horizontalize_vec
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl, &
       VMTR_W2Cfact,    &
       VMTR_W2Cfact_pl
    use mod_time, only: &
       TIME_SPLIT
    use mod_runconf, only: &
       I_RHOG,          &
       I_RHOGVX,        &
       I_RHOGVY,        &
       I_RHOGVZ,        &
       I_RHOGW,         &
       I_RHOGE,         &
       NON_HYDRO_ALPHA
    use mod_bndcnd, only: &
       BNDCND_rhovxvyvz
    use mod_numfilter, only: &
       numfilter_divdamp,   &
       numfilter_divdamp_2d
    use mod_src, only: &
       src_advection_convergence, &
       src_flux_convergence,      &
       src_pres_gradient,         &
       src_buoyancy,              &
       I_SRC_default,             &
       I_SRC_horizontal
    implicit none

    real(RP), intent(inout) :: PROG         (ADM_gall   ,ADM_kall,ADM_lall   ,6) ! prognostic variables
    real(RP), intent(inout) :: PROG_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    real(RP), intent(in)    :: vx           (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_x
    real(RP), intent(in)    :: vx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vy           (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_y
    real(RP), intent(in)    :: vy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vz           (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_z
    real(RP), intent(in)    :: vz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: eth          (ADM_gall   ,ADM_kall,ADM_lall   ) ! enthalpy
    real(RP), intent(in)    :: eth_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog_prim    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog_prim_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg_prim    (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: preg_prim_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)    :: g_TEND0      (ADM_gall   ,ADM_kall,ADM_lall   ,6) ! tendency of prognostic variables
    real(RP), intent(in)    :: g_TEND0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    real(RP), intent(inout) :: PROG_split   (ADM_gall   ,ADM_kall,ADM_lall   ,6) ! prognostic variables (split)
    real(RP), intent(inout) :: PROG_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    real(RP), intent(out)   :: PROG_mean    (ADM_gall   ,ADM_kall,ADM_lall   ,5) ! mean_flux for tracer advection
    real(RP), intent(out)   :: PROG_mean_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,5)

    integer,  intent(in)    :: num_of_itr
    real(RP), intent(in)    :: dt

    ! merged array for communication
    real(RP) :: diff_vh      (ADM_gall   ,ADM_kall,ADM_lall   ,3)
    real(RP) :: diff_vh_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,3)
    real(RP) :: diff_we      (ADM_gall   ,ADM_kall,ADM_lall   ,3)
    real(RP) :: diff_we_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,3)

    ! tendency term (large step)
    real(RP) :: grhogetot0   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhogetot0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    ! tendency term (large step + small step)
    real(RP) :: g_TEND       (ADM_gall   ,ADM_kall,ADM_lall   ,6) ! tendency of prognostic variables
    real(RP) :: g_TEND_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    ! tendency term 2
    real(RP) :: drhog        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhog_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhogvx
    real(RP) :: drhogvy
    real(RP) :: drhogvz
    real(RP) :: drhogw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhogw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhoge       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhoge_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    ! divergence damping
    real(RP) :: ddivdvx      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdvy      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdvz      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: ddivdvx_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvx_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdvy_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvy_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdvz_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvz_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    ! pressure gradient force
    real(RP) :: dpgrad       (ADM_gall   ,ADM_kall,ADM_lall   ,XDIR:ZDIR)
    real(RP) :: dpgrad_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,XDIR:ZDIR)
    real(RP) :: dpgradw      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dpgradw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    ! buoyancy force
    real(RP) :: dbuoiw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dbuoiw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    ! pressure work
    real(RP) :: drhoge_pw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhoge_pw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhoge_pwh   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhoge_pwh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: gz_tilde     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: gz_tilde_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhog_h       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_h_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: eth_h        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: eth_h_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: preg_prim_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: preg_prim_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rweight_itr

    integer  :: gall, kall, kmin, kmax, lall
    real(RP) :: grav, RovCV, alpha

    integer  :: g, k, l, iv, ns
    !---------------------------------------------------------------------------

    call PROF_rapstart('____vi_path0',2)

    gall = ADM_gall
    kall = ADM_kall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    grav  = CONST_GRAV
    RovCV = CONST_Rdry / CONST_CVdry
    alpha = real(NON_HYDRO_ALPHA,kind=RP)

!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,grhogetot0,g_TEND0), &
    !$omp collapse(2)
    do l = 1, lall
    do k = 1, kall
    do g = 1, gall
       grhogetot0(g,k,l) = g_TEND0(g,k,l,I_RHOGE)
    enddo
    enddo
    enddo
    !$omp end parallel do

    grhogetot0_pl(:,:,:) = g_TEND0_pl(:,:,:,I_RHOGE)

    ! full level -> half level
    !$omp parallel default(none),private(g,k,l), &
    !$omp shared(gall,kall,kmin,kmax,lall,rhog_h,eth_h,   &
    !$omp        PROG,eth,GRD_afact,GRD_bfact,VMTR_C2Wfact)
    do l = 1, lall
!OCL XFILL
       !$omp do
       do k = kmin, kmax+1
       do g = 1, gall
          rhog_h(g,k,l) = ( VMTR_C2Wfact(g,k,1,l) * PROG(g,k,  l,I_RHOG) &
                          + VMTR_C2Wfact(g,k,2,l) * PROG(g,k-1,l,I_RHOG) )
          eth_h (g,k,l) = GRD_afact(k) * eth(g,k,  l) &
                        + GRD_bfact(k) * eth(g,k-1,l)
       enddo
       enddo
       !$omp end do

       !$omp do
       do g = 1, gall
          rhog_h(g,kmin-1,l) = rhog_h(g,kmin,l)
          eth_h (g,kmin-1,l) = eth_h (g,kmin,l)
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,k,1,l) * PROG_pl(g,k,  l,I_RHOG) &
                                + VMTR_C2Wfact_pl(g,k,2,l) * PROG_pl(g,k-1,l,I_RHOG) )
             eth_h_pl (g,k,l) = GRD_afact(k) * eth_pl(g,k,  l) &
                              + GRD_bfact(k) * eth_pl(g,k-1,l)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,ADM_kmin-1,l) = rhog_h_pl(g,ADM_kmin,l)
             eth_h_pl (g,ADM_kmin-1,l) = eth_h_pl (g,ADM_kmin,l)
          enddo
       enddo
    endif

    !---< Calculation of source term for rhog >

    call src_flux_convergence( PROG (:,:,:,I_RHOGVX), PROG_pl (:,:,:,I_RHOGVX), & ! [IN]
                               PROG (:,:,:,I_RHOGVY), PROG_pl (:,:,:,I_RHOGVY), & ! [IN]
                               PROG (:,:,:,I_RHOGVZ), PROG_pl (:,:,:,I_RHOGVZ), & ! [IN]
                               PROG (:,:,:,I_RHOGW),  PROG_pl (:,:,:,I_RHOGW),  & ! [IN]
                               drhog(:,:,:),          drhog_pl(:,:,:),          & ! [OUT]
                               I_SRC_default                                    ) ! [IN]

    !---< Calculation of source term for Vh(vx,vy,vz) and W >

    ! divergence damping
    call numfilter_divdamp( PROG   (:,:,:,I_RHOGVX), PROG_pl   (:,:,:,I_RHOGVX), & ! [IN]
                            PROG   (:,:,:,I_RHOGVY), PROG_pl   (:,:,:,I_RHOGVY), & ! [IN]
                            PROG   (:,:,:,I_RHOGVZ), PROG_pl   (:,:,:,I_RHOGVZ), & ! [IN]
                            PROG   (:,:,:,I_RHOGW),  PROG_pl   (:,:,:,I_RHOGW),  & ! [IN]
                            ddivdvx(:,:,:),          ddivdvx_pl(:,:,:),          & ! [OUT]
                            ddivdvy(:,:,:),          ddivdvy_pl(:,:,:),          & ! [OUT]
                            ddivdvz(:,:,:),          ddivdvz_pl(:,:,:),          & ! [OUT]
                            ddivdw (:,:,:),          ddivdw_pl (:,:,:)           ) ! [OUT]

    call numfilter_divdamp_2d( PROG      (:,:,:,I_RHOGVX), PROG_pl      (:,:,:,I_RHOGVX), & ! [IN]
                               PROG      (:,:,:,I_RHOGVY), PROG_pl      (:,:,:,I_RHOGVY), & ! [IN]
                               PROG      (:,:,:,I_RHOGVZ), PROG_pl      (:,:,:,I_RHOGVZ), & ! [IN]
                               ddivdvx_2d(:,:,:),          ddivdvx_2d_pl(:,:,:),          & ! [OUT]
                               ddivdvy_2d(:,:,:),          ddivdvy_2d_pl(:,:,:),          & ! [OUT]
                               ddivdvz_2d(:,:,:),          ddivdvz_2d_pl(:,:,:)           ) ! [OUT]

    ! pressure force
    call src_pres_gradient( preg_prim(:,:,:),   preg_prim_pl(:,:,:),   & ! [IN]
                            dpgrad   (:,:,:,:), dpgrad_pl   (:,:,:,:), & ! [OUT]
                            dpgradw  (:,:,:),   dpgradw_pl  (:,:,:),   & ! [OUT]
                            I_SRC_default                              ) ! [IN]

    ! buoyancy force
    call src_buoyancy( rhog_prim(:,:,:), rhog_prim_pl(:,:,:), & ! [IN]
                       dbuoiw   (:,:,:), dbuoiw_pl   (:,:,:)  ) ! [OUT]

    !---< Calculation of source term for rhoge >

    ! advection convergence for eth
    call src_advection_convergence( PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & ! [IN]
                                    PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & ! [IN]
                                    PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & ! [IN]
                                    PROG  (:,:,:,I_RHOGW),  PROG_pl  (:,:,:,I_RHOGW),  & ! [IN]
                                    eth   (:,:,:),          eth_pl   (:,:,:),          & ! [IN]
                                    drhoge(:,:,:),          drhoge_pl(:,:,:),          & ! [OUT]
                                    I_SRC_default                                      ) ! [IN]

    !$omp parallel default(none),private(g,k,l), &
    !$omp shared(gall,kall,kmin,kmax,lall,drhoge_pw,drhoge_pwh,gz_tilde,        &
    !$omp        PROG,vx,vy,vz,dpgrad,dpgradw,dbuoiw,rhog_h,VMTR_W2Cfact,grav)
    do l = 1, lall
       !$omp do
       do k = 1, kall
       do g = 1, gall
          gz_tilde  (g,k,l) = grav - ( dpgradw(g,k,l)-dbuoiw(g,k,l) ) / rhog_h(g,k,l)
          drhoge_pwh(g,k,l) = -gz_tilde(g,k,l) * PROG(g,k,l,I_RHOGW)
       enddo
       enddo
       !$omp end do

!OCL XFILL
       !$omp do
       do k = kmin, kmax
       do g = 1, gall
          drhoge_pw(g,k,l) = ( vx(g,k,l) * dpgrad(g,k,l,XDIR)              &
                             + vy(g,k,l) * dpgrad(g,k,l,YDIR)              &
                             + vz(g,k,l) * dpgrad(g,k,l,ZDIR)              ) &
                           + ( VMTR_W2Cfact(g,k,1,l) * drhoge_pwh(g,k+1,l) &
                             + VMTR_W2Cfact(g,k,2,l) * drhoge_pwh(g,k,  l) )
       enddo
       enddo
       !$omp end do nowait

!OCL XFILL
       !$omp do
       do g = 1, gall
          drhoge_pw(g,kmin-1,l) = 0.0_RP
          drhoge_pw(g,kmax+1,l) = 0.0_RP
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             gz_tilde_pl  (g,k,l) = GRAV - ( dpgradw_pl(g,k,l)-dbuoiw_pl(g,k,l) ) / rhog_h_pl(g,k,l)
             drhoge_pwh_pl(g,k,l) = -gz_tilde_pl(g,k,l) * PROG_pl(g,k,l,I_RHOGW)
          enddo
          enddo

          do k  = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             drhoge_pw_pl(g,k,l) = ( vx_pl(g,k,l) * dpgrad_pl(g,k,l,XDIR)              &
                                   + vy_pl(g,k,l) * dpgrad_pl(g,k,l,YDIR)              &
                                   + vz_pl(g,k,l) * dpgrad_pl(g,k,l,ZDIR)              ) &
                                 + ( VMTR_W2Cfact_pl(g,k,1,l) * drhoge_pwh_pl(g,k+1,l) &
                                   + VMTR_W2Cfact_pl(g,k,2,l) * drhoge_pwh_pl(g,k,  l) )
          enddo
          enddo
          do g = 1, ADM_gall_pl
             drhoge_pw_pl(g,ADM_kmin-1,l) = 0.0_RP
             drhoge_pw_pl(g,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    endif

    !---< sum of tendencies ( large step + pres-grad + div-damp + div-damp_2d + buoyancy ) >
!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,g_TEND,g_TEND0,drhog,dpgradw,dbuoiw,drhoge,drhoge_pw,        &
    !$omp        dpgrad,ddivdvx,ddivdvy,ddivdvz,ddivdw,ddivdvx_2d,ddivdvy_2d,ddivdvz_2d,alpha), &
    !$omp collapse(2)
    do l = 1, lall
    do k = 1, kall
    do g = 1, gall
       g_TEND(g,k,l,I_RHOG)   = g_TEND0(g,k,l,I_RHOG)   + drhog  (g,k,l)
       g_TEND(g,k,l,I_RHOGVX) = g_TEND0(g,k,l,I_RHOGVX) - dpgrad (g,k,l,XDIR) + ddivdvx(g,k,l) + ddivdvx_2d(g,k,l)
       g_TEND(g,k,l,I_RHOGVY) = g_TEND0(g,k,l,I_RHOGVY) - dpgrad (g,k,l,YDIR) + ddivdvy(g,k,l) + ddivdvy_2d(g,k,l)
       g_TEND(g,k,l,I_RHOGVZ) = g_TEND0(g,k,l,I_RHOGVZ) - dpgrad (g,k,l,ZDIR) + ddivdvz(g,k,l) + ddivdvz_2d(g,k,l)
       g_TEND(g,k,l,I_RHOGW)  = g_TEND0(g,k,l,I_RHOGW)                        + ddivdw (g,k,l) * alpha &
                                                        - dpgradw(g,k,l)      + dbuoiw(g,k,l)
       g_TEND(g,k,l,I_RHOGE)  = g_TEND0(g,k,l,I_RHOGE)  + drhoge (g,k,l)      + drhoge_pw(g,k,l)
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          g_TEND_pl(g,k,l,I_RHOG)   = g_TEND0_pl(g,k,l,I_RHOG)   + drhog_pl  (g,k,l)
          g_TEND_pl(g,k,l,I_RHOGVX) = g_TEND0_pl(g,k,l,I_RHOGVX) - dpgrad_pl (g,k,l,XDIR) + ddivdvx_pl(g,k,l) + ddivdvx_2d_pl(g,k,l)
          g_TEND_pl(g,k,l,I_RHOGVY) = g_TEND0_pl(g,k,l,I_RHOGVY) - dpgrad_pl (g,k,l,YDIR) + ddivdvy_pl(g,k,l) + ddivdvy_2d_pl(g,k,l)
          g_TEND_pl(g,k,l,I_RHOGVZ) = g_TEND0_pl(g,k,l,I_RHOGVZ) - dpgrad_pl (g,k,l,ZDIR) + ddivdvz_pl(g,k,l) + ddivdvz_2d_pl(g,k,l)
          g_TEND_pl(g,k,l,I_RHOGW)  = g_TEND0_pl(g,k,l,I_RHOGW)                           + ddivdw_pl (g,k,l) * alpha &
                                                                 - dpgradw_pl(g,k,l)      + dbuoiw_pl(g,k,l)
          g_TEND_pl(g,k,l,I_RHOGE)  = g_TEND0_pl(g,k,l,I_RHOGE)  + drhoge_pl (g,k,l)      + drhoge_pw_pl(g,k,l)
       enddo
       enddo
       enddo
    endif

    ! initialization of mean mass flux
    rweight_itr = 1.0_RP / real(num_of_itr,kind=RP)

!OCL XFILL
    !$omp parallel do default(none),private(g,k,l,iv), &
    !$omp shared(gall,kall,lall,PROG_mean,PROG), &
    !$omp collapse(3)
    do iv = I_RHOG, I_RHOGW
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       PROG_mean(g,k,l,iv) = PROG(g,k,l,iv)
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel do

    PROG_mean_pl(:,:,:,I_RHOG:I_RHOGW)   = PROG_pl(:,:,:,I_RHOG:I_RHOGW)

    ! update working matrix for vertical implicit solver
    call vi_rhow_update_matrix( eth_h   (:,:,:), eth_h_pl   (:,:,:), & ! [IN]
                                gz_tilde(:,:,:), gz_tilde_pl(:,:,:), & ! [IN]
                                dt                                   ) ! [IN]

    call PROF_rapend  ('____vi_path0',2)
    !---------------------------------------------------------------------------
    !
    !> Start small step iteration
    !
    !---------------------------------------------------------------------------
    do ns = 1, num_of_itr

       call PROF_rapstart('____vi_path1',2)

       !---< calculation of preg_prim(*) from rhog(*) & rhoge(*) >

       !$omp parallel default(none),private(g,k,l), &
       !$omp shared(gall,kall,kmin,kmax,lall,preg_prim_split,PROG_split,RovCV)
       do l = 1, lall
!OCL XFILL
          !$omp do
          do k = 1, kall
          do g = 1, gall
             preg_prim_split(g,k,l) = PROG_split(g,k,l,I_RHOGE) * RovCV
          enddo
          enddo
          !$omp end do

          !$omp do
          do g = 1, gall
             preg_prim_split(g,kmin-1,l) = preg_prim_split(g,kmin,l)
             preg_prim_split(g,kmax+1,l) = preg_prim_split(g,kmax,l)
          enddo
          !$omp end do nowait

          !$omp do
          do g = 1, gall
             PROG_split(g,kmin-1,l,I_RHOGE) = PROG_split(g,kmin,l,I_RHOGE)
             PROG_split(g,kmax+1,l,I_RHOGE) = PROG_split(g,kmax,l,I_RHOGE)
          enddo
          !$omp end do
       enddo
       !$omp end parallel

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                preg_prim_split_pl(g,k,l) = PROG_split_pl(g,k,l,I_RHOGE) * RovCV
             enddo
             enddo
             do g = 1, ADM_gall_pl
                preg_prim_split_pl(g,ADM_kmin-1,l) = preg_prim_split_pl(g,ADM_kmin,l)
                preg_prim_split_pl(g,ADM_kmax+1,l) = preg_prim_split_pl(g,ADM_kmax,l)
             enddo

             do g = 1, ADM_gall_pl
                PROG_split_pl(g,ADM_kmin-1,l,I_RHOGE) = PROG_split_pl(g,ADM_kmin,l,I_RHOGE)
                PROG_split_pl(g,ADM_kmax+1,l,I_RHOGE) = PROG_split_pl(g,ADM_kmax,l,I_RHOGE)
             enddo
          enddo
       endif

       if ( TIME_SPLIT ) then

          !---< Calculation of source term for Vh(vx,vy,vz) and W (split) >

          ! divergence damping
          call numfilter_divdamp( PROG_split(:,:,:,I_RHOGVX), PROG_split_pl(:,:,:,I_RHOGVX), & ! [IN]
                                  PROG_split(:,:,:,I_RHOGVY), PROG_split_pl(:,:,:,I_RHOGVY), & ! [IN]
                                  PROG_split(:,:,:,I_RHOGVZ), PROG_split_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                  PROG_split(:,:,:,I_RHOGW),  PROG_split_pl(:,:,:,I_RHOGW),  & ! [IN]
                                  ddivdvx   (:,:,:),          ddivdvx_pl   (:,:,:),          & ! [OUT]
                                  ddivdvy   (:,:,:),          ddivdvy_pl   (:,:,:),          & ! [OUT]
                                  ddivdvz   (:,:,:),          ddivdvz_pl   (:,:,:),          & ! [OUT]
                                  ddivdw    (:,:,:),          ddivdw_pl    (:,:,:)           ) ! [OUT]

          ! 2d divergence damping
          call numfilter_divdamp_2d( PROG_split(:,:,:,I_RHOGVX), PROG_split_pl(:,:,:,I_RHOGVX), & ! [IN]
                                     PROG_split(:,:,:,I_RHOGVY), PROG_split_pl(:,:,:,I_RHOGVY), & ! [IN]
                                     PROG_split(:,:,:,I_RHOGVZ), PROG_split_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                     ddivdvx_2d(:,:,:),          ddivdvx_2d_pl(:,:,:),          & ! [OUT]
                                     ddivdvy_2d(:,:,:),          ddivdvy_2d_pl(:,:,:),          & ! [OUT]
                                     ddivdvz_2d(:,:,:),          ddivdvz_2d_pl(:,:,:)           ) ! [OUT]

          ! pressure force
          ! dpgradw=0.0_RP becaude of f_type='HORIZONTAL'.
          call src_pres_gradient( preg_prim_split(:,:,:),   preg_prim_split_pl(:,:,:),   & ! [IN]
                                  dpgrad         (:,:,:,:), dpgrad_pl         (:,:,:,:), & ! [OUT]
                                  dpgradw        (:,:,:),   dpgradw_pl        (:,:,:),   & ! [OUT] not used
                                  I_SRC_horizontal                                       ) ! [IN]

          ! buoyancy force
          ! not calculated, because this term is implicit.

          !---< sum of tendencies ( large step + split{ pres-grad + div-damp + div-damp_2d } ) >
!OCL XFILL
          !$omp parallel do default(none),private(g,k,l,drhogvx,drhogvy,drhogvz),                  &
          !$omp shared(gall,kall,lall,drhogw,diff_vh,PROG_split,dt,g_TEND,                         &
          !$omp        dpgrad,ddivdvx,ddivdvy,ddivdvz,ddivdw,ddivdvx_2d,ddivdvy_2d,ddivdvz_2d,alpha), &
          !$omp collapse(2)
          do l = 1, lall
          do k = 1, kall
          do g = 1, gall
             drhogvx       = g_TEND(g,k,l,I_RHOGVX) - dpgrad(g,k,l,XDIR) + ddivdvx(g,k,l) + ddivdvx_2d(g,k,l)
             drhogvy       = g_TEND(g,k,l,I_RHOGVY) - dpgrad(g,k,l,YDIR) + ddivdvy(g,k,l) + ddivdvy_2d(g,k,l)
             drhogvz       = g_TEND(g,k,l,I_RHOGVZ) - dpgrad(g,k,l,ZDIR) + ddivdvz(g,k,l) + ddivdvz_2d(g,k,l)
             drhogw(g,k,l) = g_TEND(g,k,l,I_RHOGW)                       + ddivdw (g,k,l) * alpha

             diff_vh(g,k,l,1) = PROG_split(g,k,l,I_RHOGVX) + drhogvx * dt
             diff_vh(g,k,l,2) = PROG_split(g,k,l,I_RHOGVY) + drhogvy * dt
             diff_vh(g,k,l,3) = PROG_split(g,k,l,I_RHOGVZ) + drhogvz * dt
          enddo
          enddo
          enddo
          !$omp end parallel do

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                drhogvx          = g_TEND_pl(g,k,l,I_RHOGVX) - dpgrad_pl(g,k,l,XDIR) + ddivdvx_pl(g,k,l) + ddivdvx_2d_pl(g,k,l)
                drhogvy          = g_TEND_pl(g,k,l,I_RHOGVY) - dpgrad_pl(g,k,l,YDIR) + ddivdvy_pl(g,k,l) + ddivdvy_2d_pl(g,k,l)
                drhogvz          = g_TEND_pl(g,k,l,I_RHOGVZ) - dpgrad_pl(g,k,l,ZDIR) + ddivdvz_pl(g,k,l) + ddivdvz_2d_pl(g,k,l)
                drhogw_pl(g,k,l) = g_TEND_pl(g,k,l,I_RHOGW)                          + ddivdw_pl (g,k,l) * alpha

                diff_vh_pl(g,k,l,1) = PROG_split_pl(g,k,l,I_RHOGVX) + drhogvx * dt
                diff_vh_pl(g,k,l,2) = PROG_split_pl(g,k,l,I_RHOGVY) + drhogvy * dt
                diff_vh_pl(g,k,l,3) = PROG_split_pl(g,k,l,I_RHOGVZ) + drhogvz * dt
             enddo
             enddo
             enddo
          endif

       else ! NO-SPLITING

          !---< sum of tendencies ( large step ) >
!OCL XFILL
          !$omp parallel do default(none),private(g,k,l,drhogvx,drhogvy,drhogvz), &
          !$omp shared(gall,kall,lall,drhogw,diff_vh,PROG_split,dt,g_TEND), &
          !$omp collapse(2)
          do l = 1, lall
          do k = 1, kall
          do g = 1, gall
             drhogvx       = g_TEND(g,k,l,I_RHOGVX)
             drhogvy       = g_TEND(g,k,l,I_RHOGVY)
             drhogvz       = g_TEND(g,k,l,I_RHOGVZ)
             drhogw(g,k,l) = g_TEND(g,k,l,I_RHOGW)

             diff_vh(g,k,l,1) = PROG_split(g,k,l,I_RHOGVX) + drhogvx * dt
             diff_vh(g,k,l,2) = PROG_split(g,k,l,I_RHOGVY) + drhogvy * dt
             diff_vh(g,k,l,3) = PROG_split(g,k,l,I_RHOGVZ) + drhogvz * dt
          enddo
          enddo
          enddo
          !$omp end parallel do

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                drhogvx          = g_TEND_pl(g,k,l,I_RHOGVX)
                drhogvy          = g_TEND_pl(g,k,l,I_RHOGVY)
                drhogvz          = g_TEND_pl(g,k,l,I_RHOGVZ)
                drhogw_pl(g,k,l) = g_TEND_pl(g,k,l,I_RHOGW)

                diff_vh_pl(g,k,l,1) = PROG_split_pl(g,k,l,I_RHOGVX) + drhogvx * dt
                diff_vh_pl(g,k,l,2) = PROG_split_pl(g,k,l,I_RHOGVY) + drhogvy * dt
                diff_vh_pl(g,k,l,3) = PROG_split_pl(g,k,l,I_RHOGVZ) + drhogvz * dt
             enddo
             enddo
             enddo
          endif

       endif ! Split/Non-split

       ! treatment for boundary condition
       call BNDCND_rhovxvyvz( ADM_gall,              & ! [IN]
                              ADM_kall,              & ! [IN]
                              ADM_lall,              & ! [IN]
                              PROG   (:,:,:,I_RHOG), & ! [IN]
                              diff_vh(:,:,:,1),      & ! [INOUT]
                              diff_vh(:,:,:,2),      & ! [INOUT]
                              diff_vh(:,:,:,3)       ) ! [INOUT]

       if ( ADM_have_pl ) then
          call BNDCND_rhovxvyvz( ADM_gall_pl,              & ! [IN]
                                 ADM_kall,                 & ! [IN]
                                 ADM_lall_pl,              & ! [IN]
                                 PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                 diff_vh_pl(:,:,:,1),      & ! [INOUT]
                                 diff_vh_pl(:,:,:,2),      & ! [INOUT]
                                 diff_vh_pl(:,:,:,3)       ) ! [INOUT]
       endif

       call COMM_data_transfer( diff_vh, diff_vh_pl )

       call PROF_rapend  ('____vi_path1',2)
       call PROF_rapstart('____vi_path2',2)

       !---< vertical implicit scheme >

       call vi_main( diff_we        (:,:,:,1),        diff_we_pl        (:,:,:,1),        & ! [OUT]
                     diff_we        (:,:,:,2),        diff_we_pl        (:,:,:,2),        & ! [OUT]
                     diff_we        (:,:,:,3),        diff_we_pl        (:,:,:,3),        & ! [OUT]
                     diff_vh        (:,:,:,1),        diff_vh_pl        (:,:,:,1),        & ! [IN]
                     diff_vh        (:,:,:,2),        diff_vh_pl        (:,:,:,2),        & ! [IN]
                     diff_vh        (:,:,:,3),        diff_vh_pl        (:,:,:,3),        & ! [IN]
                     PROG_split     (:,:,:,I_RHOG),   PROG_split_pl     (:,:,:,I_RHOG),   & ! [IN]
                     PROG_split     (:,:,:,I_RHOGVX), PROG_split_pl     (:,:,:,I_RHOGVX), & ! [IN]
                     PROG_split     (:,:,:,I_RHOGVY), PROG_split_pl     (:,:,:,I_RHOGVY), & ! [IN]
                     PROG_split     (:,:,:,I_RHOGVZ), PROG_split_pl     (:,:,:,I_RHOGVZ), & ! [IN]
                     PROG_split     (:,:,:,I_RHOGW),  PROG_split_pl     (:,:,:,I_RHOGW),  & ! [IN]
                     PROG_split     (:,:,:,I_RHOGE),  PROG_split_pl     (:,:,:,I_RHOGE),  & ! [IN]
                     preg_prim_split(:,:,:),          preg_prim_split_pl(:,:,:),          & ! [IN]
                     PROG           (:,:,:,I_RHOG),   PROG_pl           (:,:,:,I_RHOG),   & ! [IN]
                     PROG           (:,:,:,I_RHOGVX), PROG_pl           (:,:,:,I_RHOGVX), & ! [IN]
                     PROG           (:,:,:,I_RHOGVY), PROG_pl           (:,:,:,I_RHOGVY), & ! [IN]
                     PROG           (:,:,:,I_RHOGVZ), PROG_pl           (:,:,:,I_RHOGVZ), & ! [IN]
                     PROG           (:,:,:,I_RHOGW),  PROG_pl           (:,:,:,I_RHOGW),  & ! [IN]
                     eth            (:,:,:),          eth_pl            (:,:,:),          & ! [IN]
                     g_TEND         (:,:,:,I_RHOG),   g_TEND_pl         (:,:,:,I_RHOG),   & ! [IN]
                     drhogw         (:,:,:),          drhogw_pl         (:,:,:),          & ! [IN]
                     g_TEND         (:,:,:,I_RHOGE),  g_TEND_pl         (:,:,:,I_RHOGE),  & ! [IN]
                     grhogetot0     (:,:,:),          grhogetot0_pl     (:,:,:),          & ! [IN]
                     dt                                     ) ! [IN]

       ! treatment for boundary condition
       call COMM_data_transfer( diff_we, diff_we_pl )

       ! update split value and mean mass flux
!OCL XFILL
       !$omp parallel do default(none),private(g,k,l) &
       !$omp shared(gall,kall,lall,PROG_split,diff_vh,diff_we), &
       !$omp collapse(2)
       do l = 1, lall
       do k = 1, kall
       do g = 1, gall
          PROG_split(g,k,l,I_RHOGVX) = diff_vh(g,k,l,1)
          PROG_split(g,k,l,I_RHOGVY) = diff_vh(g,k,l,2)
          PROG_split(g,k,l,I_RHOGVZ) = diff_vh(g,k,l,3)
          PROG_split(g,k,l,I_RHOG)   = diff_we(g,k,l,1)
          PROG_split(g,k,l,I_RHOGW)  = diff_we(g,k,l,2)
          PROG_split(g,k,l,I_RHOGE)  = diff_we(g,k,l,3)
       enddo
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(g,k,l,iv) &
       !$omp shared(gall,kall,lall,PROG_mean,PROG_split,rweight_itr), &
       !$omp collapse(3)
       do iv = I_RHOG, I_RHOGW
       do l  = 1, lall
       do k  = 1, kall
       do g  = 1, gall
          PROG_mean(g,k,l,iv) = PROG_mean(g,k,l,iv) + PROG_split(g,k,l,iv) * rweight_itr
       enddo
       enddo
       enddo
       enddo
       !$omp end parallel do

       if ( ADM_have_pl ) then
          PROG_split_pl(:,:,:,I_RHOGVX) = diff_vh_pl(:,:,:,1)
          PROG_split_pl(:,:,:,I_RHOGVY) = diff_vh_pl(:,:,:,2)
          PROG_split_pl(:,:,:,I_RHOGVZ) = diff_vh_pl(:,:,:,3)
          PROG_split_pl(:,:,:,I_RHOG)   = diff_we_pl(:,:,:,1)
          PROG_split_pl(:,:,:,I_RHOGW)  = diff_we_pl(:,:,:,2)
          PROG_split_pl(:,:,:,I_RHOGE)  = diff_we_pl(:,:,:,3)

          PROG_mean_pl(:,:,:,I_RHOG:I_RHOGW) = PROG_mean_pl (:,:,:,I_RHOG:I_RHOGW) &
                                             + PROG_split_pl(:,:,:,I_RHOG:I_RHOGW) * rweight_itr
       endif

       call PROF_rapend  ('____vi_path2',2)

    enddo ! small step end
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------
    call PROF_rapstart('____vi_path3',2)

    ! update prognostic variables

    !$omp parallel do default(none),private(g,k,l,iv) &
    !$omp shared(gall,kall,lall,PROG,PROG_split,rweight_itr), &
    !$omp collapse(3)
    do iv = I_RHOG, I_RHOGE
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       PROG(g,k,l,iv) = PROG(g,k,l,iv) + PROG_split(g,k,l,iv)
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       PROG_pl(:,:,:,:) = PROG_pl(:,:,:,:) + PROG_split_pl(:,:,:,:)
    endif

    call OPRT_horizontalize_vec( PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                                 PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                                 PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ)  ) ! [INOUT]

    ! communication of mean velocity
    call COMM_data_transfer( PROG_mean, PROG_mean_pl )

    call PROF_rapend  ('____vi_path3',2)

    return
  end subroutine vi_small_step

  !-----------------------------------------------------------------------------
  !> Main part of the vertical implicit scheme
  subroutine vi_main( &
       rhog_split1,      rhog_split1_pl,      &
       rhogw_split1,     rhogw_split1_pl,     &
       rhoge_split1,     rhoge_split1_pl,     &
       rhogvx_split1,    rhogvx_split1_pl,    &
       rhogvy_split1,    rhogvy_split1_pl,    &
       rhogvz_split1,    rhogvz_split1_pl,    &
       rhog_split0,      rhog_split0_pl,      &
       rhogvx_split0,    rhogvx_split0_pl,    &
       rhogvy_split0,    rhogvy_split0_pl,    &
       rhogvz_split0,    rhogvz_split0_pl,    &
       rhogw_split0,     rhogw_split0_pl,     &
       rhoge_split0,     rhoge_split0_pl,     &
       preg_prim_split0, preg_prim_split0_pl, &
       rhog0,            rhog0_pl,            &
       rhogvx0,          rhogvx0_pl,          &
       rhogvy0,          rhogvy0_pl,          &
       rhogvz0,          rhogvz0_pl,          &
       rhogw0,           rhogw0_pl,           &
       eth0,             eth0_pl,             &
       grhog,            grhog_pl,            &
       grhogw,           grhogw_pl,           &
       grhoge,           grhoge_pl,           &
       grhogetot,        grhogetot_pl,        &
       dt                                     )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_const, only: &
       CONST_Rdry, &
       CONST_CVdry
    use mod_vmtr, only: &
       VMTR_PHI,          &
       VMTR_PHI_pl,       &
       VMTR_C2WfactGz,    &
       VMTR_C2WfactGz_pl
    use mod_time, only: &
       TIME_SPLIT
    use mod_bndcnd, only: &
       BNDCND_rhow
    use mod_cnvvar, only: &
       cnvvar_rhogkin
    use mod_src, only: &
       src_flux_convergence,      &
       src_advection_convergence, &
       I_SRC_horizontal,          &
       I_SRC_default
    implicit none

    real(RP), intent(out) :: rhog_split1        (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars (split, at n+1 step)
    real(RP), intent(out) :: rhog_split1_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogw_split1       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogw_split1_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhoge_split1       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhoge_split1_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx_split1      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx_split1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy_split1      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy_split1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz_split1      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz_split1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: rhog_split0        (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars (split, at n step)
    real(RP), intent(in)  :: rhog_split0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx_split0      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx_split0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy_split0      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy_split0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz_split0      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz_split0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw_split0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogw_split0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhoge_split0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhoge_split0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: preg_prim_split0   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: preg_prim_split0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: rhog0              (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars ( previous )
    real(RP), intent(in)  :: rhog0_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx0            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx0_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy0            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy0_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz0            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz0_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw0             (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogw0_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: eth0               (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: eth0_pl            (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: grhog              (ADM_gall   ,ADM_kall,ADM_lall   ) ! large step tendency
    real(RP), intent(in)  :: grhog_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: grhogw             (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: grhogw_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: grhoge             (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: grhoge_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: grhogetot          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: grhogetot_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: dt

    real(RP) :: drhog       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term at t=n+1
    real(RP) :: drhog_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhoge      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhoge_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhogetot   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhogetot_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: grhog1      (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term ( large step + t=n+1 )
    real(RP) :: grhog1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: grhoge1     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhoge1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: gpre        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: gpre_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhog1       (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars ( previous + t=n,t=n+1 )
    real(RP) :: rhog1_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx1     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy1     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz1     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw1      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhogkin0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! kinetic energy ( previous                )
    real(RP) :: rhogkin0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogkin10   (ADM_gall   ,ADM_kall,ADM_lall   ) ! kinetic energy ( previous + split(t=n)   )
    real(RP) :: rhogkin10_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogkin11   (ADM_gall   ,ADM_kall,ADM_lall   ) ! kinetic energy ( previous + split(t=n+1) )
    real(RP) :: rhogkin11_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ethtot0     (ADM_gall   ,ADM_kall,ADM_lall   ) ! total enthalpy ( h + v^{2}/2 + phi, previous )
    real(RP) :: ethtot0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer  :: gall, kall, lall
    real(RP) :: Rdry, CVdry

    integer  :: g, k, l
    !---------------------------------------------------------------------------

    gall = ADM_gall
    kall = ADM_kall
    lall = ADM_lall

    Rdry  = CONST_Rdry
    CVdry = CONST_CVdry

    !---< update grhog & grhoge >

    if ( TIME_SPLIT ) then

       ! horizontal flux convergence
       call src_flux_convergence( rhogvx_split1, rhogvx_split1_pl, & ! [IN]
                                  rhogvy_split1, rhogvy_split1_pl, & ! [IN]
                                  rhogvz_split1, rhogvz_split1_pl, & ! [IN]
                                  rhogw_split0,  rhogw_split0_pl,  & ! [IN]
                                  drhog,         drhog_pl,         & ! [OUT]
                                  I_SRC_horizontal                 ) ! [IN]

       ! horizontal advection convergence
       call src_advection_convergence( rhogvx_split1, rhogvx_split1_pl, & ! [IN]
                                       rhogvy_split1, rhogvy_split1_pl, & ! [IN]
                                       rhogvz_split1, rhogvz_split1_pl, & ! [IN]
                                       rhogw_split0,  rhogw_split0_pl,  & ! [IN]
                                       eth0,          eth0_pl,          & ! [IN]
                                       drhoge,        drhoge_pl,        & ! [OUT]
                                       I_SRC_horizontal                 ) ! [IN]

    else

!OCL XFILL
       !$omp parallel do default(none),private(g,k,l) &
       !$omp shared(gall,kall,lall,drhog,drhoge), &
       !$omp collapse(2)
       do l  = 1, lall
       do k  = 1, kall
       do g  = 1, gall
          drhog (g,k,l) = 0.0_RP
          drhoge(g,k,l) = 0.0_RP
       enddo
       enddo
       enddo
       !$omp end parallel do

       drhog_pl (:,:,:) = 0.0_RP
       drhoge_pl(:,:,:) = 0.0_RP

    endif

    ! update grhog, grhoge and calc source term of pressure
!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,grhog1,grhoge1,gpre,grhog,grhoge,drhog,drhoge,Rdry,CVdry), &
    !$omp collapse(2)
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       grhog1 (g,k,l) = grhog  (g,k,l) + drhog (g,k,l)
       grhoge1(g,k,l) = grhoge (g,k,l) + drhoge(g,k,l)
       gpre   (g,k,l) = grhoge1(g,k,l) * Rdry / CVdry
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       grhog1_pl (:,:,:) = grhog_pl  (:,:,:) + drhog_pl (:,:,:)
       grhoge1_pl(:,:,:) = grhoge_pl (:,:,:) + drhoge_pl(:,:,:)
       gpre_pl   (:,:,:) = grhoge1_pl(:,:,:) * Rdry / CVdry
    endif

    !---------------------------------------------------------------------------
    ! verical implict calculation core
    !---------------------------------------------------------------------------

    ! boundary condition for rhogw_split1
!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,rhogw_split1), &
    !$omp collapse(2)
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       rhogw_split1(g,k,l) = 0.0_RP
    enddo
    enddo
    enddo
    !$omp end parallel do

    do l = 1, ADM_lall
       call BNDCND_rhow( ADM_gall,               & ! [IN]
                         rhogvx_split1 (:,:,l),  & ! [IN]
                         rhogvy_split1 (:,:,l),  & ! [IN]
                         rhogvz_split1 (:,:,l),  & ! [IN]
                         rhogw_split1  (:,:,l),  & ! [INOUT]
                         VMTR_C2WfactGz(:,:,:,l) ) ! [IN]
    enddo

    if ( ADM_have_pl ) then
       rhogw_split1_pl(:,:,:) = 0.0_RP

       do l = 1, ADM_lall_pl
          call BNDCND_rhow( ADM_gall_pl,               & ! [IN]
                            rhogvx_split1_pl (:,:,l),  & ! [IN]
                            rhogvy_split1_pl (:,:,l),  & ! [IN]
                            rhogvz_split1_pl (:,:,l),  & ! [IN]
                            rhogw_split1_pl  (:,:,l),  & ! [INOUT]
                            VMTR_C2WfactGz_pl(:,:,:,l) ) ! [IN]
       enddo
    endif

    ! update rhogw_split1
    call vi_rhow_solver( rhogw_split1,     rhogw_split1_pl,     & ! [INOUT]
                         rhogw_split0,     rhogw_split0_pl,     & ! [IN]
                         preg_prim_split0, preg_prim_split0_pl, & ! [IN]
                         rhog_split0,      rhog_split0_pl,      & ! [IN]
                         grhog1,           grhog1_pl,           & ! [IN]
                         grhogw,           grhogw_pl,           & ! [IN]
                         gpre,             gpre_pl,             & ! [IN]
                         dt                                     ) ! [IN]

    ! update rhog_split1
    call src_flux_convergence( rhogvx_split1, rhogvx_split1_pl, & ! [IN]
                               rhogvy_split1, rhogvy_split1_pl, & ! [IN]
                               rhogvz_split1, rhogvz_split1_pl, & ! [IN]
                               rhogw_split1,  rhogw_split1_pl,  & ! [IN]
                               drhog,         drhog_pl,         & ! [OUT]
                               I_SRC_default                    ) ! [IN]

!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,rhog_split1,rhog_split0,grhog,drhog,dt), &
    !$omp collapse(2)
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       rhog_split1(g,k,l) = rhog_split0(g,k,l) + ( grhog(g,k,l) + drhog(g,k,l) ) * dt
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       rhog_split1_pl(:,:,:) = rhog_split0_pl(:,:,:) + ( grhog_pl(:,:,:) + drhog_pl(:,:,:) ) * dt
    endif

    !---------------------------------------------------------------------------
    ! energy correction by Etotal (Satoh,2002)
    !---------------------------------------------------------------------------

    ! calc rhogkin ( previous )
    call cnvvar_rhogkin( rhog0,    rhog0_pl,   & ! [IN]
                         rhogvx0,  rhogvx0_pl, & ! [IN]
                         rhogvy0,  rhogvy0_pl, & ! [IN]
                         rhogvz0,  rhogvz0_pl, & ! [IN]
                         rhogw0,   rhogw0_pl,  & ! [IN]
                         rhogkin0, rhogkin0_pl ) ! [OUT]

    ! prognostic variables ( previous + split (t=n) )
!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,rhog1,rhogvx1,rhogvy1,rhogvz1,rhogw1,              &
    !$omp        rhog0,rhogvx0,rhogvy0,rhogvz0,rhogw0,                             &
    !$omp        rhog_split0,rhogvx_split0,rhogvy_split0,rhogvz_split0,rhogw_split0), &
    !$omp collapse(2)
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       rhog1  (g,k,l) = rhog0  (g,k,l) + rhog_split0  (g,k,l)
       rhogvx1(g,k,l) = rhogvx0(g,k,l) + rhogvx_split0(g,k,l)
       rhogvy1(g,k,l) = rhogvy0(g,k,l) + rhogvy_split0(g,k,l)
       rhogvz1(g,k,l) = rhogvz0(g,k,l) + rhogvz_split0(g,k,l)
       rhogw1 (g,k,l) = rhogw0 (g,k,l) + rhogw_split0 (g,k,l)
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       rhog1_pl  (:,:,:) = rhog0_pl  (:,:,:) + rhog_split0_pl  (:,:,:)
       rhogvx1_pl(:,:,:) = rhogvx0_pl(:,:,:) + rhogvx_split0_pl(:,:,:)
       rhogvy1_pl(:,:,:) = rhogvy0_pl(:,:,:) + rhogvy_split0_pl(:,:,:)
       rhogvz1_pl(:,:,:) = rhogvz0_pl(:,:,:) + rhogvz_split0_pl(:,:,:)
       rhogw1_pl (:,:,:) = rhogw0_pl (:,:,:) + rhogw_split0_pl (:,:,:)
    endif

    ! calc rhogkin ( previous + split(t=n) )
    call cnvvar_rhogkin( rhog1,     rhog1_pl,    & ! [IN]
                         rhogvx1,   rhogvx1_pl,  & ! [IN]
                         rhogvy1,   rhogvy1_pl,  & ! [IN]
                         rhogvz1,   rhogvz1_pl,  & ! [IN]
                         rhogw1,    rhogw1_pl,   & ! [IN]
                         rhogkin10, rhogkin10_pl ) ! [OUT]

    ! prognostic variables ( previous + split (t=n+1) )
!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,rhog1,rhogvx1,rhogvy1,rhogvz1,rhogw1,              &
    !$omp        rhog0,rhogvx0,rhogvy0,rhogvz0,rhogw0,                             &
    !$omp        rhog_split1,rhogvx_split1,rhogvy_split1,rhogvz_split1,rhogw_split1), &
    !$omp collapse(2)
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       rhog1  (g,k,l) = rhog0  (g,k,l) + rhog_split1  (g,k,l)
       rhogvx1(g,k,l) = rhogvx0(g,k,l) + rhogvx_split1(g,k,l)
       rhogvy1(g,k,l) = rhogvy0(g,k,l) + rhogvy_split1(g,k,l)
       rhogvz1(g,k,l) = rhogvz0(g,k,l) + rhogvz_split1(g,k,l)
       rhogw1 (g,k,l) = rhogw0 (g,k,l) + rhogw_split1 (g,k,l)
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       rhog1_pl  (:,:,:) = rhog0_pl  (:,:,:) + rhog_split1_pl  (:,:,:)
       rhogvx1_pl(:,:,:) = rhogvx0_pl(:,:,:) + rhogvx_split1_pl(:,:,:)
       rhogvy1_pl(:,:,:) = rhogvy0_pl(:,:,:) + rhogvy_split1_pl(:,:,:)
       rhogvz1_pl(:,:,:) = rhogvz0_pl(:,:,:) + rhogvz_split1_pl(:,:,:)
       rhogw1_pl (:,:,:) = rhogw0_pl (:,:,:) + rhogw_split1_pl (:,:,:)
    endif

    ! calc rhogkin ( previous + split(t=n+1) )
    call cnvvar_rhogkin( rhog1,     rhog1_pl,    & ! [IN]
                         rhogvx1,   rhogvx1_pl,  & ! [IN]
                         rhogvy1,   rhogvy1_pl,  & ! [IN]
                         rhogvz1,   rhogvz1_pl,  & ! [IN]
                         rhogw1,    rhogw1_pl,   & ! [IN]
                         rhogkin11, rhogkin11_pl ) ! [OUT]

    ! calculate total enthalpy ( h + v^{2}/2 + phi, previous )
!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,ethtot0,eth0,rhogkin0,rhog0,VMTR_PHI), &
    !$omp collapse(2)
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       ethtot0(g,k,l) = eth0    (g,k,l)                &
                      + rhogkin0(g,k,l) / rhog0(g,k,l) &
                      + VMTR_PHI(g,k,l)
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       ethtot0_pl(:,:,:) = eth0_pl(:,:,:)                       &
                         + rhogkin0_pl(:,:,:) / rhog0_pl(:,:,:) &
                         + VMTR_PHI_pl(:,:,:)
    endif

    ! advection convergence for eth + kin + phi
    call src_advection_convergence( rhogvx1,   rhogvx1_pl,   & ! [IN]
                                    rhogvy1,   rhogvy1_pl,   & ! [IN]
                                    rhogvz1,   rhogvz1_pl,   & ! [IN]
                                    rhogw1,    rhogw1_pl,    & ! [IN]
                                    ethtot0,   ethtot0_pl,   & ! [IN]
                                    drhogetot, drhogetot_pl, & ! [OUT]
                                    I_SRC_default            ) ! [IN]

!OCL XFILL
    !$omp parallel do default(none),private(g,k,l), &
    !$omp shared(gall,kall,lall,rhoge_split1,rhoge_split0,grhogetot,drhogetot,dt, &
    !$omp        rhogkin10,rhogkin11,rhog_split0,rhog_split1,VMTR_PHI), &
    !$omp collapse(2)
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       rhoge_split1(g,k,l) = rhoge_split0 (g,k,l)                     & ! t=n
                           + ( grhogetot  (g,k,l)                     & ! tendency of total energy (num.diff+smg+nudge)
                             + drhogetot  (g,k,l) ) * dt              & ! tendency of total energy (adv.conv.)
                           + ( rhogkin10  (g,k,l)                     & ! kinetic   energy (t=n)
                             - rhogkin11  (g,k,l) )                   & ! kinetic   energy (t=n+1)
                           + ( rhog_split0(g,k,l)                     & ! potential energy (diff,t=n)
                             - rhog_split1(g,k,l) ) * VMTR_PHI(g,k,l)   ! potential energy (diff,t=n+1)
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       rhoge_split1_pl(:,:,:) = rhoge_split0_pl (:,:,:)                        &
                              + ( grhogetot_pl  (:,:,:)                        &
                                + drhogetot_pl  (:,:,:) ) * dt                 &
                              + ( rhogkin10_pl  (:,:,:)                        &
                                - rhogkin11_pl  (:,:,:) )                      &
                              + ( rhog_split0_pl(:,:,:)                        &
                                - rhog_split1_pl(:,:,:) ) * VMTR_PHI_pl(:,:,:)
    endif

    return
  end subroutine vi_main

  !-----------------------------------------------------------------------------
  !> Update tridiagonal matrix
  subroutine vi_rhow_update_matrix( &
       eth,     eth_pl,     &
       g_tilde, g_tilde_pl, &
       dt                   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_const, only: &
       GRAV  => CONST_GRAV, &
       Rdry  => CONST_Rdry, &
       CVdry => CONST_CVdry
    use mod_grd, only: &
       GRD_rdgzh, &
       GRD_rdgz,  &
       GRD_cfact, &
       GRD_dfact
    use mod_vmtr, only: &
       VMTR_GAM2H,      &
       VMTR_GAM2H_pl,   &
       VMTR_RGSQRTH,    &
       VMTR_RGSQRTH_pl, &
       VMTR_RGAMH,      &
       VMTR_RGAMH_pl,   &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    implicit none

    real(RP), intent(in) :: eth       (ADM_gall   ,ADM_kall,ADM_lall   ) ! enthalpy at the half lev
    real(RP), intent(in) :: eth_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: g_tilde   (ADM_gall   ,ADM_kall,ADM_lall   ) ! effective gravitation at the half lev
    real(RP), intent(in) :: g_tilde_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: dt

    integer  :: gall, kmin, kmax, lall
    real(RP) :: GCVovR   ! g * Cv / R
    real(RP) :: ACVovRt2 ! alpha * Cv / R / dt**2

    integer  :: g, k, l
    !---------------------------------------------------------------------------
    ! Original concept
    !
    ! A_o(:,:,:) = VMTR_RGSGAM2(:,:,:)
    ! A_i(:,:,:) = VMTR_GAM2H(:,:,:) * eth(:,:,:) ! [debug] 20120727 H.Yashiro
    ! B  (:,:,:) = g_tilde(:,:,:)
    ! C_o(:,:,:) = VMTR_RGAM2H (:,:,:) * ( CONST_CVdry / CONST_Rdry * CONST_GRAV )
    ! C_i(:,:,:) = 1.0_RP / VMTR_RGAM2H(:,:,:)
    ! D  (:,:,:) = CONST_CVdry / CONST_Rdry / ( dt*dt ) / VMTR_RGSQRTH(:,:,:)

    ! do k = ADM_kmin+1, ADM_kmax
    !    Mc(:,k,:) = dble(NON_HYDRO_ALPHA) *D(:,k,:)              &
    !              + GRD_rdgzh(k)                                 &
    !              * ( GRD_rdgz (k)   * A_o(:,k  ,:) * A_i(:,k,:) &
    !                + GRD_rdgz (k-1) * A_o(:,k-1,:) * A_i(:,k,:) &
    !                - 0.5_RP * ( GRD_dfact(k) - GRD_cfact(k-1) ) &
    !                * ( B(:,k,:) + C_o(:,k,:) * C_i(:,k,:) )     &
    !                )
    !    Mu(:,k,:) = - GRD_rdgzh(k) * GRD_rdgz(k) * A_o(:,k,:) * A_i(:,k+1,:) &
    !                - GRD_rdgzh(k) * 0.5_RP * GRD_cfact(k)                   &
    !                * ( B(:,k+1,:) + C_o(:,k,:) * C_i(:,k+1,:) )
    !    Ml(:,k,:) = - GRD_rdgzh(k) * GRD_rdgz(k) * A_o(:,k,:) * A_i(:,k-1,:) &
    !                + GRD_rdgzh(k) * 0.5_RP * GRD_dfact(k-1)                 &
    !                * ( B(:,k-1,:) + C_o(:,k,:) * C_i(:,k-1,:) )
    ! enddo

    call PROF_rapstart('____vi_rhow_update_matrix',2)

    gall = ADM_gall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    GCVovR   = GRAV * CVdry / Rdry
    ACVovRt2 = real(NON_HYDRO_ALPHA,kind=RP) * CVdry / Rdry / ( dt*dt )

    !$omp parallel do default(none),private(g,k,l),                                               &
    !$omp shared(gall,kmin,kmax,lall,eth,g_tilde,dt,GCVovR,ACVovRt2,Mu,Mc,Ml,GRD_cfact,GRD_dfact, &
    !$omp        GRD_rdgzh,GRD_rdgz,VMTR_GAM2H,VMTR_RGSQRTH,VMTR_RGAMH,VMTR_RGSGAM2), &
    !$omp collapse(2)
    do l = 1, lall
    do k = kmin+1, kmax
!OCL XFILL
       do g = 1, gall
          Mc(g,k,l) = ACVovRt2 / VMTR_RGSQRTH(g,k,l)                             &
                    + GRD_rdgzh(k) * ( ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k  )   &
                                       + VMTR_RGSGAM2(g,k-1,l) * GRD_rdgz(k-1) ) &
                                       * VMTR_GAM2H  (g,k  ,l) * eth(g,k,l)      &
                                     - ( GRD_dfact(k) - GRD_cfact(k-1) )         &
                                     * ( g_tilde(g,k,l) + GCVovR )               )
       enddo

!OCL XFILL
       do g = 1, gall
          Mu(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k)                     &
                                      * VMTR_GAM2H  (g,k+1,l) * eth(g,k+1,l)                    &
                                      + GRD_cfact(k)                                            &
                                      * ( g_tilde   (g,k+1,l)                                   &
                                        + VMTR_GAM2H(g,k+1,l)* VMTR_RGAMH(g,k,l)**2 * GCVovR  ) )
       enddo

!OCL XFILL
       do g = 1, gall
          Ml(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k)                     &
                                      * VMTR_GAM2H  (g,k-1,l) * eth(g,k-1,l)                    &
                                      - GRD_dfact(k-1)                                          &
                                      * ( g_tilde   (g,k-1,l)                                   &
                                        + VMTR_GAM2H(g,k-1,l) * VMTR_RGAMH(g,k,l)**2 * GCVovR ) )
       enddo
    enddo
    enddo
    !$omp end parallel do

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Mc_pl(g,k,l) = ACVovRt2 / VMTR_RGSQRTH_pl(g,k,l)                             &
                          + GRD_rdgzh(k) * ( ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k  )   &
                                             + VMTR_RGSGAM2_pl(g,k-1,l) * GRD_rdgz(k-1) ) &
                                             * VMTR_GAM2H_pl  (g,k  ,l) * eth_pl(g,k,l)   &
                                           - ( GRD_dfact(k) - GRD_cfact(k-1) )            &
                                           * ( g_tilde_pl(g,k,l) + GCVovR )               )
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Mu_pl(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k)                        &
                                            * VMTR_GAM2H_pl  (g,k+1,l) * eth_pl(g,k+1,l)                    &
                                            + GRD_cfact(k)                                                  &
                                            * ( g_tilde_pl   (g,k+1,l)                                      &
                                              + VMTR_GAM2H_pl(g,k+1,l)* VMTR_RGAMH_pl(g,k,l)**2 * GCVovR  ) )
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Ml_pl(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k)                        &
                                            * VMTR_GAM2H_pl  (g,k-1,l) * eth_pl(g,k-1,l)                    &
                                            - GRD_dfact(k-1)                                                &
                                            * ( g_tilde_pl   (g,k-1,l)                                      &
                                              + VMTR_GAM2H_pl(g,k-1,l) * VMTR_RGAMH_pl(g,k,l)**2 * GCVovR ) )
          enddo
          enddo
       enddo
    endif

    call PROF_rapend('____vi_rhow_update_matrix',2)

    return
  end subroutine vi_rhow_update_matrix

  !-----------------------------------------------------------------------------
  !> Tridiagonal matrix solver
  subroutine vi_rhow_solver( &
       rhogw,  rhogw_pl,  &
       rhogw0, rhogw0_pl, &
       preg0,  preg0_pl,  &
       rhog0,  rhog0_pl,  &
       Srho,   Srho_pl,   &
       Sw,     Sw_pl,     &
       Spre,   Spre_pl,   &
       dt                 )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_const, only: &
       CONST_GRAV, &
       CONST_Rdry, &
       CONST_CVdry
    use mod_grd, only: &
       GRD_rdgzh, &
       GRD_afact, &
       GRD_bfact
    use mod_vmtr, only: &
       VMTR_GSGAM2H,     &
       VMTR_GSGAM2H_pl,  &
       VMTR_RGAM,        &
       VMTR_RGAM_pl,     &
       VMTR_RGAMH,       &
       VMTR_RGAMH_pl,    &
       VMTR_RGSGAM2,     &
       VMTR_RGSGAM2_pl,  &
       VMTR_RGSGAM2H,    &
       VMTR_RGSGAM2H_pl
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    !$ use omp_lib
    implicit none

    real(RP), intent(inout) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 ), n+1
    real(RP), intent(inout) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)    :: rhogw0   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho            ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Srho     (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rho  at the full level
    real(RP), intent(in)    :: Srho_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sw       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rhow at the half level
    real(RP), intent(in)    :: Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Spre     (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for pres at the full level
    real(RP), intent(in)    :: Spre_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: dt

    real(RP) :: Sall    (ADM_gall,   ADM_kall)
    real(RP) :: Sall_pl (ADM_gall_pl,ADM_kall)
    real(RP) :: beta    (ADM_gall   )
    real(RP) :: beta_pl (ADM_gall_pl)
    real(RP) :: gamma   (ADM_gall,   ADM_kall)
    real(RP) :: gamma_pl(ADM_gall_pl,ADM_kall)

    integer  :: gall, kmin, kmax, lall
    real(RP) :: grav
    real(RP) :: CVovRt2 ! Cv / R / dt**2
    real(RP) :: alpha

    integer  :: g, k, l
    integer  :: gstr, gend
    !$ integer  :: n_per_thread
    !$ integer  :: n_thread
    !---------------------------------------------------------------------------

    call PROF_rapstart('____vi_rhow_solver',2)

    gall = ADM_gall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    grav    = CONST_GRAV
    CVovRt2 = CONST_CVdry / CONST_Rdry / (dt*dt)
    alpha   = real(NON_HYDRO_ALPHA,kind=RP)

    !$omp parallel default(none),private(g,k,l), &
    !$omp private(gstr,gend,n_thread,n_per_thread) &
    !$omp shared(gall,kmin,kmax,lall,rhogw,rhogw0,preg0,rhog0,Srho,Sw,Spre,dt,Sall,beta,gamma,Mu,Mc,Ml, &
    !$omp        GRD_afact,GRD_bfact,GRD_rdgzh,VMTR_GSGAM2H,VMTR_RGAM,VMTR_RGAMH,VMTR_RGSGAM2,VMTR_RGSGAM2H,grav,alpha,CVovRt2)
    gstr = 1
    gend = gall
    !$ n_thread     = omp_get_num_threads()
    !$ n_per_thread = gall / n_thread + int( 0.5_RP + sign(0.5_RP,mod(gall,n_thread)-0.5_RP) )
    !$ gstr         = n_per_thread * omp_get_thread_num() + 1
    !$ gend         = min( gstr+n_per_thread-1, gall )

    do l = 1, lall
       ! calc Sall
       do k = kmin+1, kmax
       do g = gstr, gend
          Sall(g,k) = (   ( rhogw0(g,k,  l)*alpha + dt * Sw  (g,k,  l) ) * VMTR_RGAMH  (g,k,  l)**2             &
                      - ( ( preg0 (g,k,  l)       + dt * Spre(g,k,  l) ) * VMTR_RGSGAM2(g,k,  l)                &
                        - ( preg0 (g,k-1,l)       + dt * Spre(g,k-1,l) ) * VMTR_RGSGAM2(g,k-1,l)                &
                        ) * dt * GRD_rdgzh(k)                                                                   &
                      - ( ( rhog0 (g,k,  l)       + dt * Srho(g,k,  l) ) * VMTR_RGAM(g,k,  l)**2 * GRD_afact(k) &
                        + ( rhog0 (g,k-1,l)       + dt * Srho(g,k-1,l) ) * VMTR_RGAM(g,k-1,l)**2 * GRD_bfact(k) &
                        ) * dt * grav                                                                           &
                      ) * CVovRt2
       enddo
       enddo

       ! boundary conditions
       do g = gstr, gend
          rhogw(g,kmin,  l) = rhogw(g,kmin,  l) * VMTR_RGSGAM2H(g,kmin,  l)
          rhogw(g,kmax+1,l) = rhogw(g,kmax+1,l) * VMTR_RGSGAM2H(g,kmax+1,l)
          Sall (g,kmin+1)   = Sall (g,kmin+1) - Ml(g,kmin+1,l) * rhogw(g,kmin,  l)
          Sall (g,kmax  )   = Sall (g,kmax  ) - Mu(g,kmax,  l) * rhogw(g,kmax+1,l)
       enddo

       !---< solve tri-daigonal matrix >

       ! condition at kmin+1
       k = kmin+1
       do g = gstr, gend
          beta (g)     = Mc(g,k,l)
          rhogw(g,k,l) = Sall(g,k) / beta(g)
       enddo

       ! forward
       do k = kmin+2, kmax
       do g = gstr, gend
          gamma(g,k)   = Mu(g,k-1,l) / beta(g)
          beta (g)     = Mc(g,k,l) - Ml(g,k,l) * gamma(g,k) ! update beta
          rhogw(g,k,l) = ( Sall(g,k) - Ml(g,k,l) * rhogw(g,k-1,l) ) / beta(g)
       enddo
       enddo

       ! backward
       do k = kmax-1, kmin+1, -1
       do g = gstr, gend
          rhogw(g,k  ,l) = rhogw(g,k  ,l) - gamma(g,k+1) * rhogw(g,k+1,l)
          rhogw(g,k+1,l) = rhogw(g,k+1,l) * VMTR_GSGAM2H(g,k+1,l) ! return value ( G^1/2 x gam2 )
       enddo
       enddo

       ! boundary treatment
       do g = gstr, gend
          rhogw(g,kmin  ,l) = rhogw(g,kmin  ,l) * VMTR_GSGAM2H(g,kmin  ,l)
          rhogw(g,kmin+1,l) = rhogw(g,kmin+1,l) * VMTR_GSGAM2H(g,kmin+1,l)
          rhogw(g,kmax+1,l) = rhogw(g,kmax+1,l) * VMTR_GSGAM2H(g,kmax+1,l)
       enddo
    enddo
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k  = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Sall_pl(g,k) = (   ( rhogw0_pl(g,k,  l)*alpha + dt * Sw_pl  (g,k,  l) ) * VMTR_RGAMH_pl  (g,k,  l)**2             &
                            - ( ( preg0_pl (g,k,  l)       + dt * Spre_pl(g,k,  l) ) * VMTR_RGSGAM2_pl(g,k,  l)                &
                              - ( preg0_pl (g,k-1,l)       + dt * Spre_pl(g,k-1,l) ) * VMTR_RGSGAM2_pl(g,k-1,l)                &
                              ) * dt * GRD_rdgzh(k)                                                                            &
                            - ( ( rhog0_pl (g,k,  l)       + dt * Srho_pl(g,k,  l) ) * VMTR_RGAM_pl(g,k,  l)**2 * GRD_afact(k) &
                              + ( rhog0_pl (g,k-1,l)       + dt * Srho_pl(g,k-1,l) ) * VMTR_RGAM_pl(g,k-1,l)**2 * GRD_bfact(k) &
                              ) * dt * grav                                                                                    &
                            ) * CVovRt2
          enddo
          enddo

          do g = 1, ADM_gall_pl
             rhogw_pl(g,ADM_kmin,  l) = rhogw_pl(g,ADM_kmin,  l) * VMTR_RGSGAM2H_pl(g,ADM_kmin,  l)
             rhogw_pl(g,ADM_kmax+1,l) = rhogw_pl(g,ADM_kmax+1,l) * VMTR_RGSGAM2H_pl(g,ADM_kmax+1,l)
             Sall_pl (g,ADM_kmin+1)   = Sall_pl (g,ADM_kmin+1) - Ml_pl(g,ADM_kmin+1,l) * rhogw_pl(g,ADM_kmin,  l)
             Sall_pl (g,ADM_kmax  )   = Sall_pl (g,ADM_kmax  ) - Mu_pl(g,ADM_kmax,  l) * rhogw_pl(g,ADM_kmax+1,l)
          enddo

          k = ADM_kmin+1
          do g = 1, ADM_gall_pl
             beta_pl (g)     = Mc_pl(g,k,l)
             rhogw_pl(g,k,l) = Sall_pl(g,k) / beta_pl(g)
          enddo

          do k = ADM_kmin+2, ADM_kmax
          do g = 1, ADM_gall_pl
             gamma_pl(g,k)   = Mu_pl(g,k-1,l) / beta_pl(g)
             beta_pl (g)     = Mc_pl(g,k,l) - Ml_pl(g,k,l) * gamma_pl(g,k) ! update beta
             rhogw_pl(g,k,l) = ( Sall_pl(g,k) - Ml_pl(g,k,l) * rhogw_pl(g,k-1,l) ) / beta_pl(g)
          enddo
          enddo

          ! backward
          do k = ADM_kmax-1, ADM_kmin+1, -1
          do g = 1, ADM_gall_pl
             rhogw_pl(g,k  ,l) = rhogw_pl(g,k  ,l) - gamma_pl(g,k+1) * rhogw_pl(g,k+1,l)
             rhogw_pl(g,k+1,l) = rhogw_pl(g,k+1,l) * VMTR_GSGAM2H_pl(g,k+1,l) ! return value ( G^1/2 x gam2 )
          enddo
          enddo

          ! boundary treatment
          do g = 1, ADM_gall_pl
             rhogw_pl(g,ADM_kmin  ,l) = rhogw_pl(g,ADM_kmin  ,l) * VMTR_GSGAM2H_pl(g,ADM_kmin  ,l)
             rhogw_pl(g,ADM_kmin+1,l) = rhogw_pl(g,ADM_kmin+1,l) * VMTR_GSGAM2H_pl(g,ADM_kmin+1,l)
             rhogw_pl(g,ADM_kmax+1,l) = rhogw_pl(g,ADM_kmax+1,l) * VMTR_GSGAM2H_pl(g,ADM_kmax+1,l)
          enddo
       enddo
    endif

    call PROF_rapend('____vi_rhow_solver',2)

    return
  end subroutine vi_rhow_solver

end module mod_vi
