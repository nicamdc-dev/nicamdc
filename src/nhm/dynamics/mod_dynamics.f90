!-------------------------------------------------------------------------------
!> Module dynamics
!!
!! @par Description
!!          This module contains the core component of fluid dynamics on icosahedral grid system
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_dynamics
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
  public :: dynamics_setup
  public :: dynamics_step

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
  integer,  private :: num_of_iteration_lstep    ! number of large steps ( 0-4 )
  integer,  private :: num_of_iteration_sstep(4) ! number of small steps in each of large steps

  logical,  private :: trcadv_out_dyndiv
  real(RP), private :: rweight_dyndiv

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine dynamics_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_time, only: &
       TIME_INTEG_TYPE, &
       TIME_SSTEP_MAX
    use mod_runconf, only: &
       TRC_ADV_TYPE,    &
       DYN_DIV_NUM,     &
       TRC_ADV_LOCATION
    use mod_bndcnd, only: &
       BNDCND_setup
    use mod_bsstate, only: &
       bsstate_setup
    use mod_numfilter, only: &
       numfilter_setup
    use mod_vi, only: &
       vi_setup
!TENTATIVE!    use mod_sgs, only: &
!TENTATIVE!       sgs_setup
    use mod_nudge, only: &
       NDG_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[dynamics]/Category[nhm]'

    if( IO_L ) write(IO_FID_LOG,*) '+++ Time integration type: ', trim(TIME_INTEG_TYPE)
    select case(TIME_INTEG_TYPE)
    case('RK2')
       if( IO_L ) write(IO_FID_LOG,*) '+++ 2-stage Runge-Kutta'

       num_of_iteration_lstep    = 2
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX

    case('RK3')
       if( IO_L ) write(IO_FID_LOG,*) '+++ 3-stage Runge-Kutta'

       num_of_iteration_lstep    = 3
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 3
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(3) = TIME_SSTEP_MAX

    case('RK4')
       if( IO_L ) write(IO_FID_LOG,*) '+++ 4-stage Runge-Kutta'

       num_of_iteration_lstep    = 4
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 4
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 3
       num_of_iteration_sstep(3) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(4) = TIME_SSTEP_MAX

    case('TRCADV')
       if( IO_L ) write(IO_FID_LOG,*) '+++ Offline tracer experiment'

       num_of_iteration_lstep    = 0

       if ( TRC_ADV_TYPE == 'DEFAULT' ) then
          write(*,*) 'xxx [dynamics_setup] unsupported advection scheme for TRCADV test! STOP.', &
                     trim(TRC_ADV_TYPE)
          call PRC_MPIstop
       endif

    case default
       write(*,*) 'xxx [dynamics_setup] unsupported integration type! STOP.', trim(TIME_INTEG_TYPE)
       call PRC_MPIstop
    endselect

    trcadv_out_dyndiv = .false.
    if ( TRC_ADV_LOCATION == 'OUT_DYN_DIV_LOOP' ) then
       if ( TRC_ADV_TYPE == 'MIURA2004' ) then
          trcadv_out_dyndiv = .true.
       else
          write(*,*) 'xxx [dynamics_setup] unsupported TRC_ADV_TYPE for OUT_DYN_DIV_LOOP. STOP.', &
                     trim(TRC_ADV_TYPE)
          call PRC_MPIstop
       endif
    endif

    rweight_dyndiv = 1.0_DP / real(DYN_DIV_NUM,kind=DP)

    !---< boundary condition module setup >---
    call BNDCND_setup

    !---< basic state module setup >---
    call bsstate_setup

    !---< numerical filter module setup >---
    call numfilter_setup

    !---< vertical implicit module setup >---
    call vi_setup

    !---< sub-grid scale dynamics module setup >---
!TENTATIVE!    call sgs_setup

    !---< nudging module setup >---
    call NDG_setup

    return
  end subroutine dynamics_setup

  !-----------------------------------------------------------------------------
  !> Main driver of the dynamical step
  subroutine dynamics_step
    use mod_const, only: &
       CONST_Rdry,  &
       CONST_CVdry, &
       CONST_Rvap
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kmax,    &
       ADM_kmin,    &
       ADM_have_sgp   ! [OpenACC]
    use mod_comm, only: &
       COMM_data_transfer, &
       sendbuf_r2r_SP, & ! [OpenACC]
       sendbuf_p2r_SP, & ! [OpenACC]
       sendbuf_r2p_SP, & ! [OpenACC]
       recvbuf_r2r_SP, & ! [OpenACC]
       recvbuf_p2r_SP, & ! [OpenACC]
       recvbuf_r2p_SP, & ! [OpenACC]
       sendbuf_r2r_DP, & ! [OpenACC]
       sendbuf_p2r_DP, & ! [OpenACC]
       sendbuf_r2p_DP, & ! [OpenACC]
       recvbuf_r2r_DP, & ! [OpenACC]
       recvbuf_p2r_DP, & ! [OpenACC]
       recvbuf_r2p_DP, & ! [OpenACC]
       Send_list_r2r,  & ! [OpenACC]
       Send_list_p2r,  & ! [OpenACC]
       Send_list_r2p,  & ! [OpenACC]
       Recv_list_r2r,  & ! [OpenACC]
       Recv_list_p2r,  & ! [OpenACC]
       Recv_list_r2p,  & ! [OpenACC]
       Copy_list_r2r,  & ! [OpenACC]
       Copy_list_p2r,  & ! [OpenACC]
       Copy_list_r2p,  & ! [OpenACC]
       Singular_list     ! [OpenACC]
    use mod_grd, only: &
       GRD_rdgz,  & ! [OpenACC]
       GRD_rdgzh, & ! [OpenACC]
       GRD_afact, & ! [OpenACC]
       GRD_bfact, & ! [OpenACC]
       GRD_cfact, & ! [OpenACC]
       GRD_dfact, & ! [OpenACC]
       GRD_x,     & ! [OpenACC]
       GRD_xr       ! [OpenACC]
    use mod_gmtr, only: &
       GMTR_p, & ! [OpenACC]
       GMTR_t, & ! [OpenACC]
       GMTR_a    ! [OpenACC]
    use mod_oprt, only: &
       OPRT_coef_div,  & ! [OpenACC]
       OPRT_coef_grad, & ! [OpenACC]
       OPRT_coef_lap,  & ! [OpenACC]
       OPRT_coef_intp, & ! [OpenACC]
       OPRT_coef_diff    ! [OpenACC]
    use mod_vmtr, only: &
       VMTR_GSGAM2,       &
       VMTR_GSGAM2_pl,    &
       VMTR_C2Wfact,      &
       VMTR_C2Wfact_pl,   &
       VMTR_C2WfactGz,    &
       VMTR_C2WfactGz_pl, &
       VMTR_PHI,          &
       VMTR_PHI_pl,       &
       VMTR_GAM2H,        & ! [OpenACC]
       VMTR_GSGAM2H,      & ! [OpenACC]
       VMTR_RGSQRTH,      & ! [OpenACC]
       VMTR_RGAM,         & ! [OpenACC]
       VMTR_RGAMH,        & ! [OpenACC]
       VMTR_RGSGAM2,      & ! [OpenACC]
       VMTR_RGSGAM2H,     & ! [OpenACC]
       VMTR_W2Cfact         ! [OpenACC]
    use mod_time, only: &
       TIME_INTEG_TYPE, &
       TIME_DTL,        &
       TIME_DTS,        &
       TIME_SPLIT,      &
       TIME_CTIME
    use mod_runconf, only: &
       I_RHOG,         &
       I_RHOGVX,       &
       I_RHOGVY,       &
       I_RHOGVZ,       &
       I_RHOGW,        &
       I_RHOGE,        &
       I_pre,          &
       I_tem,          &
       I_vx,           &
       I_vy,           &
       I_vz,           &
       I_w,            &
       TRC_VMAX,       &
       I_QV,           &
       I_TKE,          &
       NQW_STR,        &
       NQW_END,        &
       CVW,            &
       DYN_DIV_NUM,    &
       NDIFF_LOCATION, &
       TRC_ADV_TYPE,   &
       FLAG_NUDGING,   &
       THUBURN_LIM
    use mod_prgvar, only: &
       prgvar_get, &
       prgvar_set, &
       PRG_var       ! [OpenACC]
    use mod_bndcnd, only: &
       BNDCND_all
    use mod_bsstate, only: &
       rho_bs,    &
       rho_bs_pl, &
       pre_bs,    &
       pre_bs_pl
    use mod_thrmdyn, only: &
       THRMDYN_th, &
       THRMDYN_eth
    use mod_numfilter, only: &
       NUMFILTER_DOrayleigh,       &
       NUMFILTER_DOverticaldiff,   &
       numfilter_rayleigh_damping, &
       numfilter_hdiffusion,       &
       numfilter_vdiffusion,       &
       Kh_coef,                    & ! [OpenACC]
       Kh_coef_lap1,               & ! [OpenACC]
       divdamp_coef,               & ! [OpenACC]
       divdamp_2d_coef               ! [OpenACC]
    use mod_src, only: &
       src_advection_convergence_momentum, &
       src_advection_convergence,          &
       I_SRC_default
    use mod_vi, only: &
       vi_small_step, &
       Mc, Ml, Mu       ! [OpenACC]
    use mod_src_tracer, only: &
       src_tracer_advection
    use mod_forcing_driver, only: &
       forcing_update
!TENTATIVE!    use mod_sgs, only: &
!TENTATIVE!       sgs_smagorinsky
    use mod_nudge, only: &
       NDG_update_reference, &
       NDG_apply_uvtp
    implicit none

    real(RP) :: PROG         (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! prognostic variables
    real(RP) :: PROG_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: PROGq        (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! tracer variables
    real(RP) :: PROGq_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: g_TEND       (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! tendency of prognostic variables
    real(RP) :: g_TEND_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: g_TENDq      (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! tendency of tracer variables
    real(RP) :: g_TENDq_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: f_TEND       (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! forcing tendency of prognostic variables
    real(RP) :: f_TEND_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: f_TENDq      (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! forcing tendency of tracer variables
    real(RP) :: f_TENDq_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: PROG00       (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! prognostic variables (save)
    real(RP) :: PROG00_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: PROGq00      (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! tracer variables (save)
    real(RP) :: PROGq00_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(RP) :: PROG0        (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! prognostic variables (save)
    real(RP) :: PROG0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    real(RP) :: PROG_split   (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! prognostic variables (split)
    real(RP) :: PROG_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    real(RP) :: PROG_mean    (ADM_gall   ,ADM_kall,ADM_lall   ,5)
    real(RP) :: PROG_mean_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,5)

    ! for tracer advection out of the large step
    real(RP) :: f_TENDrho_mean   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: f_TENDrho_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: f_TENDq_mean     (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP) :: f_TENDq_mean_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(RP) :: PROG_mean_mean   (ADM_gall   ,ADM_kall,ADM_lall   ,5)
    real(RP) :: PROG_mean_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,5)

    real(RP) :: DIAG         (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! diagnostic variables
    real(RP) :: DIAG_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: q            (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! tracer variables
    real(RP) :: q_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !--- density
    real(RP) :: rho   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- internal energy  ( physical )
    real(RP) :: ein   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ein_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- enthalpy ( physical )
    real(RP) :: eth   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: eth_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- potential temperature ( physical )
    real(RP) :: th   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: th_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- density deviation from the base state ( G^1/2 X gamma2 )
    real(RP) :: rhogd   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- pressure deviation from the base state ( G^1/2 X gamma2 )
    real(RP) :: pregd   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: pregd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- temporary variables
    real(RP) :: qd      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qd_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: cv      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: cv_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: TKEg_corr

    integer  :: small_step_ite
    real(RP) :: dyn_step_dt
    real(RP) :: large_step_dt
    real(RP) :: small_step_dt

    logical  :: ndg_TEND_out
    logical  :: do_tke_correction

    integer  :: gall, kall, kmin, kmax, lall, nall, nmin, nmax, iqv, itke
    real(RP) :: Rdry, CVdry, Rvap

    integer  :: g, k ,l, nq, iv, nl, ndyn
    !---------------------------------------------------------------------------
    !$acc wait

    !$acc data &
    !$acc pcreate(g_TEND,g_TENDq,f_TEND,f_TENDq,PROG00,PROGq00,PROG0,PROG_split,PROG_mean) &
    !$acc pcreate(f_TENDrho_mean,f_TENDq_mean,PROG_mean_mean) &
    !$acc pcreate(DIAG,q,rho,ein,eth,th,rhogd,pregd,qd,cv) &
    !$acc pcopyin(ADM_have_sgp) &
    !$acc pcopyin(GRD_rdgz,GRD_rdgzh,GRD_afact,GRD_bfact,GRD_cfact,GRD_dfact,GRD_x,GRD_xr) &
    !$acc pcopy(sendbuf_r2r_SP,sendbuf_p2r_SP,sendbuf_r2p_SP,recvbuf_r2r_SP,recvbuf_p2r_SP,recvbuf_r2p_SP) &
    !$acc pcopy(sendbuf_r2r_DP,sendbuf_p2r_DP,sendbuf_r2p_DP,recvbuf_r2r_DP,recvbuf_p2r_DP,recvbuf_r2p_DP) &
    !$acc pcopyin(Send_list_r2r,Send_list_p2r,Send_list_r2p,Recv_list_r2r,Recv_list_p2r,Recv_list_r2p) &
    !$acc pcopyin(Copy_list_r2r,Copy_list_p2r,Copy_list_r2p,Singular_list) &
    !$acc pcopyin(GMTR_p,GMTR_t,GMTR_a) &
    !$acc pcopyin(OPRT_coef_div,OPRT_coef_grad,OPRT_coef_lap,OPRT_coef_intp,OPRT_coef_diff) &
    !$acc pcopyin(VMTR_GAM2H,VMTR_GSGAM2,VMTR_GSGAM2H,VMTR_RGSQRTH,VMTR_RGAM,VMTR_RGAMH) &
    !$acc pcopyin(VMTR_RGSGAM2,VMTR_RGSGAM2H,VMTR_W2Cfact,VMTR_C2Wfact,VMTR_C2WfactGz,VMTR_PHI) &
    !$acc pcopyin(rho_bs,pre_bs) &
    !$acc pcopyin(Kh_coef,Kh_coef_lap1,divdamp_coef,divdamp_2d_coef,Mc,Ml,Mu) &
    !$acc pcopyin(num_of_iteration_sstep)

    call PROF_rapstart('__Dynamics',1)

    gall = ADM_gall
    kall = ADM_kall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall
    nall = TRC_VMAX
    nmin = NQW_STR
    nmax = NQW_END

    iqv  = I_QV
    itke = I_TKE

    Rdry  = CONST_Rdry
    CVdry = CONST_CVdry
    Rvap  = CONST_Rvap

    call PROF_rapstart('___Pre_Post',1)

    dyn_step_dt   = real(TIME_DTL,kind=RP)
    large_step_dt = real(TIME_DTL,kind=RP) * rweight_dyndiv

#ifdef MORETIMER
call PROF_rapstart('___Pre_Post_01')
#endif

    !--- get from prg0
    call prgvar_get( PROG (:,:,:,I_RHOG),   PROG_pl (:,:,:,I_RHOG),   & ! [OUT]
                     PROG (:,:,:,I_RHOGVX), PROG_pl (:,:,:,I_RHOGVX), & ! [OUT]
                     PROG (:,:,:,I_RHOGVY), PROG_pl (:,:,:,I_RHOGVY), & ! [OUT]
                     PROG (:,:,:,I_RHOGVZ), PROG_pl (:,:,:,I_RHOGVZ), & ! [OUT]
                     PROG (:,:,:,I_RHOGW),  PROG_pl (:,:,:,I_RHOGW),  & ! [OUT]
                     PROG (:,:,:,I_RHOGE),  PROG_pl (:,:,:,I_RHOGE),  & ! [OUT]
                     PROGq(:,:,:,:),        PROGq_pl(:,:,:,:)         ) ! [OUT]

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_01')
#endif

    call PROF_rapend  ('___Pre_Post',1)

    do ndyn = 1, DYN_DIV_NUM

    call PROF_rapstart('___Pre_Post',1)

#ifdef MORETIMER
call PROF_rapstart('___Pre_Post_02')
#endif

    !--- save the value before tracer advection
    if ( ( .NOT. trcadv_out_dyndiv ) .OR. ndyn == 1 ) then

!OCL XFILL
       !$acc kernels pcopy(PROG00) pcopyin(PROG)
       !$omp parallel do default(none),private(g,k,l,iv), &
       !$omp shared(gall,kall,lall,PROG00,PROG), &
       !$omp collapse(3)
       do iv = 1, 6
       do l  = 1, lall
       do k  = 1, kall
       do g  = 1, gall
          PROG00(g,k,l,iv) = PROG(g,k,l,iv)
       enddo
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       PROG00_pl(:,:,:,:) = PROG_pl(:,:,:,:)

       if ( TRC_ADV_TYPE == 'DEFAULT' ) then

!OCL XFILL
          !$acc kernels pcopy(PROGq00) pcopyin(PROGq)
          !$omp parallel do default(none),private(g,k,l,nq), &
          !$omp shared(gall,kall,lall,nall,PROGq00,PROGq), &
          !$omp collapse(3)
          do nq = 1, nall
          do l  = 1, lall
          do k  = 1, kall
          do g  = 1, gall
             PROGq00(g,k,l,nq) = PROGq(g,k,l,nq)
          enddo
          enddo
          enddo
          enddo
          !$omp end parallel do
          !$acc end kernels

          PROGq00_pl(:,:,:,:) = PROGq_pl(:,:,:,:)

       endif
    endif

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_02')
call PROF_rapstart('___Pre_Post_03')
#endif

    !--- save the value before RK loop
!OCL XFILL
    !$acc kernels pcopy(PROG0) pcopyin(PROG)
    !$omp parallel do default(none),private(g,k,l,iv), &
    !$omp shared(gall,kall,lall,PROG0,PROG), &
    !$omp collapse(3)
    do iv = 1, 6
    do l  = 1, lall
    do k  = 1, kall
    do g  = 1, gall
       PROG0(g,k,l,iv) = PROG(g,k,l,iv)
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    PROG0_pl(:,:,:,:) = PROG_pl(:,:,:,:)

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_03')
#endif

    call PROF_rapend  ('___Pre_Post',1)

    if ( TIME_INTEG_TYPE == 'TRCADV' ) then  ! TRC-ADV Test Bifurcation

       call PROF_rapstart('___Tracer_Advection',1)

!OCL XFILL
       !$acc kernels pcopy(f_TEND)
       !$omp parallel do default(none),private(g,k,l,iv), &
       !$omp shared(gall,kall,lall,f_TEND), &
       !$omp collapse(3)
       do iv = 1, 6
       do l  = 1, lall
       do k  = 1, kall
       do g  = 1, gall
          f_TEND(g,k,l,iv) = 0.0_RP
       enddo
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       f_TEND_pl(:,:,:,:) = 0.0_RP

       call src_tracer_advection( TRC_VMAX,                                          & ! [IN]
                                  PROGq (:,:,:,:),        PROGq_pl (:,:,:,:),        & ! [INOUT]
                                  PROG0 (:,:,:,I_RHOG),   PROG0_pl (:,:,:,I_RHOG),   & ! [IN]
                                  PROG  (:,:,:,I_RHOG),   PROG_pl  (:,:,:,I_RHOG),   & ! [IN]
                                  PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & ! [IN]
                                  PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & ! [IN]
                                  PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & ! [IN]
                                  PROG  (:,:,:,I_RHOGW),  PROG_pl  (:,:,:,I_RHOGW),  & ! [IN]
                                  f_TEND(:,:,:,I_RHOG),   f_TEND_pl(:,:,:,I_RHOG),   & ! [IN]
                                  large_step_dt,                                     & ! [IN]
                                  THUBURN_LIM                                        ) ! [IN]

       call PROF_rapend  ('___Tracer_Advection',1)

       call forcing_update( PROG(:,:,:,:), PROG_pl(:,:,:,:) ) ! [INOUT]
    endif

    !---------------------------------------------------------------------------
    !
    !> Start large time step integration
    !
    !---------------------------------------------------------------------------
    do nl = 1, num_of_iteration_lstep

       call PROF_rapstart('___Pre_Post',1)

#ifdef MORETIMER
call PROF_rapstart('___Pre_Post_04')
#endif

       !---< Generate diagnostic values and set the boudary conditions
!OCL XFILL
       !$acc kernels pcopy(DIAG,rho,ein) pcopyin(PROG,VMTR_GSGAM2)
       !$omp parallel do default(none),private(g,k,l), &
       !$omp shared(gall,kall,lall,rho,ein,PROG,DIAG,VMTR_GSGAM2), &
       !$omp collapse(2)
       do l = 1, lall
       do k = 1, kall
       do g = 1, gall
          rho (g,k,l)      = PROG(g,k,l,I_RHOG)   / VMTR_GSGAM2(g,k,l)
          DIAG(g,k,l,I_vx) = PROG(g,k,l,I_RHOGVX) / PROG(g,k,l,I_RHOG)
          DIAG(g,k,l,I_vy) = PROG(g,k,l,I_RHOGVY) / PROG(g,k,l,I_RHOG)
          DIAG(g,k,l,I_vz) = PROG(g,k,l,I_RHOGVZ) / PROG(g,k,l,I_RHOG)
          ein (g,k,l)      = PROG(g,k,l,I_RHOGE)  / PROG(g,k,l,I_RHOG)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
!OCL XFILL
       !$acc kernels pcopy(q) pcopyin(PROGq,PROG)
       !$omp parallel do default(none),private(g,k,l,nq), &
       !$omp shared(gall,kall,lall,nall,q,PROG,PROGq), &
       !$omp collapse(3)
       do nq = 1, nall
       do l  = 1, lall
       do k  = 1, kall
       do g  = 1, gall
          q(g,k,l,nq) = PROGq(g,k,l,nq) / PROG(g,k,l,I_RHOG)
       enddo
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       !$acc kernels pcopy(DIAG,cv,qd) pcopyin(q,ein,rho,CVW)
       !$omp parallel do default(none),private(g,k,l,nq), &
       !$omp shared(gall,kall,lall,nall,nmin,nmax,DIAG,rho,ein,q,qd,cv,CVW,Rdry,Rvap,CVdry,iqv), &
       !$omp collapse(2)
       do l = 1, lall
       do k = 1, kall
       do g = 1, gall
          cv(g,k,l) = 0.0_RP
          qd(g,k,l) = 1.0_RP

          !$acc loop seq
          do nq = nmin, nmax
             cv(g,k,l) = cv(g,k,l) + q(g,k,l,nq) * CVW(nq)
             qd(g,k,l) = qd(g,k,l) - q(g,k,l,nq)
          enddo
          !$acc end loop

          cv(g,k,l) = cv(g,k,l) + qd(g,k,l) * CVdry

          DIAG(g,k,l,I_tem) = ein(g,k,l) / cv(g,k,l)
          DIAG(g,k,l,I_pre) = rho(g,k,l) * DIAG(g,k,l,I_tem) * ( qd(g,k,l)*Rdry + q(g,k,l,iqv)*Rvap )
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
!OCL XFILL
       !$acc kernels pcopy(DIAG) pcopyin(PROG,VMTR_C2Wfact)
       !$omp parallel do default(none),private(g,k,l), &
       !$omp shared(gall,kmin,kmax,lall,DIAG,PROG,VMTR_C2Wfact), &
       !$omp collapse(2)
       do l = 1, lall
       do k = kmin+1, kmax
       do g = 1, gall
          DIAG(g,k,l,I_w) = PROG(g,k,l,I_RHOGW) / ( VMTR_C2Wfact(g,k,1,l) * PROG(g,k  ,l,I_RHOG) &
                                                  + VMTR_C2Wfact(g,k,2,l) * PROG(g,k-1,l,I_RHOG) )
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_04')
call PROF_rapstart('___Pre_Post_05')
#endif

       call BNDCND_all( ADM_gall,                       & ! [IN]
                        ADM_kall,                       & ! [IN]
                        ADM_lall,                       & ! [IN]
                        rho           (:,:,:),          & ! [INOUT]
                        DIAG          (:,:,:,I_vx),     & ! [INOUT]
                        DIAG          (:,:,:,I_vy),     & ! [INOUT]
                        DIAG          (:,:,:,I_vz),     & ! [INOUT]
                        DIAG          (:,:,:,I_w),      & ! [INOUT]
                        ein           (:,:,:),          & ! [INOUT]
                        DIAG          (:,:,:,I_tem),    & ! [INOUT]
                        DIAG          (:,:,:,I_pre),    & ! [INOUT]
                        PROG          (:,:,:,I_RHOG),   & ! [INOUT]
                        PROG          (:,:,:,I_RHOGVX), & ! [INOUT]
                        PROG          (:,:,:,I_RHOGVY), & ! [INOUT]
                        PROG          (:,:,:,I_RHOGVZ), & ! [INOUT]
                        PROG          (:,:,:,I_RHOGW),  & ! [INOUT]
                        PROG          (:,:,:,I_RHOGE),  & ! [INOUT]
                        VMTR_GSGAM2   (:,:,:),          & ! [IN]
                        VMTR_PHI      (:,:,:),          & ! [IN]
                        VMTR_C2Wfact  (:,:,:,:),        & ! [IN]
                        VMTR_C2WfactGz(:,:,:,:)         ) ! [IN]

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_05')
call PROF_rapstart('___Pre_Post_06')
#endif

       call THRMDYN_th ( ADM_gall,          & ! [IN]
                         ADM_kall,          & ! [IN]
                         ADM_lall,          & ! [IN]
                         DIAG(:,:,:,I_tem), & ! [IN]
                         DIAG(:,:,:,I_pre), & ! [IN]
                         th  (:,:,:)        ) ! [OUT]

       call THRMDYN_eth( ADM_gall,          & ! [IN]
                         ADM_kall,          & ! [IN]
                         ADM_lall,          & ! [IN]
                         ein (:,:,:),       & ! [IN]
                         DIAG(:,:,:,I_pre), & ! [IN]
                         rho (:,:,:),       & ! [IN]
                         eth (:,:,:)        ) ! [OUT]

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_06')
call PROF_rapstart('___Pre_Post_07')
#endif

       ! perturbations ( pre, rho with metrics )
!OCL XFILL
       !$acc kernels pcopy(pregd,rhogd) pcopyin(DIAG,pre_bs,rho,rho_bs,VMTR_GSGAM2)
       !$omp parallel do default(none),private(g,k,l), &
       !$omp shared(gall,kall,lall,pregd,rhogd,DIAG,pre_bs,rho,rho_bs,VMTR_GSGAM2), &
       !$omp collapse(2)
       do l = 1, lall
       do k = 1, kall
       do g = 1, gall
          pregd(g,k,l) = ( DIAG(g,k,l,I_pre) - pre_bs(g,k,l) ) * VMTR_GSGAM2(g,k,l)
          rhogd(g,k,l) = ( rho (g,k,l)       - rho_bs(g,k,l) ) * VMTR_GSGAM2(g,k,l)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_07')
call PROF_rapstart('___Pre_Post_08')
#endif

       if ( ADM_have_pl ) then

          rho_pl (:,:,:)      = PROG_pl(:,:,:,I_RHOG)   / VMTR_GSGAM2_pl(:,:,:)
          DIAG_pl(:,:,:,I_vx) = PROG_pl(:,:,:,I_RHOGVX) / PROG_pl(:,:,:,I_RHOG)
          DIAG_pl(:,:,:,I_vy) = PROG_pl(:,:,:,I_RHOGVY) / PROG_pl(:,:,:,I_RHOG)
          DIAG_pl(:,:,:,I_vz) = PROG_pl(:,:,:,I_RHOGVZ) / PROG_pl(:,:,:,I_RHOG)
          ein_pl (:,:,:)      = PROG_pl(:,:,:,I_RHOGE)  / PROG_pl(:,:,:,I_RHOG)

          do nq = 1, TRC_VMAX
             q_pl(:,:,:,nq) = PROGq_pl(:,:,:,nq) / PROG_pl(:,:,:,I_RHOG)
          enddo

          cv_pl(:,:,:) = 0.0_RP
          qd_pl(:,:,:) = 1.0_RP
          do nq = NQW_STR, NQW_END
             cv_pl(:,:,:) = cv_pl(:,:,:) + q_pl(:,:,:,nq) * CVW(nq)
             qd_pl(:,:,:) = qd_pl(:,:,:) - q_pl(:,:,:,nq)
          enddo
          cv_pl(:,:,:) = cv_pl(:,:,:) + qd_pl(:,:,:) * CVdry

          DIAG_pl(:,:,:,I_tem) = ein_pl(:,:,:) / cv_pl(:,:,:)
          DIAG_pl(:,:,:,I_pre) = rho_pl(:,:,:) * DIAG_pl(:,:,:,I_tem) * ( qd_pl(:,:,:)*Rdry + q_pl(:,:,:,I_QV)*Rvap )

          do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             DIAG_pl(g,k,l,I_w) = PROG_pl(g,k,l,I_RHOGW) / ( VMTR_C2Wfact_pl(g,k,1,l) * PROG_pl(g,k  ,l,I_RHOG) &
                                                           + VMTR_C2Wfact_pl(g,k,2,l) * PROG_pl(g,k-1,l,I_RHOG) )
          enddo
          enddo
          enddo

          call BNDCND_all( ADM_gall_pl,                       & ! [IN]
                           ADM_kall,                          & ! [IN]
                           ADM_lall_pl,                       & ! [IN]
                           rho_pl           (:,:,:),          & ! [INOUT]
                           DIAG_pl          (:,:,:,I_vx),     & ! [INOUT]
                           DIAG_pl          (:,:,:,I_vy),     & ! [INOUT]
                           DIAG_pl          (:,:,:,I_vz),     & ! [INOUT]
                           DIAG_pl          (:,:,:,I_w),      & ! [INOUT]
                           ein_pl           (:,:,:),          & ! [INOUT]
                           DIAG_pl          (:,:,:,I_tem),    & ! [INOUT]
                           DIAG_pl          (:,:,:,I_pre),    & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOG),   & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGVX), & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGVY), & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGVZ), & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGW),  & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGE),  & ! [INOUT]
                           VMTR_GSGAM2_pl   (:,:,:),          & ! [IN]
                           VMTR_PHI_pl      (:,:,:),          & ! [IN]
                           VMTR_C2Wfact_pl  (:,:,:,:),        & ! [IN]
                           VMTR_C2WfactGz_pl(:,:,:,:)         ) ! [IN]

          call THRMDYN_th ( ADM_gall_pl,          & ! [IN]
                            ADM_kall,             & ! [IN]
                            ADM_lall_pl,          & ! [IN]
                            DIAG_pl(:,:,:,I_tem), & ! [IN]
                            DIAG_pl(:,:,:,I_pre), & ! [IN]
                            th_pl  (:,:,:)        ) ! [OUT]

          call THRMDYN_eth( ADM_gall_pl,          & ! [IN]
                            ADM_kall,             & ! [IN]
                            ADM_lall_pl,          & ! [IN]
                            ein_pl (:,:,:),       & ! [IN]
                            DIAG_pl(:,:,:,I_pre), & ! [IN]
                            rho_pl (:,:,:),       & ! [IN]
                            eth_pl (:,:,:)        ) ! [OUT]

          pregd_pl(:,:,:) = ( DIAG_pl(:,:,:,I_pre) - pre_bs_pl(:,:,:) ) * VMTR_GSGAM2_pl(:,:,:)
          rhogd_pl(:,:,:) = ( rho_pl (:,:,:)       - rho_bs_pl(:,:,:) ) * VMTR_GSGAM2_pl(:,:,:)
       else

          PROG_pl (:,:,:,:) = 0.0_RP
          DIAG_pl (:,:,:,:) = 0.0_RP

          rho_pl  (:,:,:)   = 0.0_RP
          q_pl    (:,:,:,:) = 0.0_RP

          th_pl   (:,:,:)   = 0.0_RP
          eth_pl  (:,:,:)   = 0.0_RP
          pregd_pl(:,:,:)   = 0.0_RP
          rhogd_pl(:,:,:)   = 0.0_RP

       endif

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_08')
#endif

       !$acc wait
       call PROF_rapend  ('___Pre_Post',1)
       !------------------------------------------------------------------------
       !> LARGE step
       !------------------------------------------------------------------------
       call PROF_rapstart('___Large_step',1)

       !--- calculation of advection tendency including Coriolis force
       call src_advection_convergence_momentum( DIAG  (:,:,:,I_vx),     DIAG_pl  (:,:,:,I_vx),     & ! [IN]
                                                DIAG  (:,:,:,I_vy),     DIAG_pl  (:,:,:,I_vy),     & ! [IN]
                                                DIAG  (:,:,:,I_vz),     DIAG_pl  (:,:,:,I_vz),     & ! [IN]
                                                DIAG  (:,:,:,I_w),      DIAG_pl  (:,:,:,I_w),      & ! [IN]
                                                PROG  (:,:,:,I_RHOG),   PROG_pl  (:,:,:,I_RHOG),   & ! [IN]
                                                PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & ! [IN]
                                                PROG  (:,:,:,I_RHOGW),  PROG_pl  (:,:,:,I_RHOGW),  & ! [IN]
                                                g_TEND(:,:,:,I_RHOGVX), g_TEND_pl(:,:,:,I_RHOGVX), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGVY), g_TEND_pl(:,:,:,I_RHOGVY), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGVZ), g_TEND_pl(:,:,:,I_RHOGVZ), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGW),  g_TEND_pl(:,:,:,I_RHOGW)   ) ! [OUT]

#ifdef MORETIMER
call PROF_rapstart('___Large_step_01')
#endif

!OCL XFILL
       !$acc kernels pcopy(g_TEND)
       !$omp parallel do default(none),private(g,k,l), &
       !$omp shared(gall,kall,lall,g_TEND), &
       !$omp collapse(2)
       do l = 1, lall
       do k = 1, kall
       do g = 1, gall
          g_TEND(g,k,l,I_RHOG)  = 0.0_RP
          g_TEND(g,k,l,I_RHOGE) = 0.0_RP
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       g_TEND_pl(:,:,:,I_RHOG)  = 0.0_RP
       g_TEND_pl(:,:,:,I_RHOGE) = 0.0_RP

#ifdef MORETIMER
call PROF_rapend  ('___Large_step_01')
#endif

       !---< numerical diffusion term
       if ( NDIFF_LOCATION == 'IN_LARGE_STEP' ) then

          if ( nl == 1 ) then ! only first step
             !------ numerical diffusion
             call numfilter_hdiffusion( PROG   (:,:,:,I_RHOG), PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                        rho    (:,:,:),        rho_pl    (:,:,:),        & ! [IN]
                                        DIAG   (:,:,:,I_vx),   DIAG_pl   (:,:,:,I_vx),   & ! [IN]
                                        DIAG   (:,:,:,I_vy),   DIAG_pl   (:,:,:,I_vy),   & ! [IN]
                                        DIAG   (:,:,:,I_vz),   DIAG_pl   (:,:,:,I_vz),   & ! [IN]
                                        DIAG   (:,:,:,I_w),    DIAG_pl   (:,:,:,I_w),    & ! [IN]
                                        DIAG   (:,:,:,I_tem),  DIAG_pl   (:,:,:,I_tem),  & ! [IN]
                                        q      (:,:,:,:),      q_pl      (:,:,:,:),      & ! [IN]
                                        f_TEND (:,:,:,:),      f_TEND_pl (:,:,:,:),      & ! [OUT]
                                        f_TENDq(:,:,:,:),      f_TENDq_pl(:,:,:,:)       ) ! [OUT]

             if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
                call numfilter_vdiffusion( PROG   (:,:,:,I_RHOG), PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                           rho    (:,:,:),        rho_pl    (:,:,:),        & ! [IN]
                                           DIAG   (:,:,:,I_vx),   DIAG_pl   (:,:,:,I_vx),   & ! [IN]
                                           DIAG   (:,:,:,I_vy),   DIAG_pl   (:,:,:,I_vy),   & ! [IN]
                                           DIAG   (:,:,:,I_vz),   DIAG_pl   (:,:,:,I_vz),   & ! [IN]
                                           DIAG   (:,:,:,I_w),    DIAG_pl   (:,:,:,I_w),    & ! [IN]
                                           DIAG   (:,:,:,I_tem),  DIAG_pl   (:,:,:,I_tem),  & ! [IN]
                                           q      (:,:,:,:),      q_pl      (:,:,:,:),      & ! [IN]
                                           f_TEND (:,:,:,:),      f_TEND_pl (:,:,:,:),      & ! [INOUT]
                                           f_TENDq(:,:,:,:),      f_TENDq_pl(:,:,:,:)       ) ! [INOUT]
             endif

             if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
                call numfilter_rayleigh_damping( PROG  (:,:,:,I_RHOG),   PROG_pl  (:,:,:,I_RHOG),   & ! [IN]
                                                 DIAG  (:,:,:,I_vx),     DIAG_pl  (:,:,:,I_vx),     & ! [IN]
                                                 DIAG  (:,:,:,I_vy),     DIAG_pl  (:,:,:,I_vy),     & ! [IN]
                                                 DIAG  (:,:,:,I_vz),     DIAG_pl  (:,:,:,I_vz),     & ! [IN]
                                                 DIAG  (:,:,:,I_w),      DIAG_pl  (:,:,:,I_w),      & ! [IN]
                                                 f_TEND(:,:,:,I_RHOGVX), f_TEND_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGVY), f_TEND_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGVZ), f_TEND_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGW),  f_TEND_pl(:,:,:,I_RHOGW)   ) ! [INOUT]
             endif
          endif

       elseif( NDIFF_LOCATION == 'IN_LARGE_STEP2' ) then

          !------ numerical diffusion
          call numfilter_hdiffusion( PROG   (:,:,:,I_RHOG), PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                     rho    (:,:,:),        rho_pl    (:,:,:),        & ! [IN]
                                     DIAG   (:,:,:,I_vx),   DIAG_pl   (:,:,:,I_vx),   & ! [IN]
                                     DIAG   (:,:,:,I_vy),   DIAG_pl   (:,:,:,I_vy),   & ! [IN]
                                     DIAG   (:,:,:,I_vz),   DIAG_pl   (:,:,:,I_vz),   & ! [IN]
                                     DIAG   (:,:,:,I_w),    DIAG_pl   (:,:,:,I_w),    & ! [IN]
                                     DIAG   (:,:,:,I_tem),  DIAG_pl   (:,:,:,I_tem),  & ! [IN]
                                     q      (:,:,:,:),      q_pl      (:,:,:,:),      & ! [IN]
                                     f_TEND (:,:,:,:),      f_TEND_pl (:,:,:,:),      & ! [OUT]
                                     f_TENDq(:,:,:,:),      f_TENDq_pl(:,:,:,:)       ) ! [OUT]

          if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
             call numfilter_vdiffusion( PROG   (:,:,:,I_RHOG), PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                        rho    (:,:,:),        rho_pl    (:,:,:),        & ! [IN]
                                        DIAG   (:,:,:,I_vx),   DIAG_pl   (:,:,:,I_vx),   & ! [IN]
                                        DIAG   (:,:,:,I_vy),   DIAG_pl   (:,:,:,I_vy),   & ! [IN]
                                        DIAG   (:,:,:,I_vz),   DIAG_pl   (:,:,:,I_vz),   & ! [IN]
                                        DIAG   (:,:,:,I_w),    DIAG_pl   (:,:,:,I_w),    & ! [IN]
                                        DIAG   (:,:,:,I_tem),  DIAG_pl   (:,:,:,I_tem),  & ! [IN]
                                        q      (:,:,:,:),      q_pl      (:,:,:,:),      & ! [IN]
                                        f_TEND (:,:,:,:),      f_TEND_pl (:,:,:,:),      & ! [INOUT]
                                        f_TENDq(:,:,:,:),      f_TENDq_pl(:,:,:,:)       ) ! [INOUT]
          endif

          if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
             call numfilter_rayleigh_damping( PROG  (:,:,:,I_RHOG),   PROG_pl  (:,:,:,I_RHOG),   & ! [IN]
                                              DIAG  (:,:,:,I_vx),     DIAG_pl  (:,:,:,I_vx),     & ! [IN]
                                              DIAG  (:,:,:,I_vy),     DIAG_pl  (:,:,:,I_vy),     & ! [IN]
                                              DIAG  (:,:,:,I_vz),     DIAG_pl  (:,:,:,I_vz),     & ! [IN]
                                              DIAG  (:,:,:,I_w),      DIAG_pl  (:,:,:,I_w),      & ! [IN]
                                              f_TEND(:,:,:,I_RHOGVX), f_TEND_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGVY), f_TEND_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGVZ), f_TEND_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGW),  f_TEND_pl(:,:,:,I_RHOGW)   ) ! [INOUT]
          endif

       endif

!TENTATIVE!       if ( TB_TYPE == 'SMG' ) then ! Smagorinksy-type SGS model
!TENTATIVE!          call sgs_smagorinsky( nl,                                                  &
!TENTATIVE!                                rho    (:,:,:),          rho_pl    (:,:,:),          & ! [IN]
!TENTATIVE!                                PROG   (:,:,:,I_RHOG),   PROG_pl   (:,:,:,I_RHOG),   & ! [IN]
!TENTATIVE!                                PROGq  (:,:,:,:),        PROGq_pl  (:,:,:,:  ),      & ! [IN]
!TENTATIVE!                                DIAG   (:,:,:,I_vx),     DIAG_pl   (:,:,:,I_vx),     & ! [IN]
!TENTATIVE!                                DIAG   (:,:,:,I_vy),     DIAG_pl   (:,:,:,I_vy),     & ! [IN]
!TENTATIVE!                                DIAG   (:,:,:,I_vz),     DIAG_pl   (:,:,:,I_vz),     & ! [IN]
!TENTATIVE!                                DIAG   (:,:,:,I_w),      DIAG_pl   (:,:,:,I_w),      & ! [IN]
!TENTATIVE!                                DIAG   (:,:,:,I_tem),    DIAG_pl   (:,:,:,I_tem),    & ! [IN]
!TENTATIVE!                                q      (:,:,:,:),        q_pl      (:,:,:,:),        & ! [IN]
!TENTATIVE!                                th     (:,:,:),          th_pl     (:,:,:),          & ! [IN]
!TENTATIVE!                                f_TEND (:,:,:,I_RHOG),   f_TEND_pl (:,:,:,I_RHOG),   & ! [INOUT]
!TENTATIVE!                                f_TEND (:,:,:,I_RHOGVX), f_TEND_pl (:,:,:,I_RHOGVX), & ! [INOUT]
!TENTATIVE!                                f_TEND (:,:,:,I_RHOGVY), f_TEND_pl (:,:,:,I_RHOGVY), & ! [INOUT]
!TENTATIVE!                                f_TEND (:,:,:,I_RHOGVZ), f_TEND_pl (:,:,:,I_RHOGVZ), & ! [INOUT]
!TENTATIVE!                                f_TEND (:,:,:,I_RHOGW),  f_TEND_pl (:,:,:,I_RHOGW),  & ! [INOUT]
!TENTATIVE!                                f_TEND (:,:,:,I_RHOGE),  f_TEND_pl (:,:,:,I_RHOGE),  & ! [INOUT]
!TENTATIVE!                                f_TENDq(:,:,:,:),        f_TENDq_pl(:,:,:,:)         ) ! [INOUT]
!TENTATIVE!       endif

       if ( FLAG_NUDGING ) then

          if ( nl == 1 ) then
             call NDG_update_reference( TIME_CTIME )
          endif

          if ( nl == num_of_iteration_lstep ) then
             ndg_TEND_out = .true.
          else
             ndg_TEND_out = .false.
          endif

          call NDG_apply_uvtp( rho   (:,:,:),          rho_pl   (:,:,:),          & ! [IN]
                               DIAG  (:,:,:,I_vx),     DIAG_pl  (:,:,:,I_vx),     & ! [IN]
                               DIAG  (:,:,:,I_vy),     DIAG_pl  (:,:,:,I_vy),     & ! [IN]
                               DIAG  (:,:,:,I_vz),     DIAG_pl  (:,:,:,I_vz),     & ! [IN]
                               DIAG  (:,:,:,I_w),      DIAG_pl  (:,:,:,I_w),      & ! [IN]
                               DIAG  (:,:,:,I_tem),    DIAG_pl  (:,:,:,I_tem),    & ! [IN]
                               DIAG  (:,:,:,I_pre),    DIAG_pl  (:,:,:,I_pre),    & ! [IN]
                               f_TEND(:,:,:,I_RHOGVX), f_TEND_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGVY), f_TEND_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGVZ), f_TEND_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGW),  f_TEND_pl(:,:,:,I_RHOGW),  & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGE),  f_TEND_pl(:,:,:,I_RHOGE),  & ! [INOUT]
                               ndg_TEND_out                                       ) ! [IN] ( TEND out )
       endif

#ifdef MORETIMER
call PROF_rapstart('___Large_step_02')
#endif

       !--- sum the large step TEND ( advection + coriolis + num.diff.,SGS,nudge )
       !$acc kernels pcopy(g_TEND) pcopyin(f_TEND)
       !$omp parallel do default(none),private(g,k,l,iv), &
       !$omp shared(gall,kall,lall,g_TEND,f_TEND), &
       !$omp collapse(3)
       do iv = 1, 6
       do l  = 1, lall
       do k  = 1, kall
       do g  = 1, gall
          g_TEND(g,k,l,iv) = g_TEND(g,k,l,iv) + f_TEND(g,k,l,iv)
       enddo
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       g_TEND_pl(:,:,:,:) = g_TEND_pl(:,:,:,:) + f_TEND_pl(:,:,:,:)

#ifdef MORETIMER
call PROF_rapend  ('___Large_step_02')
#endif

       !$acc wait
       call PROF_rapend  ('___Large_step',1)
       !------------------------------------------------------------------------
       !> SMALL step
       !------------------------------------------------------------------------
       call PROF_rapstart('___Small_step',1)

#ifdef MORETIMER
call PROF_rapstart('___Small_step_01')
#endif

       if ( nl /= 1 ) then ! update split values

!OCL XFILL
          !$acc kernels pcopy(PROG_split) pcopyin(PROG0,PROG)
          !$omp parallel do default(none),private(g,k,l,iv), &
          !$omp shared(gall,kall,lall,PROG_split,PROG0,PROG), &
          !$omp collapse(3)
          do iv = 1, 6
          do l  = 1, lall
          do k  = 1, kall
          do g  = 1, gall
             PROG_split(g,k,l,iv) = PROG0(g,k,l,iv) - PROG(g,k,l,iv)
          enddo
          enddo
          enddo
          enddo
          !$omp end parallel do
          !$acc end kernels

          PROG_split_pl(:,:,:,:) = PROG0_pl(:,:,:,:) - PROG_pl(:,:,:,:)

       else

!OCL XFILL
          !$acc kernels pcopy(PROG_split)
          !$omp parallel do default(none),private(g,k,l,iv), &
          !$omp shared(gall,kall,lall,PROG_split), &
          !$omp collapse(3)
          do iv = 1, 6
          do l  = 1, lall
          do k  = 1, kall
          do g  = 1, gall
             PROG_split(g,k,l,iv) = 0.0_RP
          enddo
          enddo
          enddo
          enddo
          !$omp end parallel do
          !$acc end kernels

          PROG_split_pl(:,:,:,:) = 0.0_RP

       endif

#ifdef MORETIMER
call PROF_rapend  ('___Small_step_01')
#endif

       !------ Core routine for small step
       !------    1. By this subroutine, prognostic variables ( rho,.., rhoge ) are calculated through
       !------       the 'num_of_iteration_sstep(nl)'-th times small step.
       !------    2. grho, grhogvx, ..., and  grhoge has the large step
       !------       tendencies initially, however, they are re-used in this subroutine.
       !------
       if ( TIME_SPLIT ) then
          small_step_ite = num_of_iteration_sstep(nl)
          small_step_dt  = TIME_DTS * rweight_dyndiv
       else
          small_step_ite = 1
          small_step_dt  = large_step_dt / (num_of_iteration_lstep-nl+1)
       endif

       call vi_small_step( PROG      (:,:,:,:),    PROG_pl      (:,:,:,:),    & ! [INOUT] prognostic variables
                           DIAG      (:,:,:,I_vx), DIAG_pl      (:,:,:,I_vx), & ! [IN] diagnostic value
                           DIAG      (:,:,:,I_vy), DIAG_pl      (:,:,:,I_vy), & ! [IN]
                           DIAG      (:,:,:,I_vz), DIAG_pl      (:,:,:,I_vz), & ! [IN]
                           eth       (:,:,:),      eth_pl       (:,:,:),      & ! [IN]
                           rhogd     (:,:,:),      rhogd_pl     (:,:,:),      & ! [IN]
                           pregd     (:,:,:),      pregd_pl     (:,:,:),      & ! [IN]
                           g_TEND    (:,:,:,:),    g_TEND_pl    (:,:,:,:),    & ! [IN] large step TEND
                           PROG_split(:,:,:,:),    PROG_split_pl(:,:,:,:),    & ! [INOUT] split value
                           PROG_mean (:,:,:,:),    PROG_mean_pl(:,:,:,:),     & ! [OUT] mean value
                           small_step_ite,                                    & ! [IN]
                           small_step_dt                                      ) ! [IN]

       !$acc wait
       call PROF_rapend  ('___Small_step',1)
       !------------------------------------------------------------------------
       !>  Tracer advection (in the large step)
       !------------------------------------------------------------------------
       call PROF_rapstart('___Tracer_Advection',1)

       do_tke_correction = .false.

       if ( .NOT. trcadv_out_dyndiv ) then ! calc here or not

       if ( TRC_ADV_TYPE == 'MIURA2004' ) then

          if ( nl == num_of_iteration_lstep ) then

             call src_tracer_advection( TRC_VMAX,                                                & ! [IN]
                                        PROGq    (:,:,:,:),        PROGq_pl    (:,:,:,:),        & ! [INOUT]
                                        PROG00   (:,:,:,I_RHOG),   PROG00_pl   (:,:,:,I_RHOG),   & ! [IN]
                                        PROG_mean(:,:,:,I_RHOG),   PROG_mean_pl(:,:,:,I_RHOG),   & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVX), PROG_mean_pl(:,:,:,I_RHOGVX), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVY), PROG_mean_pl(:,:,:,I_RHOGVY), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVZ), PROG_mean_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGW),  PROG_mean_pl(:,:,:,I_RHOGW),  & ! [IN]
                                        f_TEND   (:,:,:,I_RHOG),   f_TEND_pl   (:,:,:,I_RHOG),   & ! [IN]
                                        large_step_dt,                                           & ! [IN]
                                        THUBURN_LIM                                              ) ! [IN]

#ifdef MORETIMER
call PROF_rapstart('___Tracer_Advection_01')
#endif

             !$acc kernels pcopy(PROGq) pcopyin(f_TENDq)
             !$omp parallel do default(none),private(g,k,l,nq), &
             !$omp shared(gall,kall,lall,nall,PROGq,f_TENDq,large_step_dt), &
             !$omp collapse(3)
             do nq = 1, nall
             do l  = 1, lall
             do k  = 1, kall
             do g  = 1, gall
                PROGq(g,k,l,nq) = PROGq(g,k,l,nq) + large_step_dt * f_TENDq(g,k,l,nq) ! update rhogq by viscosity
             enddo
             enddo
             enddo
             enddo
             !$omp end parallel do
             !$acc end kernels

             if ( ADM_have_pl ) then
                PROGq_pl(:,:,:,:) = PROGq_pl(:,:,:,:) + large_step_dt * f_TENDq_pl(:,:,:,:)
             endif

#ifdef MORETIMER
call PROF_rapend  ('___Tracer_Advection_01')
#endif

             ! [comment] H.Tomita: I don't recommend adding the hyperviscosity term because of numerical instability in this case.
             if( I_TKE >= 0 ) do_tke_correction = .true.

          endif ! Last large step only

       elseif( TRC_ADV_TYPE == 'DEFAULT' ) then

          do nq = 1, TRC_VMAX

             call src_advection_convergence( PROG_mean(:,:,:,I_RHOGVX), PROG_mean_pl(:,:,:,I_RHOGVX), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGVY), PROG_mean_pl(:,:,:,I_RHOGVY), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGVZ), PROG_mean_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGW),  PROG_mean_pl(:,:,:,I_RHOGW),  & ! [IN]
                                             q        (:,:,:,nq),       q_pl        (:,:,:,nq),       & ! [IN]
                                             g_TENDq  (:,:,:,nq),       g_TENDq_pl  (:,:,:,nq),       & ! [OUT]
                                             I_SRC_default                                            ) ! [IN]

          enddo ! tracer LOOP

!OCL XFILL
          !$acc kernels pcopy(PROGq) pcopyin(PROGq00,g_TENDq,f_TENDq,num_of_iteration_sstep)
          !$omp parallel default(none),private(g,k,l,nq), &
          !$omp shared(gall,kall,lall,nall,nl,kmin,kmax,PROGq,PROGq00,g_TENDq,f_TENDq,num_of_iteration_sstep,small_step_dt)
          do nq = 1, nall
          do l  = 1, lall
             !$omp do
             do k = 1, kall
             do g = 1, gall
                PROGq(g,k,l,nq) = PROGq00(g,k,l,nq)                                                                        &
                                + ( num_of_iteration_sstep(nl) * small_step_dt ) * ( g_TENDq(g,k,l,nq) + f_TENDq(g,k,l,nq) )
             enddo
             enddo
             !$omp end do

             !$omp do
             do g = 1, gall
                PROGq(g,kmin-1,l,nq) = 0.0_RP
                PROGq(g,kmax+1,l,nq) = 0.0_RP
             enddo
             !$omp end do
          enddo
          enddo
          !$omp end parallel
          !$acc end kernels

          if ( ADM_have_pl ) then
             PROGq_pl(:,:,:,:) = PROGq00_pl(:,:,:,:)                                                                           &
                               + ( num_of_iteration_sstep(nl) * small_step_dt ) * ( g_TENDq_pl(:,:,:,:) + f_TENDq_pl(:,:,:,:) )

             PROGq_pl(:,ADM_kmin-1,:,:) = 0.0_RP
             PROGq_pl(:,ADM_kmax+1,:,:) = 0.0_RP
          endif

          if( I_TKE >= 0 ) do_tke_correction = .true.

       endif

       ! TKE fixer
       if ( do_tke_correction ) then

#ifdef MORETIMER
call PROF_rapstart('___Tracer_Advection_02')
#endif

          !$acc kernels pcopy(PROG,PROGq) pcopyin(VMTR_GSGAM2)
          do l = 1, ADM_lall
             !$omp parallel do default(none),private(g,k,TKEG_corr), &
             !$omp shared(l,gall,kall,PROG,PROGq,itke)
             do k = 1, kall
             do g = 1, gall
                TKEG_corr = max( -PROGq(g,k,l,itke), 0.0_RP )

                PROG (g,k,l,I_RHOGE) = PROG (g,k,l,I_RHOGE) - TKEG_corr
                PROGq(g,k,l,itke)    = PROGq(g,k,l,itke)    + TKEG_corr
             enddo
             enddo
             !$omp end parallel do
          enddo
          !$acc end kernels

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                TKEg_corr = max( -PROGq_pl(g,k,l,I_TKE), 0.0_RP )

                PROG_pl (g,k,l,I_RHOGE) = PROG_pl (g,k,l,I_RHOGE) - TKEG_corr
                PROGq_pl(g,k,l,I_TKE)   = PROGq_pl(g,k,l,I_TKE)   + TKEG_corr
             enddo
             enddo
             enddo
          endif

#ifdef MORETIMER
call PROF_rapend  ('___Tracer_Advection_02')
#endif

       endif

       else

          !--- calculation of mean ( mean mass flux and tendency )
          if ( nl == num_of_iteration_lstep ) then
             if ( ndyn == 1 ) then
!OCL XFILL
                !$acc kernels pcopy(PROG_mean_mean) pcopyin(PROG_mean)
                !$omp parallel do default(none),private(g,k,l,iv), &
                !$omp shared(gall,kall,lall,PROG_mean_mean,rweight_dyndiv,PROG_mean), &
                !$omp collapse(3)
                do iv = 1, 5
                do l  = 1, lall
                do k  = 1, kall
                do g  = 1, gall
                   PROG_mean_mean(g,k,l,iv) = rweight_dyndiv * PROG_mean(g,k,l,iv)
                enddo
                enddo
                enddo
                enddo
                !$omp end parallel do
                !$acc end kernels
!OCL XFILL
                !$acc kernels pcopy(f_TENDrho_mean) pcopyin(f_TEND)
                !$omp parallel do default(none),private(g,k,l), &
                !$omp shared(gall,kall,lall,f_TENDrho_mean,rweight_dyndiv,f_TEND), &
                !$omp collapse(2)
                do l  = 1, lall
                do k  = 1, kall
                do g  = 1, gall
                   f_TENDrho_mean(g,k,l) = rweight_dyndiv * f_TEND(g,k,l,I_RHOG)
                enddo
                enddo
                enddo
                !$omp end parallel do
                !$acc end kernels
!OCL XFILL
                !$acc kernels pcopy(f_TENDq_mean) pcopyin(f_TENDq)
                !$omp parallel do default(none),private(g,k,l,nq), &
                !$omp shared(gall,kall,lall,nall,f_TENDq_mean,rweight_dyndiv,f_TENDq), &
                !$omp collapse(3)
                do nq = 1, nall
                do l  = 1, lall
                do k  = 1, kall
                do g  = 1, gall
                   f_TENDq_mean(g,k,l,nq) = rweight_dyndiv * f_TENDq(g,k,l,nq)
                enddo
                enddo
                enddo
                enddo
                !$omp end parallel do
                !$acc end kernels

                PROG_mean_mean_pl(:,:,:,:) = rweight_dyndiv * PROG_mean_pl(:,:,:,:)
                f_TENDrho_mean_pl(:,:,:)   = rweight_dyndiv * f_TEND_pl   (:,:,:,I_RHOG)
                f_TENDq_mean_pl  (:,:,:,:) = rweight_dyndiv * f_TENDq_pl  (:,:,:,:)
             else
                !$acc kernels pcopy(PROG_mean_mean) pcopyin(PROG_mean)
                !$omp parallel do default(none),private(g,k,l,iv), &
                !$omp shared(gall,kall,lall,PROG_mean_mean,rweight_dyndiv,PROG_mean), &
                !$omp collapse(3)
                do iv = 1, 5
                do l  = 1, lall
                do k  = 1, kall
                do g  = 1, gall
                   PROG_mean_mean(g,k,l,iv) = PROG_mean_mean(g,k,l,iv) + rweight_dyndiv * PROG_mean(g,k,l,iv)
                enddo
                enddo
                enddo
                enddo
                !$omp end parallel do
                !$acc end kernels

                !$acc kernels pcopy(f_TENDrho_mean) pcopyin(f_TEND)
                !$omp parallel do default(none),private(g,k,l), &
                !$omp shared(gall,kall,lall,f_TENDrho_mean,rweight_dyndiv,f_TEND), &
                !$omp collapse(2)
                do l  = 1, lall
                do k  = 1, kall
                do g  = 1, gall
                   f_TENDrho_mean(g,k,l) = f_TENDrho_mean(g,k,l) + rweight_dyndiv * f_TEND(g,k,l,I_RHOG)
                enddo
                enddo
                enddo
                !$omp end parallel do
                !$acc end kernels

                !$acc kernels pcopy(f_TENDq_mean) pcopyin(f_TENDq)
                !$omp parallel do default(none),private(g,k,l,nq), &
                !$omp shared(gall,kall,lall,nall,f_TENDq_mean,rweight_dyndiv,f_TENDq), &
                !$omp collapse(3)
                do nq = 1, nall
                do l  = 1, lall
                do k  = 1, kall
                do g  = 1, gall
                   f_TENDq_mean(g,k,l,nq) = f_TENDq_mean(g,k,l,nq) + rweight_dyndiv * f_TENDq(g,k,l,nq)
                enddo
                enddo
                enddo
                enddo
                !$omp end parallel do
                !$acc end kernels

                PROG_mean_mean_pl(:,:,:,:) = PROG_mean_mean_pl(:,:,:,:) + rweight_dyndiv * PROG_mean_pl(:,:,:,:)
                f_TENDrho_mean_pl(:,:,:)   = f_TENDrho_mean_pl(:,:,:)   + rweight_dyndiv * f_TEND_pl   (:,:,:,I_RHOG)
                f_TENDq_mean_pl  (:,:,:,:) = f_TENDq_mean_pl  (:,:,:,:) + rweight_dyndiv * f_TENDq_pl  (:,:,:,:)
             endif
          endif

       endif

       call PROF_rapend  ('___Tracer_Advection',1)
       call PROF_rapstart('___Pre_Post',1)

#ifdef MORETIMER
call PROF_rapstart('___Pre_Post_09c')
#endif

       !------ Update
       if ( nl /= num_of_iteration_lstep ) then
          call COMM_data_transfer( PROG, PROG_pl )
       endif

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_09c')
#endif

       call PROF_rapend  ('___Pre_Post',1)

    enddo !--- large step

    !---------------------------------------------------------------------------
    !>  Tracer advection (out of the large step)
    !---------------------------------------------------------------------------
    if ( trcadv_out_dyndiv .AND. ndyn == DYN_DIV_NUM ) then
       call PROF_rapstart('___Tracer_Advection',1)

       call src_tracer_advection( TRC_VMAX,                                                          & ! [IN]
                                  PROGq         (:,:,:,:),        PROGq_pl         (:,:,:,:),        & ! [INOUT]
                                  PROG00        (:,:,:,I_RHOG),   PROG00_pl        (:,:,:,I_RHOG),   & ! [IN]
                                  PROG_mean_mean(:,:,:,I_RHOG),   PROG_mean_mean_pl(:,:,:,I_RHOG),   & ! [IN]
                                  PROG_mean_mean(:,:,:,I_RHOGVX), PROG_mean_mean_pl(:,:,:,I_RHOGVX), & ! [IN]
                                  PROG_mean_mean(:,:,:,I_RHOGVY), PROG_mean_mean_pl(:,:,:,I_RHOGVY), & ! [IN]
                                  PROG_mean_mean(:,:,:,I_RHOGVZ), PROG_mean_mean_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                  PROG_mean_mean(:,:,:,I_RHOGW),  PROG_mean_mean_pl(:,:,:,I_RHOGW),  & ! [IN]
                                  f_TENDrho_mean(:,:,:),          f_TENDrho_mean_pl(:,:,:),          & ! [IN]
                                  dyn_step_dt,                                                       & ! [IN]
                                  THUBURN_LIM                                                        ) ! [IN]

#ifdef MORETIMER
call PROF_rapstart('___Tracer_Advection_03')
#endif

             !$acc kernels pcopy(PROGq) pcopyin(f_TENDq)
             !$omp parallel do default(none),private(g,k,l,nq), &
             !$omp shared(gall,kall,lall,nall,PROGq,f_TENDq_mean,dyn_step_dt), &
             !$omp collapse(3)
             do nq = 1, nall
             do l  = 1, lall
             do k  = 1, kall
             do g  = 1, gall
                PROGq(g,k,l,nq) = PROGq(g,k,l,nq) + dyn_step_dt * f_TENDq_mean(g,k,l,nq) ! update rhogq by viscosity
             enddo
             enddo
             enddo
             enddo
             !$omp end parallel do
             !$acc end kernels

             if ( ADM_have_pl ) then
                PROGq_pl(:,:,:,:) = PROGq_pl(:,:,:,:) + dyn_step_dt * f_TENDq_mean_pl(:,:,:,:)
             endif

#ifdef MORETIMER
call PROF_rapend  ('___Tracer_Advection_03')
call PROF_rapstart('___Tracer_Advection_04')
#endif

          !$acc kernels pcopy(PROG,PROGq) pcopyin(VMTR_GSGAM2)
          do l = 1, ADM_lall
             !$omp parallel do default(none),private(g,k,TKEG_corr), &
             !$omp shared(l,gall,kall,PROG,PROGq,itke)
             do k = 1, kall
             do g = 1, gall
                TKEG_corr = max( -PROGq(g,k,l,itke), 0.0_RP )

                PROG (g,k,l,I_RHOGE) = PROG (g,k,l,I_RHOGE) - TKEG_corr
                PROGq(g,k,l,itke)    = PROGq(g,k,l,itke)    + TKEG_corr
             enddo
             enddo
             !$omp end parallel do
          enddo
          !$acc end kernels

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                TKEg_corr = max( -PROGq_pl(g,k,l,I_TKE), 0.0_RP )

                PROG_pl (g,k,l,I_RHOGE) = PROG_pl (g,k,l,I_RHOGE) - TKEG_corr
                PROGq_pl(g,k,l,I_TKE)   = PROGq_pl(g,k,l,I_TKE)   + TKEG_corr
             enddo
             enddo
             enddo
          endif

#ifdef MORETIMER
call PROF_rapend  ('___Tracer_Advection_04')
#endif

       call PROF_rapend  ('___Tracer_Advection',1)
    endif

    enddo !--- divided step for dynamics

    call PROF_rapstart('___Pre_Post',1)

#ifdef MORETIMER
call PROF_rapstart('___Pre_Post_10')
#endif

    call prgvar_set( PROG(:,:,:,I_RHOG),   PROG_pl(:,:,:,I_RHOG),   & ! [IN]
                     PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & ! [IN]
                     PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & ! [IN]
                     PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ), & ! [IN]
                     PROG(:,:,:,I_RHOGW),  PROG_pl(:,:,:,I_RHOGW),  & ! [IN]
                     PROG(:,:,:,I_RHOGE),  PROG_pl(:,:,:,I_RHOGE),  & ! [IN]
                     PROGq(:,:,:,:),       PROGq_pl(:,:,:,:)        ) ! [IN]

#ifdef MORETIMER
call PROF_rapend  ('___Pre_Post_10')
#endif

    call PROF_rapend  ('___Pre_Post',1)

!     call PROF_rapstart('___Pre_Post',1)
!
!     !=> Niwa [TM]
!     if (trim(TM_RUN) == 'ON' .AND. TM_OUT_METVAR) then
!        call tm_metvar_setmfx( &
!             v_mean_c(:,:,:,I_rhogvx),  v_mean_c_pl(:,:,:,I_rhogvx), &
!             v_mean_c(:,:,:,I_rhogvy),  v_mean_c_pl(:,:,:,I_rhogvy), &
!             v_mean_c(:,:,:,I_rhogvz),  v_mean_c_pl(:,:,:,I_rhogvz), &
!             v_mean_c(:,:,:,I_rhogw),   v_mean_c_pl(:,:,:,I_rhogw)   )
!        !
!        call tm_metvar_setmet( &
!             rhog,                   rhog_pl, &
!             rhoge,                  rhoge_pl, &
!             rhogq(:,:,:,I_QV:I_QC), rhogq_pl(:,:,:,I_QV:I_QC) &
!             )
!     endif
!     !<= Niwa [TM]
!
!     call PROF_rapend  ('___Pre_Post',1)

    call PROF_rapend  ('__Dynamics',1)

    !$acc end data
    !$acc wait

    return
  end subroutine dynamics_step

end module mod_dynamics
