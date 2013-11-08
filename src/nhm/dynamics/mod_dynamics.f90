!-------------------------------------------------------------------------------
!>
!! Dynamical core
!!
!! @par Description
!!         This module is the core module of fluid dynamics
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita ) Imported from igdc-4.33
!! @li      2006-04-17 (H.Tomita ) Add IN_LARGE_STEP2
!! @li      2007-05-08 (H.Tomita ) Change the treatment of I_TKE.
!! @li      2008-01-24 (Y.Niwa   ) add revised MIURA2004 for tracer advection
!! @li      2008-05-24 (T.Mitsui ) fix miss-conditioning for frcvar
!! @li      2008-09-09 (Y.Niwa   ) move nudging routine here
!! @li      2009-09-08 (S.Iga    ) frhog and frhog_pl in ndg are deleted ( suggested by ES staff)
!! @li      2010-05-06 (M.Satoh  ) define QV_conv only if CP_TYPE='TDK' .or. 'KUO'
!! @li      2010-07-16 (A.T.Noda ) bug fix for TDK
!! @li      2010-08-16 (A.T.Noda ) Bug fix (Qconv not diveded by density)
!! @li      2010-08-20 (A.T.Noda ) Bug fix (Qconv should be TEND, and not be multiplied by DT)
!! @li      2010-11-29 (A.T.Noda ) Introduce the Smagorinsky model
!! @li      2011-08-16 (M.Satoh  ) bug fix for TDK: conv => TEND; qv_dyn_tend = v grad q = ( div(rho v q) - div(rho v)*q )/rho
!! @li      2011-08-16 (M.Satoh  ) move codes related to CP_TYPE below the tracer calculation
!! @li      2013-08-23 (H.Yashiro) Change arguments from character to index/switch
!! @li      2013-06-13 (R.Yoshida) add tracer advection mode
!!
!<
module mod_dynamics
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
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
  integer, private, parameter :: I_RHOG     = 1 ! Density x G^{1/2} x gamma^2
  integer, private, parameter :: I_RHOGVX   = 2 ! Density x G^{1/2} x gamma^2 x Horizontal velocity (X-direction)
  integer, private, parameter :: I_RHOGVY   = 3 ! Density x G^{1/2} x gamma^2 x Horizontal velocity (Y-direction)
  integer, private, parameter :: I_RHOGVZ   = 4 ! Density x G^{1/2} x gamma^2 x Horizontal velocity (Z-direction)
  integer, private, parameter :: I_RHOGW    = 5 ! Density x G^{1/2} x gamma^2 x Vertical   velocity
  integer, private, parameter :: I_RHOGE    = 6 ! Density x G^{1/2} x gamma^2 x Internal Energy
  integer, private, parameter :: I_RHOGETOT = 7 ! Density x G^{1/2} x gamma^2 x Total Energy

  integer, private, save :: num_of_iteration_lstep    ! number of large steps ( 0-4 )
  integer, private, save :: num_of_iteration_sstep(4) ! number of small steps in each of large steps

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine dynamics_setup
    use mod_adm, only: &
       ADM_proc_stop
    use mod_time, only: &
       TIME_INTEG_TYPE, &
       TIME_SSTEP_MAX, &
     ctime => TIME_CTIME, &
     dtime => TIME_DTL
    use mod_runconf, only: &
       TRC_ADV_TYPE, &
       FLAG_NUDGING
  use mod_bsstate, only: &
     bsstate_setup
  use mod_bndcnd, only   :   &
     bndcnd_setup
  use mod_numfilter, only  : &
     numfilter_setup
  use mod_nudge, only: &
     NDG_setup
    use mod_sgs, only: &
       sgs_setup
    implicit none
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[dynamics]/Category[nhm]'

    write(ADM_LOG_FID,*) '+++ Time integration type: ', trim(TIME_INTEG_TYPE)
    select case(TIME_INTEG_TYPE)
    case('RK2')
       write(ADM_LOG_FID,*) '+++ 2-stage Runge-Kutta'

       num_of_iteration_lstep    = 2
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX

    case('RK3')
       write(ADM_LOG_FID,*) '+++ 3-stage Runge-Kutta'

       num_of_iteration_lstep    = 3
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 3
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(3) = TIME_SSTEP_MAX

    case('RK4')
       write(ADM_LOG_FID,*) '+++ 4-stage Runge-Kutta'

       num_of_iteration_lstep    = 4
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 4
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 3
       num_of_iteration_sstep(3) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(4) = TIME_SSTEP_MAX

    case('TRCADV') ! R.Yoshida 13/06/13 [add]
       write(ADM_LOG_FID,*) '+++ Offline tracer experiment'

       num_of_iteration_lstep    = 0

       if ( TRC_ADV_TYPE == 'DEFAULT' ) then
          write(ADM_LOG_FID,*) 'xxx unsupported advection scheme for TRCADV test! STOP.'
          call ADM_proc_stop
       endif

    case default
       write(ADM_LOG_FID,*) 'xxx unsupported integration type! STOP.'
       call ADM_proc_stop
    endselect


  !---< boundary condition module setup >---
  call bndcnd_setup

  !---< basic state module setup >---
  call bsstate_setup

  !---< numerical filter module setup >---
  call numfilter_setup

  !---< nudging module setup >---
  if( FLAG_NUDGING ) call NDG_setup


    call sgs_setup

    return
  end subroutine dynamics_setup

  !-----------------------------------------------------------------------------
  subroutine dynamics_step
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gmax,    &
       ADM_gmin,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_cnst, only: &
       CNST_RAIR, &
       CNST_RVAP, &
       CNST_CV
    use mod_time, only: &
       TIME_INTEG_TYPE, &
       TIME_DTL,        &
       TIME_DTS,        &
       TIME_SPLIT,      &
       TIME_CTIME
    use mod_vmtr, only: &
       VMTR_GSGAM2,        &
       VMTR_GSGAM2_pl,     &
       VMTR_GSGAM2H,       &
       VMTR_GSGAM2H_pl,    &
       VMTR_PHI,           &
       VMTR_PHI_pl,        &
       VMTR_C2Wfact_Gz,    &
       VMTR_C2Wfact_Gz_pl, &
       VMTR_C2Wfact,       &
       VMTR_C2Wfact_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_runconf, only: &
       TRC_VMAX,       &
       I_QV,           &
       I_TKE,          &
       NQW_STR,        &
       NQW_END,        &
       CVW,            &
       NDIFF_LOCATION, &
       TRC_ADV_TYPE,   &
       FLAG_NUDGING,   &
       TB_TYPE,        &
       THUBURN_LIM       ! R.Yoshida 13/06/13 [add]
    use mod_thrmdyn, only: &
       thrmdyn_th, &
       thrmdyn_eth
    use mod_bndcnd, only: &
       bndcnd_all
    use mod_bsstate, only: &
       pre_bs, pre_bs_pl, &
       tem_bs, tem_bs_pl, &
       rho_bs, rho_bs_pl
    use mod_prgvar, only: &
       prgvar_set, &
       prgvar_get
    use mod_src, only: &
       src_advection_convergence_momentum, &
       src_advection_convergence,          &
       I_SRC_default
    use mod_numfilter, only: &
       NUMFILTER_DOrayleigh,       & ! [add] H.Yashiro 20120530
       NUMFILTER_DOverticaldiff,   & ! [add] H.Yashiro 20120530
       numfilter_rayleigh_damping, &
       numfilter_hdiffusion,       &
       numfilter_vdiffusion
    use mod_nudge, only: & ! Y.Niwa add 08/09/09
       NDG_update_reference, &
       NDG_apply_uvtp
    use mod_sgs, only: &
       sgs_smagorinsky
    use mod_vi, only: &
       vi_small_step
    use mod_src_tracer, only: &
       src_tracer_advection
    use mod_forcing_driver, only: &
       forcing_update ! R.Yoshida 13/06/13 [add]
    implicit none

    integer, parameter :: nmax_TEND     = 7
    integer, parameter :: nmax_PROG     = 6
    integer, parameter :: nmax_v_mean_c = 5

    real(8) :: g_TEND    (ADM_gall,   ADM_kall,ADM_lall,   nmax_TEND) !--- tendency
    real(8) :: g_TEND_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,nmax_TEND)
    real(8) :: g_TENDq   (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)  !--- tendency of q
    real(8) :: g_TENDq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(8) :: f_TEND    (ADM_gall,   ADM_kall,ADM_lall,   nmax_TEND) !--- forcing tendency
    real(8) :: f_TEND_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,nmax_TEND)
    real(8) :: f_TENDq   (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)  !--- forcing tendency of q
    real(8) :: f_TENDq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(8) :: PROG0     (ADM_gall,   ADM_kall,ADM_lall,   nmax_PROG) !--- prognostic variables (save)
    real(8) :: PROG0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl,nmax_PROG)
    real(8) :: PROGq0    (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)  !--- tracer variables (save)
    real(8) :: PROGq0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(8) :: PROG      (ADM_gall,   ADM_kall,ADM_lall,   nmax_PROG) !--- prognostic variables
    real(8) :: PROG_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,nmax_PROG)
    real(8) :: PROGq     (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)  !--- tracer variables
    real(8) :: PROGq_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(8) :: PROG_split   (ADM_gall,   ADM_kall,ADM_lall,   nmax_PROG) !--- prognostic variables (split)
    real(8) :: PROG_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,nmax_PROG)

    real(8) :: v_mean_c   (ADM_gall,   ADM_kall,ADM_lall   ,nmax_v_mean_c)
    real(8) :: v_mean_c_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,nmax_v_mean_c)

    !--- density ( physical )
    real(8) :: rho   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rho_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- horizontal velocity_x  ( physical )
    real(8) :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- horizontal velocity_y  ( physical )
    real(8) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- horizontal velocity_z  ( physical )
    real(8) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- vertical velocity ( physical )
    real(8) :: w   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: w_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- internal energy  ( physical )
    real(8) :: ein   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: ein_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- mass concentration of water substance ( physical )
    real(8) :: q   (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)
    real(8) :: q_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !--- enthalpy ( physical )
    real(8) :: eth   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: eth_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- pressure ( physical )
    real(8) :: pre   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: pre_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- temperature ( physical )
    real(8) :: tem   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: tem_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- potential temperature ( physical )
    real(8) :: th   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: th_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- density deviation from the base state ( G^{1/2} X gamma2 )
    real(8) :: rhogd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- pressure deviation from the base state ( G^{1/2} X gamma2 )
    real(8) :: pregd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: pregd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- temperature deviation from the base state ( physical )
    real(8) :: temd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: temd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- temporary variables
    real(8) :: qd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: qd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: cv   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: cv_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), parameter :: TKE_MIN = 0.D0
    real(8)            :: TKEg_corr

    integer :: small_step_ite
    real(8) :: small_step_dt

    logical :: ndg_TEND_out
    logical :: do_tke_correction

    integer :: g, k ,l, nq, nl

    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++Dynamics')

    !--- get from prg0
    call prgvar_get( PROG(:,:,:,I_RHOG),   PROG_pl(:,:,:,I_RHOG),   & ! [OUT]
                     PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & ! [OUT]
                     PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & ! [OUT]
                     PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ), & ! [OUT]
                     PROG(:,:,:,I_RHOGW),  PROG_pl(:,:,:,I_RHOGW),  & ! [OUT]
                     PROG(:,:,:,I_RHOGE),  PROG_pl(:,:,:,I_RHOGE),  & ! [OUT]
                     PROGq(:,:,:,:),       PROGq_pl(:,:,:,:),       & ! [OUT]
                     0                                              ) ! [IN]

    !--- save
    PROG0   (:,:,:,:) = PROG   (:,:,:,:)
    PROG0_pl(:,:,:,:) = PROG_pl(:,:,:,:)

    if ( TRC_ADV_TYPE == 'DEFAULT' ) then
       PROGq0   (:,:,:,:) = PROGq   (:,:,:,:)
       PROGq0_pl(:,:,:,:) = PROGq_pl(:,:,:,:)
    endif

    if ( TIME_INTEG_TYPE == 'TRCADV' ) then  ! TRC-ADV Test Bifurcation

       call DEBUG_rapstart('+++Tracer Advection')

       f_TEND   (:,:,:,:) = 0.D0
       f_TEND_pl(:,:,:,:) = 0.D0

       call src_tracer_advection( TRC_VMAX,                                        & ! [IN]
                                  PROGq(:,:,:,:),        PROGq_pl(:,:,:,:),        & ! [INOUT]
                                  PROG0(:,:,:,I_RHOG),   PROG0_pl(:,:,:,I_RHOG),   & ! [IN]
                                  PROG0(:,:,:,I_rhog),   PROG0_pl(:,:,:,I_rhog),   & ! [IN]
                                  PROG0(:,:,:,I_rhogvx), PROG0_pl(:,:,:,I_rhogvx), & ! [IN]
                                  PROG0(:,:,:,I_rhogvy), PROG0_pl(:,:,:,I_rhogvy), & ! [IN]
                                  PROG0(:,:,:,I_rhogvz), PROG0_pl(:,:,:,I_rhogvz), & ! [IN]
                                  PROG0(:,:,:,I_rhogw),  PROG0_pl(:,:,:,I_rhogw),  & ! [IN]
                                  f_TEND(:,:,:,I_RHOG),  f_TEND_pl(:,:,:,I_RHOG),  & ! [IN]
                                  TIME_DTL,                                        & ! [IN]
                                  THUBURN_LIM                                      ) ! [IN] [add] 20130613 R.Yoshida

       call DEBUG_rapend  ('+++Tracer Advection')

       call forcing_update( PROG0(:,:,:,:), PROG0_pl(:,:,:,:), & ! [IN]
                            PROG (:,:,:,:), PROG_pl (:,:,:,:)  ) ! [INOUT]
    endif

    !---------------------------------------------------------------------------
    !
    !> Start large time step integration
    !
    !---------------------------------------------------------------------------
    do nl = 1, num_of_iteration_lstep

       !---< Generate diagnostic values and set the boudary conditions
       rho(:,:,:) = PROG(:,:,:,I_RHOG  ) / VMTR_GSGAM2(:,:,:)
       vx (:,:,:) = PROG(:,:,:,I_RHOGVX) / PROG(:,:,:,I_RHOG)
       vy (:,:,:) = PROG(:,:,:,I_RHOGVY) / PROG(:,:,:,I_RHOG)
       vz (:,:,:) = PROG(:,:,:,I_RHOGVZ) / PROG(:,:,:,I_RHOG)
       ein(:,:,:) = PROG(:,:,:,I_RHOGE ) / PROG(:,:,:,I_RHOG)

       do nq = 1, TRC_VMAX
          q(:,:,:,nq) = PROGq(:,:,:,nq) / PROG(:,:,:,I_RHOG)
       enddo

       cv(:,:,:)  = 0.D0
       qd(:,:,:)  = 1.D0
       do nq = NQW_STR, NQW_END
          cv(:,:,:) = cv(:,:,:) + q(:,:,:,nq) * CVW(nq)
          qd(:,:,:) = qd(:,:,:) - q(:,:,:,nq)
       enddo
       cv(:,:,:) = cv(:,:,:) + qd(:,:,:) * CNST_CV

       tem(:,:,:) = ein(:,:,:) / cv(:,:,:)
       pre(:,:,:) = rho(:,:,:) * tem(:,:,:) * ( qd(:,:,:)*CNST_RAIR + q(:,:,:,I_QV)*CNST_RVAP )

       do l = 1, ADM_lall
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall
             w(g,k,l) = PROG(g,k,l,I_RHOGW) / ( VMTR_C2Wfact(1,g,k,l) * PROG(g,k  ,l,I_RHOG) &
                                              + VMTR_C2Wfact(2,g,k,l) * PROG(g,k-1,l,I_RHOG) )
          enddo
          enddo

          !--- boundary conditions
          call bndcnd_all( ADM_gall,                 & ! [IN]
                           rho (:,:,l),              & ! [INOUT]
                           vx  (:,:,l),              & ! [INOUT]
                           vy  (:,:,l),              & ! [INOUT]
                           vz  (:,:,l),              & ! [INOUT]
                           w   (:,:,l),              & ! [INOUT]
                           ein (:,:,l),              & ! [INOUT]
                           tem (:,:,l),              & ! [INOUT]
                           pre (:,:,l),              & ! [INOUT]
                           PROG(:,:,l,I_RHOG),       & ! [INOUT]
                           PROG(:,:,l,I_RHOGVX),     & ! [INOUT]
                           PROG(:,:,l,I_RHOGVY),     & ! [INOUT]
                           PROG(:,:,l,I_RHOGVZ),     & ! [INOUT]
                           PROG(:,:,l,I_RHOGW),      & ! [INOUT]
                           PROG(:,:,l,I_RHOGE),      & ! [INOUT]
                           VMTR_GSGAM2    (:,:,l),   & ! [IN]
                           VMTR_PHI       (:,:,l),   & ! [IN]
                           VMTR_C2Wfact   (:,:,:,l), & ! [IN]
                           VMTR_C2Wfact_Gz(:,:,:,l)  ) ! [IN]

          call thrmdyn_th( ADM_gall, th(:,:,l), tem(:,:,l), pre(:,:,l) )

          call thrmdyn_eth( ADM_gall, eth(:,:,l), ein(:,:,l), pre(:,:,l), rho(:,:,l) )
       enddo ! region LOOP

       !--- perturbations ( pred, rhod, temd )
       pregd(:,:,:) = ( pre(:,:,:) - pre_bs(:,:,:) ) * VMTR_GSGAM2(:,:,:)
       rhogd(:,:,:) = ( rho(:,:,:) - rho_bs(:,:,:) ) * VMTR_GSGAM2(:,:,:)
       temd (:,:,:) =   tem(:,:,:) - tem_bs(:,:,:)

       if ( ADM_prc_me == ADM_prc_pl ) then

          rho_pl(:,:,:) = PROG_pl(:,:,:,I_RHOG  ) / VMTR_GSGAM2_pl(:,:,:)
          vx_pl (:,:,:) = PROG_pl(:,:,:,I_RHOGVX) / PROG_pl(:,:,:,I_RHOG)
          vy_pl (:,:,:) = PROG_pl(:,:,:,I_RHOGVY) / PROG_pl(:,:,:,I_RHOG)
          vz_pl (:,:,:) = PROG_pl(:,:,:,I_RHOGVZ) / PROG_pl(:,:,:,I_RHOG)
          ein_pl(:,:,:) = PROG_pl(:,:,:,I_RHOGE ) / PROG_pl(:,:,:,I_RHOG)

          do nq = 1, TRC_VMAX
             q_pl(:,:,:,nq) = PROGq_pl(:,:,:,nq) / PROG_pl(:,:,:,I_RHOG)
          enddo

          cv_pl(:,:,:)  = 0.D0
          qd_pl(:,:,:)  = 1.D0
          do nq = NQW_STR, NQW_END
             cv_pl(:,:,:) = cv_pl(:,:,:) + q_pl(:,:,:,nq) * CVW(nq)
             qd_pl(:,:,:) = qd_pl(:,:,:) - q_pl(:,:,:,nq)
          enddo
          cv_pl(:,:,:) = cv_pl(:,:,:) + qd_pl(:,:,:) * CNST_CV

          tem_pl(:,:,:) = ein_pl(:,:,:) / cv_pl(:,:,:)
          pre_pl(:,:,:) = rho_pl(:,:,:) * tem_pl(:,:,:) * ( qd_pl(:,:,:)*CNST_RAIR + q_pl(:,:,:,I_QV)*CNST_RVAP )

          do l = 1, ADM_lall_pl
             do k = ADM_kmin+1, ADM_kmax
             do g = 1, ADM_gall_pl
                w_pl(g,k,l) = PROG_pl(g,k,l,I_RHOGW) / ( VMTR_C2Wfact_pl(1,g,k,l) * PROG_pl(g,k  ,l,I_RHOG) &
                                                       + VMTR_C2Wfact_pl(2,g,k,l) * PROG_pl(g,k-1,l,I_RHOG) )
             enddo
             enddo

             !--- boundary conditions
             call bndcnd_all( ADM_gall_pl,                 & ! [IN]
                              rho_pl (:,:,l),              & ! [INOUT]
                              vx_pl  (:,:,l),              & ! [INOUT]
                              vy_pl  (:,:,l),              & ! [INOUT]
                              vz_pl  (:,:,l),              & ! [INOUT]
                              w_pl   (:,:,l),              & ! [INOUT]
                              ein_pl (:,:,l),              & ! [INOUT]
                              tem_pl (:,:,l),              & ! [INOUT]
                              pre_pl (:,:,l),              & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOG  ),     & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGVX),     & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGVY),     & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGVZ),     & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGW ),     & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGE ),     & ! [INOUT]
                              VMTR_GSGAM2_pl    (:,:,l),   & ! [IN]
                              VMTR_PHI_pl       (:,:,l),   & ! [IN]
                              VMTR_C2Wfact_pl   (:,:,:,l), & ! [IN]
                              VMTR_C2Wfact_Gz_pl(:,:,:,l)  ) ! [IN]

             call thrmdyn_th ( ADM_gall_pl, th_pl (:,:,l), tem_pl(:,:,l), pre_pl(:,:,l) )
             call thrmdyn_eth( ADM_gall_pl, eth_pl(:,:,l), ein_pl(:,:,l), pre_pl(:,:,l), rho_pl(:,:,l) )
          enddo

          pregd_pl(:,:,:) = ( pre_pl(:,:,:) - pre_bs_pl(:,:,:) ) * VMTR_GSGAM2_pl(:,:,:)
          rhogd_pl(:,:,:) = ( rho_pl(:,:,:) - rho_bs_pl(:,:,:) ) * VMTR_GSGAM2_pl(:,:,:)
          temd_pl (:,:,:) =   tem_pl(:,:,:) - tem_bs_pl(:,:,:)
       else

          rho_pl(:,:,:) = 0.D0
          vx_pl (:,:,:) = 0.D0
          vy_pl (:,:,:) = 0.D0
          vz_pl (:,:,:) = 0.D0
          w_pl  (:,:,:) = 0.D0
          ein_pl(:,:,:) = 0.D0

          q_pl  (:,:,:,:) = 0.D0

          tem_pl(:,:,:) = 0.D0
          pre_pl(:,:,:) = 0.D0
          th_pl (:,:,:) = 0.D0
          eth_pl(:,:,:) = 0.D0

          pregd_pl(:,:,:) = 0.D0
          rhogd_pl(:,:,:) = 0.D0
          temd_pl (:,:,:) = 0.D0

       endif

       !------------------------------------------------------------------------
       !> LARGE step
       !------------------------------------------------------------------------
       call DEBUG_rapstart('+++Large step')

       !--- calculation of advection tendency including Coriolis force
       call src_advection_convergence_momentum( vx,                     vx_pl,                     & ! [IN]
                                                vy,                     vy_pl,                     & ! [IN]
                                                vz,                     vz_pl,                     & ! [IN]
                                                w,                      w_pl,                      & ! [IN]
                                                PROG  (:,:,:,I_RHOG  ), PROG_pl  (:,:,:,I_RHOG  ), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & ! [IN]
                                                PROG  (:,:,:,I_RHOGW ), PROG_pl  (:,:,:,I_RHOGW ), & ! [IN]
                                                g_TEND(:,:,:,I_RHOGVX), g_TEND_pl(:,:,:,I_RHOGVX), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGVY), g_TEND_pl(:,:,:,I_RHOGVY), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGVZ), g_TEND_pl(:,:,:,I_RHOGVZ), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGW ), g_TEND_pl(:,:,:,I_RHOGW )  ) ! [OUT]

       g_TEND   (:,:,:,I_RHOG)     = 0.D0
       g_TEND   (:,:,:,I_RHOGE)    = 0.D0
       g_TEND   (:,:,:,I_RHOGETOT) = 0.D0
       g_TEND_pl(:,:,:,I_RHOG)     = 0.D0
       g_TEND_pl(:,:,:,I_RHOGE)    = 0.D0
       g_TEND_pl(:,:,:,I_RHOGETOT) = 0.D0

       !---< numerical diffusion term
       if ( NDIFF_LOCATION == 'IN_LARGE_STEP' ) then

          if ( nl == 1 ) then ! only first step
             f_TEND    (:,:,:,:) = 0.D0
             f_TEND_pl (:,:,:,:) = 0.D0
             f_TENDq   (:,:,:,:) = 0.D0
             f_TENDq_pl(:,:,:,:) = 0.D0

             !------ numerical diffusion
             call numfilter_hdiffusion( rho,                       rho_pl,                       & ! [IN]
                                        vx,                        vx_pl,                        & ! [IN]
                                        vy,                        vy_pl,                        & ! [IN]
                                        vz,                        vz_pl,                        & ! [IN]
                                        w,                         w_pl,                         & ! [IN]
                                        temd,                      temd_pl,                      & ! [IN]
                                        q,                         q_pl,                         & ! [IN]
                                        f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & ! [INOUT]
                                        f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) ! [INOUT]

             if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
                call numfilter_vdiffusion( rho,                       rho_pl,                       & ! [IN]
                                           vx,                        vx_pl,                        & ! [IN]
                                           vy,                        vy_pl,                        & ! [IN]
                                           vz,                        vz_pl,                        & ! [IN]
                                           w,                         w_pl,                         & ! [IN]
                                           temd,                      temd_pl,                      & ! [IN]
                                           q,                         q_pl,                         & ! [IN]
                                           f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & ! [INOUT]
                                           f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & ! [INOUT]
                                           f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & ! [INOUT]
                                           f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & ! [INOUT]
                                           f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & ! [INOUT]
                                           f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & ! [INOUT]
                                           f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & ! [INOUT]
                                           f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) ! [INOUT]
             endif

             if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
                call numfilter_rayleigh_damping( rho,                       rho_pl,                       & ! [IN]
                                                 vx,                        vx_pl,                        & ! [IN]
                                                 vy,                        vy_pl,                        & ! [IN]
                                                 vz,                        vz_pl,                        & ! [IN]
                                                 w,                         w_pl,                         & ! [IN]
                                                 f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & ! [INOUT]
                                                 f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & ! [INOUT]
                                                 f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & ! [INOUT]
                                                 f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   )  ) ! [INOUT]
             endif

          endif

       elseif( NDIFF_LOCATION == 'IN_LARGE_STEP2' ) then
          f_TEND    (:,:,:,:) = 0.D0
          f_TEND_pl (:,:,:,:) = 0.D0
          f_TENDq   (:,:,:,:) = 0.D0
          f_TENDq_pl(:,:,:,:) = 0.D0

          !------ numerical diffusion
          call numfilter_hdiffusion( rho,                       rho_pl,                       & ! [IN]
                                     vx,                        vx_pl,                        & ! [IN]
                                     vy,                        vy_pl,                        & ! [IN]
                                     vz,                        vz_pl,                        & ! [IN]
                                     w,                         w_pl,                         & ! [IN]
                                     temd,                      temd_pl,                      & ! [IN]
                                     q,                         q_pl,                         & ! [IN]
                                     f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & ! [INOUT]
                                     f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & ! [INOUT]
                                     f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & ! [INOUT]
                                     f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & ! [INOUT]
                                     f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & ! [INOUT]
                                     f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & ! [INOUT]
                                     f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & ! [INOUT]
                                     f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) ! [INOUT]

          if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
             call numfilter_vdiffusion( rho,                       rho_pl,                       & ! [IN]
                                        vx,                        vx_pl,                        & ! [IN]
                                        vy,                        vy_pl,                        & ! [IN]
                                        vz,                        vz_pl,                        & ! [IN]
                                        w,                         w_pl,                         & ! [IN]
                                        temd,                      temd_pl,                      & ! [IN]
                                        q,                         q_pl,                         & ! [IN]
                                        f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & ! [INOUT]
                                        f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & ! [INOUT]
                                        f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) ! [INOUT]
          endif

          if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
             call numfilter_rayleigh_damping( rho,                       rho_pl,                       & ! [IN]
                                              vx,                        vx_pl,                        & ! [IN]
                                              vy,                        vy_pl,                        & ! [IN]
                                              vz,                        vz_pl,                        & ! [IN]
                                              w,                         w_pl,                         & ! [IN]
                                              f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & ! [INOUT]
                                              f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & ! [INOUT]
                                              f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & ! [INOUT]
                                              f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   )  ) ! [INOUT]
          endif
       endif

       if ( TB_TYPE == 'SMG' ) then ! Smagorinksy-type SGS model [add] A.Noda 10/11/29
          call sgs_smagorinsky( nl,                                                      &
                                rho,                       rho_pl,                       & ! [IN]
                                PROG(:,:,:,I_RHOG  ),      PROG_pl(:,:,:,I_RHOG  ),      & ! [IN]
                                PROGq(:,:,:,:),            PROGq_pl(:,:,:,:  ),          & ! [IN]
                                vx,                        vx_pl,                        & ! [IN]
                                vy,                        vy_pl,                        & ! [IN]
                                vz,                        vz_pl,                        & ! [IN]
                                w,                         w_pl,                         & ! [IN]
                                tem,                       tem_pl,                       & ! [IN]
                                q,                         q_pl,                         & ! [IN]
                                th,                        th_pl,                        & ! [IN]
                                f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & ! [INOUT]
                                f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & ! [INOUT]
                                f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & ! [INOUT]
                                f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & ! [INOUT]
                                f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & ! [INOUT]
                                f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & ! [INOUT]
                                f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & ! [INOUT]
                                f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) ! [INOUT]
       endif

       !--- Nudging routines [add] Y.Niwa 08/09/09
       if ( FLAG_NUDGING ) then

          if ( nl == 1 ) then
             call NDG_update_reference( TIME_CTIME )
          endif

          if ( nl == num_of_iteration_lstep ) then
             ndg_TEND_out = .true.
          else
             ndg_TEND_out = .false.
          endif

          call NDG_apply_uvtp( rho,                      rho_pl,                      & ! [IN]
                               vx,                       vx_pl,                       & ! [IN]
                               vy,                       vy_pl,                       & ! [IN]
                               vz,                       vz_pl,                       & ! [IN]
                               w,                        w_pl,                        & ! [IN]
                               tem,                      tem_pl,                      & ! [IN]
                               pre,                      pre_pl,                      & ! [IN]
                               f_TEND(:,:,:,I_RHOGVX  ), f_TEND_pl(:,:,:,I_RHOGVX  ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGVY  ), f_TEND_pl(:,:,:,I_RHOGVY  ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGVZ  ), f_TEND_pl(:,:,:,I_RHOGVZ  ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGW   ), f_TEND_pl(:,:,:,I_RHOGW   ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGE   ), f_TEND_pl(:,:,:,I_RHOGE   ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGETOT), f_TEND_pl(:,:,:,I_RHOGETOT), & ! [INOUT]
                               ndg_TEND_out                                           ) ! [IN] ( TEND out )
       endif

       !--- sum the large step TEND ( advection + coriolis + num.diff.,SGS,nudge )
       g_TEND   (:,:,:,:) = g_TEND   (:,:,:,:) + f_TEND   (:,:,:,:)
       g_TEND_pl(:,:,:,:) = g_TEND_pl(:,:,:,:) + f_TEND_pl(:,:,:,:)

       call DEBUG_rapend  ('+++Large step')
       !------------------------------------------------------------------------
       !> SMALL step
       !------------------------------------------------------------------------
       call DEBUG_rapstart('+++Small step')

       if ( nl /= 1 ) then ! update split values
          PROG_split   (:,:,:,:) = PROG0   (:,:,:,:) - PROG   (:,:,:,:)
          PROG_split_pl(:,:,:,:) = PROG0_pl(:,:,:,:) - PROG_pl(:,:,:,:)
       else
          PROG_split   (:,:,:,:) = 0.D0
          PROG_split_pl(:,:,:,:) = 0.D0
       endif

       !------ Core routine for small step
       !------    1. By this subroutine, prognostic variables ( rho,.., rhoge ) are calculated through
       !------       the 'num_of_iteration_sstep(nl)'-th times small step.
       !------    2. grho, grhogvx, ..., and  grhoge has the large step
       !------       tendencies initially, however, they are re-used in this subroutine.
       !------
       if ( TIME_SPLIT ) then
          small_step_ite = num_of_iteration_sstep(nl)
          small_step_dt  = TIME_DTS
       else
          small_step_ite = 1
          small_step_dt  = TIME_DTL / (num_of_iteration_lstep+1-nl)
       endif

       call vi_small_step( PROG(:,:,:,I_RHOG  ),       PROG_pl(:,:,:,I_RHOG  ),       & ! [INOUT] prognostic variables
                           PROG(:,:,:,I_RHOGVX),       PROG_pl(:,:,:,I_RHOGVX),       & ! [INOUT]
                           PROG(:,:,:,I_RHOGVY),       PROG_pl(:,:,:,I_RHOGVY),       & ! [INOUT]
                           PROG(:,:,:,I_RHOGVZ),       PROG_pl(:,:,:,I_RHOGVZ),       & ! [INOUT]
                           PROG(:,:,:,I_RHOGW ),       PROG_pl(:,:,:,I_RHOGW ),       & ! [INOUT]
                           PROG(:,:,:,I_RHOGE ),       PROG_pl(:,:,:,I_RHOGE ),       & ! [INOUT]
                           vx,                         vx_pl,                         & ! [IN] diagnostic value
                           vy,                         vy_pl,                         & ! [IN]
                           vz,                         vz_pl,                         & ! [IN]
                           eth,                        eth_pl,                        & ! [IN]
                           rhogd,                      rhogd_pl,                      & ! [IN]
                           pregd,                      pregd_pl,                      & ! [IN]
                           g_TEND(:,:,:,I_RHOG    ),   g_TEND_pl(:,:,:,I_RHOG    ),   & ! [IN] large step TEND
                           g_TEND(:,:,:,I_RHOGVX  ),   g_TEND_pl(:,:,:,I_RHOGVX  ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGVY  ),   g_TEND_pl(:,:,:,I_RHOGVY  ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGVZ  ),   g_TEND_pl(:,:,:,I_RHOGVZ  ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGW   ),   g_TEND_pl(:,:,:,I_RHOGW   ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGE   ),   g_TEND_pl(:,:,:,I_RHOGE   ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGETOT),   g_TEND_pl(:,:,:,I_RHOGETOT),   & ! [IN]
                           PROG_split(:,:,:,I_RHOG  ), PROG_split_pl(:,:,:,I_RHOG  ), & ! [INOUT] split value
                           PROG_split(:,:,:,I_RHOGVX), PROG_split_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                           PROG_split(:,:,:,I_RHOGVY), PROG_split_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                           PROG_split(:,:,:,I_RHOGVZ), PROG_split_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                           PROG_split(:,:,:,I_RHOGW ), PROG_split_pl(:,:,:,I_RHOGW ), & ! [INOUT]
                           PROG_split(:,:,:,I_RHOGE ), PROG_split_pl(:,:,:,I_RHOGE ), & ! [INOUT]
                           v_mean_c,                   v_mean_c_pl,                   & ! [OUT] mean value
                           small_step_ite,                                            & ! [IN]
                           small_step_dt                                              ) ! [IN]

       call DEBUG_rapend  ('+++Small step')

       !------------------------------------------------------------------------
       !>  Tracer advection
       !------------------------------------------------------------------------
       call DEBUG_rapstart('+++Tracer Advection')
       do_tke_correction = .false.

       if ( TRC_ADV_TYPE == 'MIURA2004' ) then

          if ( nl == num_of_iteration_lstep ) then

             call src_tracer_advection( TRC_VMAX,                                              & ! [IN]
                                     PROGq(:,:,:,:),           PROGq_pl(:,:,:,:),           & ! [INOUT]
                                     PROG0(:,:,:,I_RHOG),      PROG0_pl(:,:,:,I_RHOG),      & ! [IN]
                                     v_mean_c(:,:,:,I_rhog),   v_mean_c_pl(:,:,:,I_rhog),   & ! [IN]
                                     v_mean_c(:,:,:,I_rhogvx), v_mean_c_pl(:,:,:,I_rhogvx), & ! [IN]
                                     v_mean_c(:,:,:,I_rhogvy), v_mean_c_pl(:,:,:,I_rhogvy), & ! [IN]
                                     v_mean_c(:,:,:,I_rhogvz), v_mean_c_pl(:,:,:,I_rhogvz), & ! [IN]
                                     v_mean_c(:,:,:,I_rhogw),  v_mean_c_pl(:,:,:,I_rhogw),  & ! [IN]
                                     f_TEND (:,:,:,I_RHOG),    f_TEND_pl (:,:,:,I_RHOG),    & ! [IN]
                                     TIME_DTL,                                              & ! [IN]
                                     THUBURN_LIM                                            ) ! [IN]  ![add] 20130613 R.Yoshida

             PROGq(:,:,:,:) = PROGq(:,:,:,:) + TIME_DTL * f_TENDq(:,:,:,:) ! update rhogq by viscosity

             PROGq(:,ADM_kmin-1,:,:) = 0.D0
             PROGq(:,ADM_kmax+1,:,:) = 0.D0

             if ( ADM_prc_pl == ADM_prc_me ) then
                PROGq_pl(:,:,:,:) = PROGq_pl(:,:,:,:) + TIME_DTL * f_TENDq_pl(:,:,:,:)

                PROGq_pl(:,ADM_kmin-1,:,:) = 0.D0
                PROGq_pl(:,ADM_kmax+1,:,:) = 0.D0
             endif

             ! [comment] H.Tomita: I don't recommend adding the hyperviscosity term because of numerical instability in this case.
             if( I_TKE >= 0 ) do_tke_correction = .true.

          endif ! Last large step only

       elseif( TRC_ADV_TYPE == 'DEFAULT' ) then

          do nq = 1, TRC_VMAX

             call src_advection_convergence( v_mean_c(:,:,:,I_rhogvx), v_mean_c_pl(:,:,:,I_rhogvx), & ! [IN]
                                             v_mean_c(:,:,:,I_rhogvy), v_mean_c_pl(:,:,:,I_rhogvy), & ! [IN]
                                             v_mean_c(:,:,:,I_rhogvz), v_mean_c_pl(:,:,:,I_rhogvz), & ! [IN]
                                             v_mean_c(:,:,:,I_rhogw),  v_mean_c_pl(:,:,:,I_rhogw),  & ! [IN]
                                             q(:,:,:,nq),              q_pl(:,:,:,nq),              & ! [IN]
                                             g_TENDq(:,:,:,nq),        g_TENDq_pl(:,:,:,nq),        & ! [OUT]
                                             I_SRC_default                                          ) ! [IN]  [mod] H.Yashiro 20120530

             PROGq(:,:,:,:) = PROGq0(:,:,:,:)                                                                   &
                            + ( num_of_iteration_sstep(nl) * TIME_DTS ) * ( g_TENDq(:,:,:,:) + f_TENDq(:,:,:,:) )

             PROGq(:,ADM_kmin-1,:,:) = 0.D0
             PROGq(:,ADM_kmax+1,:,:) = 0.D0

             if ( ADM_prc_pl == ADM_prc_me ) then
                      PROGq_pl(:,:,:,:) = PROGq0_pl(:,:,:,:)                          &
                                        + ( num_of_iteration_sstep(nl) * TIME_DTS )   &
                                        * ( g_TENDq_pl(:,:,:,:) + f_TENDq_pl(:,:,:,:) )

                PROGq_pl(:,ADM_kmin-1,:,:) = 0.D0
                PROGq_pl(:,ADM_kmax+1,:,:) = 0.D0
             endif

          enddo ! tracer LOOP

          if( I_TKE >= 0 ) do_tke_correction = .true.

       endif

       call DEBUG_rapend  ('+++Tracer Advection')

       !--- TKE fixer [comment] 2011/08/16 M.Satoh: need this fixer for every small time steps
       if ( do_tke_correction ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             TKEG_corr = max( TKE_MIN * VMTR_GSGAM2(g,k,l) - PROGq(g,k,l,I_TKE), 0.D0 )

             PROG (g,k,l,I_RHOGE) = PROG (g,k,l,I_RHOGE) - TKEG_corr
             PROGq(g,k,l,I_TKE)   = PROGq(g,k,l,I_TKE)   + TKEG_corr
          enddo
          enddo
          enddo

          if ( ADM_prc_pl == ADM_prc_me ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                TKEg_corr = max( TKE_MIN * VMTR_GSGAM2_pl(g,k,l) - PROGq_pl(g,k,l,I_TKE), 0.D0 )

                PROG_pl (g,k,l,I_RHOGE) = PROG_pl (g,k,l,I_RHOGE) - TKEG_corr
                PROGq_pl(g,k,l,I_TKE)   = PROGq_pl(g,k,l,I_TKE)   + TKEG_corr
             enddo
             enddo
             enddo
          endif
       endif

       !------ Update
       if ( nl /= num_of_iteration_lstep ) then
          call COMM_data_transfer( PROG, PROG_pl )

          PROG(suf(ADM_gall_1d,1),:,:,:) = PROG(suf(ADM_gmax+1,ADM_gmin),:,:,:)
          PROG(suf(1,ADM_gall_1d),:,:,:) = PROG(suf(ADM_gmin,ADM_gmax+1),:,:,:)
       endif

    enddo !--- large step

    call prgvar_set( PROG(:,:,:,I_RHOG),   PROG_pl(:,:,:,I_RHOG),   & ! [IN]
                     PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & ! [IN]
                     PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & ! [IN]
                     PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ), & ! [IN]
                     PROG(:,:,:,I_RHOGW),  PROG_pl(:,:,:,I_RHOGW),  & ! [IN]
                     PROG(:,:,:,I_RHOGE),  PROG_pl(:,:,:,I_RHOGE),  & ! [IN]
                     PROGq(:,:,:,:),       PROGq_pl(:,:,:,:),       & ! [IN]
                     0                                              ) ! [IN]

    call DEBUG_rapend  ('++Dynamics')

    return
  end subroutine dynamics_step

end module mod_dynamics
