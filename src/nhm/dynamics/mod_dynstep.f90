!-------------------------------------------------------------------------------
!
!+  Dynamical step
!
!-------------------------------------------------------------------------------
module mod_dynstep
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module is for the dynamical step 
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                06-04-17   Add IN_LARGE_STEP2
  !                06-08-11   Add the option for tracer advection.
  !                07-01-26   Add flag [rayleigh_damp_only_w] 
  !                           in numfilter_rayleigh_damping.
  !                07-05-08   H.Tomita : Change the treatment of I_TKE.
  !                08-01-24   Y.Niwa: add revised MIURA2004 for tracer advection
  !                           old: 'MIURA2004OLD', revised: 'MIURA2004'
  !                08-01-30   Y.Niwa: add rho_pl = 0.D0
  !                08-04-12   T.Mitsui save memory(prgvar, frcvar, rhog0xxxx)
  !                08-05-24   T.Mitsui fix miss-conditioning for frcvar
  !                08-09-09   Y.Niwa move nudging routine here
  !                08-10-05   T.Mitsui all_phystep_post is already needless
  !                09-09-08   S.Iga  frhog and frhog_pl in ndg are deleted ( suggested by ES staff)
  !                10-05-06   M.Satoh: define QV_conv only if CP_TYPE='TDK' .or. 'KUO'
  !                10-07-16   A.T.Noda: bug fix for TDK
  !                10-08-16   A.T.Noda: Bug fix (Qconv not diveded by density)
  !                10-08-20   A.T.Noda: Bug fix (Qconv should be TEND, and not be multiplied by DT)
  !                10-11-29   A.T.Noda: Introduce the Smagorinsky model
  !                11-08-16   M.Satoh: bug fix for TDK: conv => TEND
  !                           qv_dyn_tend = v grad q
  !                                       = ( div(rho v q) - div(rho v)*q )/rho
  !                11-08-16   M.Satoh: move codes related to CP_TYPE below the tracer calculation
  !                11-11-28   Y.Yamada: Merge Terai-san timer into the original code.
  !                12-03-09   S.Iga: tuned (phase4-1)
  !                12-04-06   T.yamaura: optimized for K
  !                12-05-30   T.Yashiro: Change arguments from character to index/switch
  !                12-10-22   R.Yoshida  : add papi instructions
  !      -----------------------------------------------------------------------
  !
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
  public :: dynstep
  public :: dynstep2

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


  !-----------------------------------------------------------------------------
contains 
  !-----------------------------------------------------------------------------
  subroutine dynstep
    use mod_debug
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
    use mod_time, only:  &
       TIME_INTEG_TYPE,  &
       TIME_SSTEP_MAX,   &
       TIME_DTL,         &
       TIME_DTS,         &
       TIME_SPLIT
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_GSGAM2,     &
       VMTR_GSGAM2_pl,  &
       VMTR_GSGAM2H,    &
       VMTR_GSGAM2H_pl, &
       VMTR_GZXH,       &
       VMTR_GZXH_pl,    &
       VMTR_GZYH,       &
       VMTR_GZYH_pl,    &
       VMTR_GZZH,       &
       VMTR_GZZH_pl,    &
       VMTR_PHI,        &
       VMTR_PHI_pl,     &
       VMTR_C2Wfact,    &
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
       FLAG_NUDGING,   & ! Y.Niwa add 08/09/09
       CP_TYPE,        & ! 2010.5.11 M.Satoh [add]
       TB_TYPE           ! [add] 10/11/29 A.Noda
    use mod_bsstate, only: &
       pre_bs, pre_bs_pl, &
       tem_bs, tem_bs_pl, &
       rho_bs, rho_bs_pl
    use mod_bndcnd, only: &
       bndcnd_all
    use mod_prgvar, only: &
       prgvar_set,    &
       prgvar_get,    &
       prgvar_get_noq
    use mod_diagvar, only: &
       diagvar,       &
       diagvar_pl,    &
       I_RHOGQV_CONV, &
       I_QV_DYN_TEND    ! 2011.08.16 M.Satoh
    use mod_thrmdyn, only: &
       thrmdyn_th, &
       thrmdyn_eth
    use mod_numfilter, only: &
       NUMFILTER_DOrayleigh,       & ! [add] H.Yashiro 20120530
       NUMFILTER_DOverticaldiff,   & ! [add] H.Yashiro 20120530
       numfilter_rayleigh_damping, &
       numfilter_numerical_hdiff,  &
       numfilter_numerical_vdiff
    use mod_vi, only :         &
       vi_small_step
    use mod_src, only: &
       src_advection_convergence_v, &
       src_advection_convergence,   &
       src_update_tracer,           &
       I_SRC_default                  ! [add] H.Yashiro 20120530
    use mod_ndg, only: & ! Y.Niwa add 08/09/09
       ndg_nudging_uvtp, &
       ndg_update_var
    use mod_tb_smg, only: & ! [add] 10/11/29 A.Noda
       tb_smg_driver
    implicit none

    integer, parameter :: nmax_TEND   = 7
    integer, parameter :: nmax_PROG = 6
    integer, parameter :: nmax_v_mean_c   = 5

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

    !--- [IN]ternal energy  ( physical )
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
    real(8) :: qd      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: qd_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: cv      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: cv_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), parameter :: TKE_MIN = 0.D0
    real(8)            :: TKEg_corr

    integer :: small_step_ite
    real(8) :: small_step_dt

    logical :: ndg_TEND_out

    logical, save :: iflag = .true.
    integer, save :: num_of_iteration_lstep    ! number of large steps ( 2-4 )
    integer, save :: num_of_iteration_sstep(4) ! number of small steps in each of large steps

    integer :: ij, k ,l, nq, nl

    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

#ifdef PAPI_OPS
    ! <-- [add] PAPI R.Yoshida 20121022
    !call PAPIF_flips( PAPI_real_time_i, PAPI_proc_time_i, PAPI_flpins, PAPI_mflins, PAPI_check )
    call PAPIF_flops( PAPI_real_time_o, PAPI_proc_time_o, PAPI_flpops, PAPI_mflops, PAPI_check )
#endif

    if ( iflag ) then
       iflag = .false.

       select case(trim(TIME_INTEG_TYPE))
       case('RK2')
          num_of_iteration_lstep = 2
          num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 2
          num_of_iteration_sstep(2) = TIME_SSTEP_MAX
       case('RK3')
          num_of_iteration_lstep = 3
          num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 3
          num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 2
          num_of_iteration_sstep(3) = TIME_SSTEP_MAX
        case('RK4')
          num_of_iteration_lstep = 4
          num_of_iteration_sstep(1) = TIME_SSTEP_MAX/4
          num_of_iteration_sstep(2) = TIME_SSTEP_MAX/3
          num_of_iteration_sstep(3) = TIME_SSTEP_MAX/2
          num_of_iteration_sstep(4) = TIME_SSTEP_MAX
       case default
          write(*,*) 'Msg : Sub[sub_dynstep]'
          write(*,*) ' --- Error : invalid TIME_INTEG_TYPE=', TIME_INTEG_TYPE
       endselect
    endif

    !--- get from prg0
    call prgvar_get( PROG(:,:,:,I_RHOG),   PROG_pl(:,:,:,I_RHOG),   & !--- [OUT]
                     PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & !--- [OUT]
                     PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & !--- [OUT]
                     PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ), & !--- [OUT]
                     PROG(:,:,:,I_RHOGW),  PROG_pl(:,:,:,I_RHOGW),  & !--- [OUT]
                     PROG(:,:,:,I_RHOGE),  PROG_pl(:,:,:,I_RHOGE),  & !--- [OUT]
                     PROGq(:,:,:,:),       PROGq_pl(:,:,:,:),       & !--- [OUT]
                     0                                              ) !--- [IN]

    !--- save
    PROG0   (:,:,:,:) = PROG   (:,:,:,:)
    PROG0_pl(:,:,:,:) = PROG_pl(:,:,:,:)

    if ( TRC_ADV_TYPE == 'DEFAULT' ) then
       PROGq0   (:,:,:,:) = PROGq   (:,:,:,:)
       PROGq0_pl(:,:,:,:) = PROGq_pl(:,:,:,:)
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
             do ij = 1, ADM_gall
                w(ij,k,l) = PROG(ij,k,l,I_RHOGW) &
                          / ( VMTR_GSGAM2H(ij,k,l) * 0.5D0 * ( GRD_afac(k) * rho(ij,k  ,l) &
                                                             + GRD_bfac(k) * rho(ij,k-1,l) ) )
             enddo
          enddo

          !--- boundary conditions
          call bndcnd_all( ADM_gall,             & !--- [IN]
                           rho   (:,:,l),        & !--- [INOUT]
                           vx    (:,:,l),        & !--- [INOUT]
                           vy    (:,:,l),        & !--- [INOUT]
                           vz    (:,:,l),        & !--- [INOUT]
                           w     (:,:,l),        & !--- [INOUT]
                           ein   (:,:,l),        & !--- [INOUT]
                           tem   (:,:,l),        & !--- [INOUT]
                           pre   (:,:,l),        & !--- [INOUT]
                           PROG(:,:,l,I_RHOG),   & !--- [INOUT]
                           PROG(:,:,l,I_RHOGVX), & !--- [INOUT]
                           PROG(:,:,l,I_RHOGVY), & !--- [INOUT]
                           PROG(:,:,l,I_RHOGVZ), & !--- [INOUT]
                           PROG(:,:,l,I_RHOGW),  & !--- [INOUT]
                           PROG(:,:,l,I_RHOGE),  & !--- [INOUT]
                           VMTR_GSGAM2 (:,:,l),  & !--- [IN]
                           VMTR_GSGAM2H(:,:,l),  & !--- [IN]
                           VMTR_PHI    (:,:,l),  & !--- [IN]
                           VMTR_C2Wfact(:,:,:,l) ) !--- [IN]

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
                do ij = 1, ADM_gall_pl
                   w_pl(ij,k,l) = PROG_pl(ij,k,l,I_RHOGW) &
                             / ( 0.5D0 * ( GRD_afac(k) * rho_pl(ij,k  ,l) &
                                         + GRD_bfac(k) * rho_pl(ij,k-1,l) ) * VMTR_GSGAM2H_pl(ij,k,l) )
                enddo
             enddo

             !--- boundary conditions
             call bndcnd_all( ADM_gall_pl,            & !--- [IN]
                              rho_pl   (:,:,l),       & !--- [INOUT]
                              vx_pl    (:,:,l),       & !--- [INOUT]
                              vy_pl    (:,:,l),       & !--- [INOUT]
                              vz_pl    (:,:,l),       & !--- [INOUT]
                              w_pl     (:,:,l),       & !--- [INOUT]
                              ein_pl   (:,:,l),       & !--- [INOUT]
                              tem_pl   (:,:,l),       & !--- [INOUT]
                              pre_pl   (:,:,l),       & !--- [INOUT]
                              PROG_pl(:,:,l,I_RHOG),   & !--- [INOUT]
                              PROG_pl(:,:,l,I_RHOGVX), & !--- [INOUT]
                              PROG_pl(:,:,l,I_RHOGVY), & !--- [INOUT]
                              PROG_pl(:,:,l,I_RHOGVZ), & !--- [INOUT]
                              PROG_pl(:,:,l,I_RHOGW),  & !--- [INOUT]
                              PROG_pl(:,:,l,I_RHOGE),  & !--- [INOUT]
                              VMTR_GSGAM2_pl (:,:,l),   & !--- [IN]
                              VMTR_GSGAM2H_pl(:,:,l),   & !--- [IN]
                              VMTR_PHI_pl    (:,:,l),   & !--- [IN]
                              VMTR_C2Wfact_pl(:,:,:,l)  ) !--- [IN]

             call thrmdyn_th( ADM_gall_pl, th_pl(:,:,l), tem_pl(:,:,l), pre_pl(:,:,l) )

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

       !--- calculation of advection tendency including Coriolis force
       call src_advection_convergence_v( vx,                     vx_pl,                     & !--- [IN]
                                         vy,                     vy_pl,                     & !--- [IN]
                                         vz,                     vz_pl,                     & !--- [IN]
                                         w,                      w_pl,                      & !--- [IN]
                                         PROG  (:,:,:,I_RHOG  ), PROG_pl  (:,:,:,I_RHOG  ), & !--- [IN]
                                         PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & !--- [IN]
                                         PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & !--- [IN]
                                         PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & !--- [IN]
                                         PROG  (:,:,:,I_RHOGW ), PROG_pl  (:,:,:,I_RHOGW ), & !--- [IN]
                                         g_TEND(:,:,:,I_RHOGVX), g_TEND_pl(:,:,:,I_RHOGVX), & !--- [OUT]
                                         g_TEND(:,:,:,I_RHOGVY), g_TEND_pl(:,:,:,I_RHOGVY), & !--- [OUT]
                                         g_TEND(:,:,:,I_RHOGVZ), g_TEND_pl(:,:,:,I_RHOGVZ), & !--- [OUT]
                                         g_TEND(:,:,:,I_RHOGW ), g_TEND_pl(:,:,:,I_RHOGW )  ) !--- [OUT]

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
             call numfilter_numerical_hdiff( rho,                       rho_pl,                       & !--- [IN]
                                             vx,                        vx_pl,                        & !--- [IN]
                                             vy,                        vy_pl,                        & !--- [IN]
                                             vz,                        vz_pl,                        & !--- [IN]
                                             w,                         w_pl,                         & !--- [IN]
                                             temd,                      temd_pl,                      & !--- [IN]
                                             q,                         q_pl,                         & !--- [IN]
                                             f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & !--- [INOUT]
                                             f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) !--- [INOUT]

             if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
                call numfilter_numerical_vdiff( rho,                       rho_pl,                       & !--- [IN]
                                                vx,                        vx_pl,                        & !--- [IN]
                                                vy,                        vy_pl,                        & !--- [IN]
                                                vz,                        vz_pl,                        & !--- [IN]
                                                w,                         w_pl,                         & !--- [IN]
                                                temd,                      temd_pl,                      & !--- [IN]
                                                q,                         q_pl,                         & !--- [IN]
                                                f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & !--- [INOUT]
                                                f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & !--- [INOUT]
                                                f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & !--- [INOUT]
                                                f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & !--- [INOUT]
                                                f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & !--- [INOUT]
                                                f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & !--- [INOUT]
                                                f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & !--- [INOUT]
                                                f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) !--- [INOUT]
             endif

             if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
                call numfilter_rayleigh_damping( rho,                       rho_pl,                       & !--- [IN]
                                                 vx,                        vx_pl,                        & !--- [IN]
                                                 vy,                        vy_pl,                        & !--- [IN]
                                                 vz,                        vz_pl,                        & !--- [IN]
                                                 w,                         w_pl,                         & !--- [IN]
                                                 f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & !--- [INOUT]
                                                 f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & !--- [INOUT]
                                                 f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & !--- [INOUT]
                                                 f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   )  ) !--- [INOUT]
             endif

          endif

       elseif( NDIFF_LOCATION == 'IN_LARGE_STEP2' ) then
          f_TEND    (:,:,:,:) = 0.D0
          f_TEND_pl (:,:,:,:) = 0.D0
          f_TENDq   (:,:,:,:) = 0.D0
          f_TENDq_pl(:,:,:,:) = 0.D0

          !------ numerical diffusion
          call numfilter_numerical_hdiff( rho,                       rho_pl,                       & !--- [IN]
                                          vx,                        vx_pl,                        & !--- [IN]
                                          vy,                        vy_pl,                        & !--- [IN]
                                          vz,                        vz_pl,                        & !--- [IN]
                                          w,                         w_pl,                         & !--- [IN]
                                          temd,                      temd_pl,                      & !--- [IN]
                                          q,                         q_pl,                         & !--- [IN]
                                          f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & !--- [INOUT]
                                          f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & !--- [INOUT]
                                          f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & !--- [INOUT]
                                          f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & !--- [INOUT]
                                          f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & !--- [INOUT]
                                          f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & !--- [INOUT]
                                          f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & !--- [INOUT]
                                          f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) !--- [INOUT]

          if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
             call numfilter_numerical_vdiff( rho,                       rho_pl,                       & !--- [IN]
                                             vx,                        vx_pl,                        & !--- [IN]
                                             vy,                        vy_pl,                        & !--- [IN]
                                             vz,                        vz_pl,                        & !--- [IN]
                                             w,                         w_pl,                         & !--- [IN]
                                             temd,                      temd_pl,                      & !--- [IN]
                                             q,                         q_pl,                         & !--- [IN]
                                             f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & !--- [INOUT]
                                             f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & !--- [INOUT]
                                             f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) !--- [INOUT]
          endif

          if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
             call numfilter_rayleigh_damping( rho,                       rho_pl,                       & !--- [IN]
                                              vx,                        vx_pl,                        & !--- [IN]
                                              vy,                        vy_pl,                        & !--- [IN]
                                              vz,                        vz_pl,                        & !--- [IN]
                                              w,                         w_pl,                         & !--- [IN]
                                              f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & !--- [INOUT]
                                              f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & !--- [INOUT]
                                              f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & !--- [INOUT]
                                              f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   )  ) !--- [INOUT]
          endif

       else !--- OUT_LARGE
          f_TEND    (:,:,:,:) = 0.D0
          f_TEND_pl (:,:,:,:) = 0.D0
          f_TENDq   (:,:,:,:) = 0.D0
          f_TENDq_pl(:,:,:,:) = 0.D0
       endif

       ! Smagorinksy-type SGS model [add] A.Noda 10/11/29
       if ( trim(TB_TYPE ) == 'SMG' ) then

          if ( ADM_prc_me /= ADM_prc_pl ) then ! for safety
             th_pl(:,:,:) = 0.D0
             vx_pl(:,:,:) = 0.D0
             vy_pl(:,:,:) = 0.D0
             vz_pl(:,:,:) = 0.D0
             w_pl (:,:,:) = 0.D0
          endif

          call tb_smg_driver( nl,                                                      &
                              rho,                       rho_pl,                       & !--- [IN] : density
                              PROG(:,:,:,I_RHOG  ),      PROG_pl(:,:,:,I_RHOG  ),      & !--- [IN]
                              PROGq(:,:,:,:),            PROGq_pl(:,:,:,:  ),          & !--- [IN]
                              vx,                        vx_pl,                        & !--- [IN] : Vx
                              vy,                        vy_pl,                        & !--- [IN] : Vy
                              vz,                        vz_pl,                        & !--- [IN] : Vz
                              w,                         w_pl,                         & !--- [IN] : w
                              tem,                       tem_pl,                       & !--- [IN] : temperature
                              q,                         q_pl,                         & !--- [IN] : q
                              th,                        th_pl,                        & !--- [IN] : potential temperature 
                              f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & !--- [INOUT]
                              f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & !--- [INOUT]
                              f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & !--- [INOUT]
                              f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & !--- [INOUT]
                              f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & !--- [INOUT]
                              f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & !--- [INOUT]
                              f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & !--- [INOUT]
                              f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) !--- [INOUT]
       endif

       !--- Nudging routines [add] Y.Niwa 08/09/09
       if ( FLAG_NUDGING ) then

          if ( nl == 1 ) then
             call ndg_update_var
          endif

          if ( nl == num_of_iteration_lstep ) then
             ndg_TEND_out = .true.
          else
             ndg_TEND_out = .false.
          endif

          call ndg_nudging_uvtp( rho,                      rho_pl,                      & !--- [IN]
                                 vx,                       vx_pl,                       & !--- [IN]
                                 vy,                       vy_pl,                       & !--- [IN]
                                 vz,                       vz_pl,                       & !--- [IN]
                                 w,                        w_pl,                        & !--- [IN]
                                 tem,                      tem_pl,                      & !--- [IN]
                                 pre,                      pre_pl,                      & !--- [IN]
                                 f_TEND(:,:,:,I_RHOGVX  ), f_TEND_pl(:,:,:,I_RHOGVX  ), & !--- [INOUT]
                                 f_TEND(:,:,:,I_RHOGVY  ), f_TEND_pl(:,:,:,I_RHOGVY  ), & !--- [INOUT]
                                 f_TEND(:,:,:,I_RHOGVZ  ), f_TEND_pl(:,:,:,I_RHOGVZ  ), & !--- [INOUT]
                                 f_TEND(:,:,:,I_RHOGW   ), f_TEND_pl(:,:,:,I_RHOGW   ), & !--- [INOUT]
                                 f_TEND(:,:,:,I_RHOGE   ), f_TEND_pl(:,:,:,I_RHOGE   ), & !--- [INOUT]
                                 f_TEND(:,:,:,I_RHOGETOT), f_TEND_pl(:,:,:,I_RHOGETOT), & !--- [INOUT]
                                 ndg_TEND_out                                           ) !--- [IN] ( TEND out )
       endif

       !--- sum the large step TEND ( advection + coriolis + num.diff.,SGS,nudge )
       g_TEND   (:,:,:,:) = g_TEND   (:,:,:,:) + f_TEND   (:,:,:,:)
       g_TEND_pl(:,:,:,:) = g_TEND_pl(:,:,:,:) + f_TEND_pl(:,:,:,:)

       !------------------------------------------------------------------------
       !> SMALL step
       !------------------------------------------------------------------------
       if ( nl /= 1 ) then ! update split values
          PROG_split   (:,:,:,:) = PROG0   (:,:,:,:) - PROG   (:,:,:,:)
          PROG_split_pl(:,:,:,:) = PROG0_pl(:,:,:,:) - PROG_pl(:,:,:,:)
       else
          PROG_split   (:,:,:,:) = 0.D0
          PROG_split_pl(:,:,:,:) = 0.D0
       endif

       !------ Core routine for small step
       !------    1. By this subroutine, prognostic variables 
       !------       ( rho,.., rhoge ) are calculated through
       !------       the 'num_of_iteration_sstep(nl)'-th times small step.
       !------    2. grho, grhogvx, ..., and  grhoge has the large step
       !------       tendencies initially, however,
       !------       they are re-used in this subroutine.
       !------
       if ( TIME_SPLIT ) then
          small_step_ite = num_of_iteration_sstep(nl)
          small_step_dt  = TIME_DTS
       else
          small_step_ite = 1
          small_step_dt  = TIME_DTL / (num_of_iteration_lstep+1-nl)
       endif

       call vi_small_step( PROG(:,:,:,I_RHOG  ),       PROG_pl(:,:,:,I_RHOG  ),       & !--- [INOUT] prognostic variables
                           PROG(:,:,:,I_RHOGVX),       PROG_pl(:,:,:,I_RHOGVX),       & !--- [INOUT]
                           PROG(:,:,:,I_RHOGVY),       PROG_pl(:,:,:,I_RHOGVY),       & !--- [INOUT]
                           PROG(:,:,:,I_RHOGVZ),       PROG_pl(:,:,:,I_RHOGVZ),       & !--- [INOUT]
                           PROG(:,:,:,I_RHOGW ),       PROG_pl(:,:,:,I_RHOGW ),       & !--- [INOUT]
                           PROG(:,:,:,I_RHOGE ),       PROG_pl(:,:,:,I_RHOGE ),       & !--- [INOUT]
                           vx,                         vx_pl,                         & !--- [IN] diagnostic value
                           vy,                         vy_pl,                         & !--- [IN]
                           vz,                         vz_pl,                         & !--- [IN]
                           eth,                        eth_pl,                        & !--- [IN]
                           rhogd,                      rhogd_pl,                      & !--- [IN]
                           pregd,                      pregd_pl,                      & !--- [IN]
                           g_TEND(:,:,:,I_RHOG    ),   g_TEND_pl(:,:,:,I_RHOG    ),   & !--- [IN] large step TEND
                           g_TEND(:,:,:,I_RHOGVX  ),   g_TEND_pl(:,:,:,I_RHOGVX  ),   & !--- [IN]
                           g_TEND(:,:,:,I_RHOGVY  ),   g_TEND_pl(:,:,:,I_RHOGVY  ),   & !--- [IN]
                           g_TEND(:,:,:,I_RHOGVZ  ),   g_TEND_pl(:,:,:,I_RHOGVZ  ),   & !--- [IN]
                           g_TEND(:,:,:,I_RHOGW   ),   g_TEND_pl(:,:,:,I_RHOGW   ),   & !--- [IN]
                           g_TEND(:,:,:,I_RHOGE   ),   g_TEND_pl(:,:,:,I_RHOGE   ),   & !--- [IN]
                           g_TEND(:,:,:,I_RHOGETOT),   g_TEND_pl(:,:,:,I_RHOGETOT),   & !--- [IN]
                           PROG_split(:,:,:,I_RHOG  ), PROG_split_pl(:,:,:,I_RHOG  ), & !--- [INOUT] split value
                           PROG_split(:,:,:,I_RHOGVX), PROG_split_pl(:,:,:,I_RHOGVX), & !--- [INOUT]
                           PROG_split(:,:,:,I_RHOGVY), PROG_split_pl(:,:,:,I_RHOGVY), & !--- [INOUT]
                           PROG_split(:,:,:,I_RHOGVZ), PROG_split_pl(:,:,:,I_RHOGVZ), & !--- [INOUT]
                           PROG_split(:,:,:,I_RHOGW ), PROG_split_pl(:,:,:,I_RHOGW ), & !--- [INOUT]
                           PROG_split(:,:,:,I_RHOGE ), PROG_split_pl(:,:,:,I_RHOGE ), & !--- [INOUT]
                           v_mean_c,                   v_mean_c_pl,                   & !--- [OUT] mean value
                           small_step_ite,                                            & !--- [IN]
                           small_step_dt                                              ) !--- [IN]

       !------------------------------------------------------------------------
       !>  Tracer advection
       !------------------------------------------------------------------------
       if ( TRC_ADV_TYPE == 'MIURA2004' ) then

          if ( nl == num_of_iteration_lstep ) then

             call src_update_tracer( TRC_VMAX,                                              & !--- [IN]
                                     PROGq(:,:,:,:),           PROGq_pl(:,:,:,:),           & !--- [INOUT]
                                     PROG0(:,:,:,I_RHOG),      PROG0_pl(:,:,:,I_RHOG),      & !--- [IN]
                                     v_mean_c(:,:,:,I_rhog),   v_mean_c_pl(:,:,:,I_rhog),   & !--- [IN]
                                     v_mean_c(:,:,:,I_rhogvx), v_mean_c_pl(:,:,:,I_rhogvx), & !--- [IN]
                                     v_mean_c(:,:,:,I_rhogvy), v_mean_c_pl(:,:,:,I_rhogvy), & !--- [IN]
                                     v_mean_c(:,:,:,I_rhogvz), v_mean_c_pl(:,:,:,I_rhogvz), & !--- [IN]
                                     v_mean_c(:,:,:,I_rhogw),  v_mean_c_pl(:,:,:,I_rhogw),  & !--- [IN]
                                     f_TEND (:,:,:,I_RHOG),    f_TEND_pl (:,:,:,I_RHOG),    & !--- [IN]
                                     TIME_DTL                                               ) !--- [IN]

             PROGq(:,:,:,:) = PROGq(:,:,:,:) + TIME_DTL * f_TENDq(:,:,:,:) ! update rhogq by viscosity

             PROGq(:,ADM_kmin-1,:,:) = 0.D0
             PROGq(:,ADM_kmax+1,:,:) = 0.D0

             if ( ADM_prc_pl == ADM_prc_me ) then
                PROGq_pl(:,:,:,:) = PROGq_pl(:,:,:,:) + TIME_DTL * f_TENDq_pl(:,:,:,:)

                PROGq_pl(:,ADM_kmin-1,:,:) = 0.D0
                PROGq_pl(:,ADM_kmax+1,:,:) = 0.D0
             endif

             !<---- I don't recommend adding the hyperviscosity term
             !<---- because of numerical instability in this case. ( H.Tomita )

          endif ! Last large step only

       elseif( TRC_ADV_TYPE == 'DEFAULT' ) then

          do nq = 1, TRC_VMAX

             call src_advection_convergence( v_mean_c(:,:,:,I_rhogvx), v_mean_c_pl(:,:,:,I_rhogvx), & !--- [IN]
                                             v_mean_c(:,:,:,I_rhogvy), v_mean_c_pl(:,:,:,I_rhogvy), & !--- [IN]
                                             v_mean_c(:,:,:,I_rhogvz), v_mean_c_pl(:,:,:,I_rhogvz), & !--- [IN]
                                             v_mean_c(:,:,:,I_rhogw),  v_mean_c_pl(:,:,:,I_rhogw),  & !--- [IN]
                                             q(:,:,:,nq),              q_pl(:,:,:,nq),              & !--- [IN]
                                             g_TENDq(:,:,:,nq),        g_TENDq_pl(:,:,:,nq),        & !--- [OUT]
                                             I_SRC_default                                          ) !--- [IN]  [mod] H.Yashiro 20120530

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

       endif

       !--- TKE fixer ( TKE >= 0.D0 )
       ! 2011/08/16 M.Satoh [comment] need this fixer for every small time steps
       if ( I_TKE >= 0 ) then
          if ( TRC_ADV_TYPE == 'DEFAULT' .OR. nl == num_of_iteration_lstep ) then
             do l  = 1, ADM_lall
                do k  = 1, ADM_kall
                do ij = 1, ADM_gall
                   TKEg_corr = TKE_MIN * VMTR_GSGAM2(ij,k,l) - PROGq(ij,k,l,I_TKE)

                   if ( TKEg_corr >= 0.D0 ) then
                      PROG(ij,k,l,I_RHOGE)       = PROG(ij,k,l,I_RHOGE)       - TKEg_corr
                      PROGq(ij,k,l,I_TKE) = PROGq(ij,k,l,I_TKE) + TKEg_corr
                   endif
                enddo
                enddo
             enddo

             if ( ADM_prc_pl == ADM_prc_me ) then
                do l  = 1, ADM_lall_pl
                   do k  = 1, ADM_kall
                   do ij = 1, ADM_gall_pl
                      TKEg_corr = TKE_MIN * VMTR_GSGAM2_pl(ij,k,l) - PROGq_pl(ij,k,l,I_TKE)

                      if ( TKEg_corr >= 0.D0 ) then
                         PROG_pl(ij,k,l,I_RHOGE)       = PROG_pl(ij,k,l,I_RHOGE)       - TKEg_corr
                         PROGq_pl(ij,k,l,I_TKE) = PROGq_pl(ij,k,l,I_TKE) + TKEg_corr
                      endif
                   enddo
                   enddo
                enddo
             endif

          endif
       endif

       !------ Update
       if ( nl /= num_of_iteration_lstep ) then
          ! communication
          call COMM_data_transfer( PROG, PROG_pl )

          PROG(suf(ADM_gall_1d,1),:,:,:) = PROG(suf(ADM_gmax+1,ADM_gmin),:,:,:)
          PROG(suf(1,ADM_gall_1d),:,:,:) = PROG(suf(ADM_gmin,ADM_gmax+1),:,:,:)
       endif

    enddo !--- large step 

    ! 2011/08/16 M.Satoh move from above =>
    ! 2010/05/06 M.Satoh (comment on TDK)
    ! 2011/08/16 M.Satoh bug fix for TDK: conv => TEND
    ! 2011/08/16 M.Satoh
    !    use rhogq-rhogq_old to calculate qconv & qTEND
    ! Kuo param.    : total convergence of water vapor [integ over time]
    ! Tiedtke scheme: TEND of water vapor [per time]
    if ( CP_TYPE == 'KUO' ) then
       diagvar(:,:,:,I_RHOGQV_CONV) = PROGq(:,:,:,I_QV) - PROGq0(:,:,:,I_QV)
       if ( ADM_prc_pl == ADM_prc_me ) then
          diagvar_pl(:,:,:,I_RHOGQV_CONV) = PROGq_pl(:,:,:,I_QV) - PROGq0_pl(:,:,:,I_QV)
       endif
    endif

    if ( CP_TYPE == 'TDK' ) then
       ! 2011.08.16 M.Satoh,bug fix: conv => TEND
       ! qv_dyn_tend = v grad q
       !             = ( div(rho v q) - div(rho v)*q )/rho
       diagvar(:,:,:,I_QV_DYN_TEND) = ( ( PROGq(:,:,:,I_QV)-PROGq0(:,:,:,I_QV) ) &
                                      - ( PROG(:,:,:,I_RHOG)-PROG0(:,:,:,I_RHOG) )             &
                                        * PROGq(:,:,:,I_QV) / PROG(:,:,:,I_RHOG)           & ! = q(:,:,:,I_QV)
                                      ) / TIME_DTL / PROG(:,:,:,I_RHOG)
       if ( ADM_prc_pl == ADM_prc_me ) then
          diagvar_pl(:,:,:,I_QV_DYN_TEND) = ( ( PROGq_pl(:,:,:,I_QV)-PROGq0_pl(:,:,:,I_QV) ) &
                                            - ( PROG_pl(:,:,:,I_RHOG)-PROG0_pl(:,:,:,I_RHOG) )             &
                                              * PROGq_pl(:,:,:,I_QV) / PROG_pl(:,:,:,I_RHOG)           & ! = q_pl(:,:,:,I_QV)
                                            ) / TIME_DTL / PROG_pl(:,:,:,I_RHOG)
       endif
    endif
    ! <= 2011/08/16 M.Satoh move from above

    call prgvar_set( PROG(:,:,:,I_RHOG),   PROG_pl(:,:,:,I_RHOG),   & !--- [IN]
                     PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & !--- [IN]
                     PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & !--- [IN]
                     PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ), & !--- [IN]
                     PROG(:,:,:,I_RHOGW),  PROG_pl(:,:,:,I_RHOGW),  & !--- [IN]
                     PROG(:,:,:,I_RHOGE),  PROG_pl(:,:,:,I_RHOGE),  & !--- [IN]
                     PROGq(:,:,:,:),       PROGq_pl(:,:,:,:),       & !--- [IN]
                     0                                              ) !--- [IN]

#ifdef PAPI_OPS
    ! <-- [add] PAPI R.Yoshida 20121022
    !call PAPIF_flips( PAPI_real_time_i, PAPI_proc_time_i, PAPI_flpins, PAPI_mflins, PAPI_check )
    call PAPIF_flops( PAPI_real_time_o, PAPI_proc_time_o, PAPI_flpops, PAPI_mflops, PAPI_check )
#endif

    return
  end subroutine dynstep

  !-----------------------------------------------------------------------------
  subroutine dynstep2
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
    use mod_cnst, only: &
       CNST_RAIR, &
       CNST_RVAP, &
       CNST_CV
    use mod_time, only: &
       TIME_DTL
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    use mod_oprt, only : &
       OPRT_horizontalize_vec
    use mod_vmtr, only: &
       VMTR_GSGAM2,     &
       VMTR_GSGAM2_pl,  &
       VMTR_GSGAM2H,    &
       VMTR_GSGAM2H_pl, &
       VMTR_GZXH,       &
       VMTR_GZXH_pl,    &
       VMTR_GZYH,       &
       VMTR_GZYH_pl,    &
       VMTR_GZZH,       &
       VMTR_GZZH_pl,    &
       VMTR_PHI,        &
       VMTR_PHI_pl
    use mod_runconf, only: &
       TRC_VMAX,       &
       I_QV,           &
       NQW_STR,        &
       NQW_END,        &
       CVW,            &
       NDIFF_LOCATION
    use mod_bsstate, only: &
       tem_bs, tem_bs_pl
    use mod_bndcnd, only: &
       bndcnd_all
    use mod_prgvar, only: &
       prgvar_set,&
       prgvar_get
    use mod_numfilter, only: &
       NUMFILTER_DOrayleigh,       & ! [add] H.Yashiro 20120530
       NUMFILTER_DOverticaldiff,   & ! [add] H.Yashiro 20120530
       numfilter_rayleigh_damping, &
       numfilter_numerical_hdiff,  &
       numfilter_numerical_vdiff
    implicit none

    !--- forcing tendensy of rhog  ( G^{1/2} X gamma2 )
    real(8) :: frhog   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhogvx  ( G^{1/2} X gamma2 )
    real(8) :: frhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhogvy  ( G^{1/2} X gamma2 )
    real(8) :: frhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhogvz  ( G^{1/2} X gamma2 )
    real(8) :: frhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhogw  ( G^{1/2} X gamma2 )
    real(8) :: frhogw(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhoge  ( G^{1/2} X gamma2 )
    real(8) :: frhoge(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhoge_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- tendensy of rhogetot  ( G^{1/2} X gamma2 )
    real(8) :: frhogetot(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogetot_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhogq  ( G^{1/2} X gamma2 )
    real(8) :: frhogq(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: frhogq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !
    !--- rho X ( G^{1/2} X gamma2 )
    real(8) :: rhog(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- rho X ( G^{1/2} X gamma2 ) X vx
    real(8) :: rhogvx(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- rho X ( G^{1/2} X gamma2 ) X vy
    real(8) :: rhogvy(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- rho X ( G^{1/2} X gamma2 ) X vz
    real(8) :: rhogvz(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- rho X ( G^{1/2} X gamma2 ) X w
    real(8) :: rhogw(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- rho X ( G^{1/2} X gamma2 ) X ein
    real(8) :: rhoge(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhoge_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- rho X ( G^{1/2} X gamma2 ) X q
    real(8) :: rhogq(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: rhogq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    !
    !--- density ( physical )
    real(8) :: rho(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rho_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- horizontal velocity_x  ( physical )
    real(8) :: vx(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- horizontal velocity_y  ( physical )
    real(8) :: vy(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- horizontal velocity_z  ( physical )
    real(8) :: vz(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- vertical velocity ( physical )
    real(8) :: w(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: w_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- [IN]ternal energy  ( physical )
    real(8) :: ein(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: ein_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- mass concentration of water substance ( physical )
    real(8) :: q(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: q_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !--- pressure ( physical )
    real(8) :: pre(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: pre_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- temperature ( physical )
    real(8) :: tem(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: tem_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- temperature deviation from the base state ( physical )
    real(8) :: temd(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: temd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- temporary variables
    real(8) :: qd      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: qd_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rrhog   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rrhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: cv      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: cv_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: ij, k ,l, nq
    !---------------------------------------------------------------------------

!    if ( NDIFF_LOCATION == 'OUT_LARGE_STEP' ) then
!
!       !--- get from prg0
!       call prgvar_get( rhog,   rhog_pl,   & !--- [OUT]
!                        rhogvx, rhogvx_pl, & !--- [OUT]
!                        rhogvy, rhogvy_pl, & !--- [OUT]
!                        rhogvz, rhogvz_pl, & !--- [OUT]
!                        rhogw,  rhogw_pl,  & !--- [OUT]
!                        rhoge,  rhoge_pl,  & !--- [OUT]
!                        rhogq,  rhogq_pl,  & !--- [OUT]
!                        0                  ) !--- [IN]
!
!       !---< Generate diagnostic values and set the boudary conditions
!       do l = 1, ADM_lall
!          do k = 1, ADM_kall
!             do ij = 1, ADM_gall
!                rrhog(ij,k,l) = 1.D0 / rhog(ij,k,l)
!             enddo
!             do ij = 1, ADM_gall
!               rho(ij,k,l) = rhog(ij,k,l) / VMTR_GSGAM2(ij,k,l)
!             enddo
!             do ij = 1, ADM_gall
!                ein(ij,k,l) = rhoge(ij,k,l) * rrhog(ij,k,l)
!             enddo
!             do ij = 1, ADM_gall
!                vx(ij,k,l) = rhogvx(ij,k,l) * rrhog(ij,k,l)
!             enddo
!             do ij = 1, ADM_gall
!                vy(ij,k,l) = rhogvy(ij,k,l) * rrhog(ij,k,l)
!             enddo
!             do ij = 1, ADM_gall
!                vz(ij,k,l) = rhogvz(ij,k,l) * rrhog(ij,k,l)
!             enddo
!          enddo
!          do nq = 1, TRC_VMAX
!             do k = 1, ADM_kall
!                do ij = 1, ADM_gall
!                   q(ij,k,l,nq) = rhogq(ij,k,l,nq) * rrhog(ij,k,l)
!                enddo
!             enddo
!          enddo
!          do k = 1, ADM_kall
!             do ij = 1, ADM_gall
!                qd(ij,k,l) = 1.D0
!             enddo
!          enddo
!          do nq = NQW_STR, NQW_END
!             do k  = 1, ADM_kall
!             do ij = 1, ADM_gall
!                qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
!             enddo
!             enddo
!          enddo
!          do k  = 1, ADM_kall
!             do ij = 1, ADM_gall
!                cv(ij,k,l) = qd(ij,k,l) * CNST_CV
!             enddo
!          enddo
!          do nq = NQW_STR, NQW_END
!             do k  = 1, ADM_kall
!             do ij = 1, ADM_gall
!                cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
!             enddo
!             enddo
!          enddo
!          do k = 1, ADM_kall
!             do ij = 1, ADM_gall
!                tem(ij,k,l) = ein(ij,k,l) / cv(ij,k,l)
!             enddo
!             do ij = 1, ADM_gall
!                pre(ij,k,l) = rho(ij,k,l) * tem(ij,k,l) * ( qd(ij,k,l)*CNST_RAIR + q(ij,k,l,I_QV)*CNST_RVAP )
!             enddo
!          enddo
!          do k = ADM_kmin+1, ADM_kmax
!             do ij = 1, ADM_gall
!                w(ij,k,l) = rhogw(ij,k,l) &
!                          / ( 0.5D0 * ( GRD_afac(k) * rho(ij,k  ,l) &
!                                      + GRD_bfac(k) * rho(ij,k-1,l) ) * VMTR_GSGAM2H(ij,k,l) )
!             enddo
!          enddo
!
!          !--- boundary conditions
!          call bndcnd_all( ADM_gall,            & !--- [IN]
!                           vx    (:,:,l),       & !--- [INOUT]
!                           vy    (:,:,l),       & !--- [INOUT]
!                           vz    (:,:,l),       & !--- [INOUT]
!                           w     (:,:,l),       & !--- [INOUT]
!                           tem   (:,:,l),       & !--- [INOUT]
!                           rho   (:,:,l),       & !--- [INOUT]
!                           pre   (:,:,l),       & !--- [INOUT]
!                           ein   (:,:,l),       & !--- [INOUT]
!                           rhog  (:,:,l),       & !--- [INOUT]
!                           rhogvx(:,:,l),       & !--- [INOUT]
!                           rhogvy(:,:,l),       & !--- [INOUT]
!                           rhogvz(:,:,l),       & !--- [INOUT]
!                           rhogw (:,:,l),       & !--- [INOUT]
!                           rhoge (:,:,l),       & !--- [INOUT]
!                           VMTR_GSGAM2 (:,:,l),   & !--- [IN]
!                           VMTR_PHI    (:,:,l),   & !--- [IN]
!                           VMTR_C2Wfact(:,:,:,l)  ) !--- [IN]
!
!             !--- perturbations ( pred, rhod, temd )
!          do k = 1, ADM_kall
!             do ij = 1, ADM_gall
!                temd(ij,k,l) = tem(ij,k,l) - tem_bs(ij,k,l)
!             enddo
!          enddo
!
!       enddo ! region LOOP
!
!       if ( ADM_prc_me == ADM_prc_pl ) then
!
!          do l = 1, ADM_lall_pl
!             do k = 1, ADM_kall
!                do ij = 1, ADM_gall_pl
!                   rrhog_pl(ij,k,l) = 1.D0 / rhog_pl(ij,k,l)
!                enddo
!                do ij = 1, ADM_gall_pl
!                  rho_pl(ij,k,l)  = rhog_pl(ij,k,l) / VMTR_GSGAM2_pl(ij,k,l)
!                enddo
!                do ij = 1, ADM_gall_pl
!                   ein_pl(ij,k,l) = rhoge_pl(ij,k,l)  * rrhog_pl(ij,k,l)
!                enddo
!                do ij = 1, ADM_gall_pl
!                   vx_pl(ij,k,l)  = rhogvx_pl(ij,k,l) * rrhog_pl(ij,k,l)
!                enddo
!                do ij = 1, ADM_gall_pl
!                   vy_pl(ij,k,l)  = rhogvy_pl(ij,k,l) * rrhog_pl(ij,k,l)
!                enddo
!                do ij = 1, ADM_gall_pl
!                   vz_pl(ij,k,l)  = rhogvz_pl(ij,k,l) * rrhog_pl(ij,k,l)
!                enddo
!             enddo
!             do nq = 1, TRC_VMAX
!                do k = 1, ADM_kall
!                   do ij = 1, ADM_gall_pl
!                      q_pl(ij,k,l,nq) = rhogq_pl(ij,k,l,nq) * rrhog_pl(ij,k,l)
!                   enddo
!                enddo
!             enddo
!             do k = 1, ADM_kall
!                do ij = 1, ADM_gall_pl
!                   qd_pl(ij,k,l) = 1.D0
!                enddo
!             enddo
!             do nq = NQW_STR, NQW_END
!                do k = 1, ADM_kall
!                   do ij = 1, ADM_gall_pl
!                      qd_pl(ij,k,l) = qd_pl(ij,k,l) - q_pl(ij,k,l,nq)
!                   enddo
!                enddo
!             enddo
!             do k = 1, ADM_kall
!                do ij = 1, ADM_gall_pl
!                   cv_pl(ij,k,l) = qd_pl(ij,k,l) * CNST_CV
!                enddo
!             enddo
!             do nq = NQW_STR, NQW_END
!                do k = 1, ADM_kall
!                   do ij = 1, ADM_gall_pl
!                      cv_pl(ij,k,l) = cv_pl(ij,k,l) + q_pl(ij,k,l,nq) * CVW(nq)
!                   enddo
!                enddo
!             enddo
!             do k = 1, ADM_kall
!                do ij = 1, ADM_gall_pl
!                   tem_pl(ij,k,l) = ein_pl(ij,k,l) / cv_pl(ij,k,l)
!                enddo
!                do ij = 1, ADM_gall_pl
!                   pre_pl(ij,k,l) = rho_pl(ij,k,l) * tem_pl(ij,k,l) &
!                        * ( qd_pl(ij,k,l)*CNST_RAIR+q_pl(ij,k,l,I_QV)*CNST_RVAP )
!                enddo
!             enddo
!             do k = ADM_kmin+1, ADM_kmax
!                do ij = 1, ADM_gall_pl
!                   w_pl(ij,k,l) = rhogw_pl(ij,k,l) &
!                                / ( 0.5D0 * ( GRD_afac(k) * rho_pl(ij,k  ,l) &
!                                            + GRD_bfac(k) * rho_pl(ij,k-1,l) ) * VMTR_GSGAM2H_pl(ij,k,l) )
!                enddo
!             enddo
!
!             call bndcnd_all( ADM_gall_pl,            & !--- [IN]
!                              vx_pl    (:,:,l),       & !--- [INOUT]
!                              vy_pl    (:,:,l),       & !--- [INOUT]
!                              vz_pl    (:,:,l),       & !--- [INOUT]
!                              w_pl     (:,:,l),       & !--- [INOUT]
!                              tem_pl   (:,:,l),       & !--- [INOUT]
!                              rho_pl   (:,:,l),       & !--- [INOUT]
!                              pre_pl   (:,:,l),       & !--- [INOUT]
!                              ein_pl   (:,:,l),       & !--- [INOUT]
!                              rhog_pl  (:,:,l),       & !--- [INOUT]
!                              rhogvx_pl(:,:,l),       & !--- [INOUT]
!                              rhogvy_pl(:,:,l),       & !--- [INOUT]
!                              rhogvz_pl(:,:,l),       & !--- [INOUT]
!                              rhogw_pl (:,:,l),       & !--- [INOUT]
!                              rhoge_pl (:,:,l),       & !--- [INOUT]
!                              VMTR_GSGAM2_pl (:,:,l),   & !--- [IN]
!                              VMTR_PHI_pl    (:,:,l),   & !--- [IN]
!                              VMTR_C2Wfact_pl(:,:,:,l)  ) !--- [IN]
!
!             do k = 1, ADM_kall
!                do ij = 1, ADM_gall_pl
!                   temd_pl(ij,k,l) = tem_pl(ij,k,l) - tem_bs_pl(ij,k,l)
!                enddo
!             enddo
!
!          enddo
!       endif
!
!             do l = 1, ADM_lall
!                do k = 1, ADM_kall
!                   do ij = 1, ADM_gall
!                      frhog  (ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall
!                      frhogvx(ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall
!                      frhogvy(ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall
!                      frhogvz(ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall
!                      frhogw (ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall
!                      frhoge (ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall
!                      frhogetot(ij,k,l) = 0.D0
!                   enddo
!                enddo
!             enddo
!             do nq = 1,TRC_VMAX
!             do l = 1, ADM_lall
!                do k = 1, ADM_kall
!                   do ij = 1, ADM_gall
!                      frhogq(ij,k,l,nq) = 0.D0
!                   enddo
!                enddo
!             enddo
!             enddo
!
!             do l = 1, ADM_lall_pl
!                do k = 1, ADM_kall
!                   do ij = 1, ADM_gall_pl
!                      frhog_pl  (ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall_pl
!                      frhogvx_pl(ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall_pl
!                      frhogvy_pl(ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall_pl
!                      frhogvz_pl(ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall_pl
!                      frhogw_pl (ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall_pl
!                      frhoge_pl (ij,k,l) = 0.D0
!                   enddo
!                   do ij = 1, ADM_gall_pl
!                      frhogetot_pl(ij,k,l) = 0.D0
!                   enddo
!                enddo
!             enddo
!             do nq = 1, TRC_VMAX
!             do l  = 1, ADM_lall_pl
!                do k = 1, ADM_kall
!                   do ij = 1, ADM_gall_pl
!                      frhogq_pl(ij,k,l,nq) = 0.D0
!                   enddo
!                enddo
!             enddo
!             enddo
!
!             if ( NUMFILTER_DOrayleigh ) then ! [add] H.Yashiro 20120530
!
!             !------ rayleigh damping
!             call numfilter_rayleigh_damping( &
!                  rho, rho_pl,                & !--- [IN]
!                  vx, vx_pl,                  & !--- [IN]
!                  vy, vy_pl,                  & !--- [IN]
!                  vz, vz_pl,                  & !--- [IN]
!                  w, w_pl,                    & !--- [IN]
!                  frhogvx, frhogvx_pl,        & !--- [INOUT]
!                  frhogvy, frhogvy_pl,        & !--- [INOUT]
!                  frhogvz, frhogvz_pl,        & !--- [INOUT]
!                  frhogw, frhogw_pl           ) !--- [INOUT]
!
!             endif
!
!             !------ numerical diffusion
!             call numfilter_numerical_hdiff(  &
!                  rho, rho_pl,                & !--- [IN]
!                  vx, vx_pl,                  & !--- [IN]
!                  vy, vy_pl,                  & !--- [IN]
!                  vz, vz_pl,                  & !--- [IN]
!                  w , w_pl,                   & !--- [IN]
!                  temd, temd_pl,              & !--- [IN]
!                  q , q_pl,                   & !--- [IN]
!                  frhog, frhog_pl,            & !--- [INOUT]
!                  frhogvx, frhogvx_pl,        & !--- [INOUT]
!                  frhogvy, frhogvy_pl,        & !--- [INOUT]
!                  frhogvz, frhogvz_pl,        & !--- [INOUT]
!                  frhogw , frhogw_pl,         & !--- [INOUT]
!                  frhoge , frhoge_pl,         & !--- [INOUT]
!                  frhogetot , frhogetot_pl,         & !--- [INOUT]
!                  frhogq , frhogq_pl )          !--- [INOUT]
!
!             if ( NUMFILTER_DOverticaldiff ) then ! [add] H.Yashiro 20120530
!
!             call numfilter_numerical_vdiff(  &
!                  rho,   rho_pl,              &  !--- [IN]
!                  vx, vx_pl,                  &  !--- [IN]
!                  vy, vy_pl,                  &  !--- [IN]
!                  vz, vz_pl,                  &  !--- [IN]
!                  w , w_pl,                   &  !--- [IN]
!                  temd,temd_pl,               &  !--- [IN]
!                  q,q_pl,                     &  !--- [IN]
!                  frhog, frhog_pl,            &  !--- [INOUT]
!                  frhogvx, frhogvx_pl,        &  !--- [INOUT]
!                  frhogvy, frhogvy_pl,        &  !--- [INOUT]
!                  frhogvz, frhogvz_pl,        &  !--- [INOUT]
!                  frhogw , frhogw_pl,         &  !--- [INOUT]
!                  frhoge , frhoge_pl,         &  !--- [INOUT]
!                  frhogetot , frhogetot_pl,         & !--- [INOUT]
!                  frhogq , frhogq_pl          )  !--- [INOUT]
!
!             endif
!
!       call OPRT_horizontalize_vec( frhogvx, frhogvx_pl, & !--- [INOUT]
!                                    frhogvy, frhogvy_pl, & !--- [INOUT]
!                                    frhogvz, frhogvz_pl  ) !--- [INOUT]
!
!       rhog   = rhog   + TIME_DTL * frhog
!       rhogvx = rhogvx + TIME_DTL * frhogvx
!       rhogvy = rhogvy + TIME_DTL * frhogvy
!       rhogvz = rhogvz + TIME_DTL * frhogvz
!       rhogw  = rhogw  + TIME_DTL * frhogw
!       rhoge  = rhoge  + TIME_DTL * frhoge
!       rhogq  = rhogq  + TIME_DTL * frhogq
!
!       if ( ADM_prc_me == ADM_prc_pl ) then
!          rhog_pl   = rhog_pl   + TIME_DTL * frhog_pl
!          rhogvx_pl = rhogvx_pl + TIME_DTL * frhogvx_pl
!          rhogvy_pl = rhogvy_pl + TIME_DTL * frhogvy_pl
!          rhogvz_pl = rhogvz_pl + TIME_DTL * frhogvz_pl
!          rhogw_pl  = rhogw_pl  + TIME_DTL * frhogw_pl
!          rhoge_pl  = rhoge_pl  + TIME_DTL * frhoge_pl
!          rhogq_pl  = rhogq_pl  + TIME_DTL * frhogq_pl
!       endif
!
!       call prgvar_set( rhog,   rhog_pl,   & !--- [IN]
!                        rhogvx, rhogvx_pl, & !--- [IN]
!                        rhogvy, rhogvy_pl, & !--- [IN]
!                        rhogvz, rhogvz_pl, & !--- [IN]
!                        rhogw,  rhogw_pl,  & !--- [IN]
!                        rhoge,  rhoge_pl,  & !--- [IN]
!                        rhogq,  rhogq_pl,  & !--- [IN]
!                        0                  )
!
!    endif

    return
  end subroutine dynstep2

end module mod_dynstep
