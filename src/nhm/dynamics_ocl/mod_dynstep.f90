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
  !                10-08-20   A.T.Noda: Bug fix (Qconv should be tendency, and not be multiplied by DT)
  !                10-11-29   A.T.Noda: Introduce the Smagorinsky model
  !                11-08-16   M.Satoh: bug fix for TDK: conv => tendency
  !                           qv_dyn_tend = v grad q
  !                                       = ( div(rho v q) - div(rho v)*q )/rho
  !                11-08-16   M.Satoh: move codes related to CP_TYPE below the tracer calculation
  !                11-11-28   Y.Yamada: Merge Terai-san timer into the original code.
  !                12-03-09   S.Iga: tuned (phase4-1)
  !                12-04-06   T.yamaura: optimized for K
  !                12-05-30   T.Yashiro: Change arguments from character to index/switch
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
  !-----------------------------------------------------------------------------

contains 

  !-----------------------------------------------------------------------------
  subroutine dynstep
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
       TIME_INTEG_TYPE, &
       TIME_SSTEP_MAX,  &
       TIME_DTL,        &
       TIME_DTS,        &
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
       VMTR_GZZH_pl
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
       rho_bs, rho_bs_pl, &
       phi, phi_pl
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
       src_update_tracer_rev,       & ! Y.Niwa add 080124
       I_SRC_default                  ! [add] H.Yashiro 20120530
    use mod_ndg, only: & ! Y.Niwa add 08/09/09
       ndg_nudging_uvtp, &
       ndg_update_var
    use mod_tb_smg, only: & ! [add] 10/11/29 A.Noda
       tb_smg_driver
    implicit none
    
    !--- forcing tendensy of rhog  ( G^{1/2} X gamma2 )
    real(8) :: frhog   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- forcing tendensy of rhogvx  ( G^{1/2} X gamma2 )
    real(8) :: frhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- forcing tendensy of rhogvy  ( G^{1/2} X gamma2 )
    real(8) :: frhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- forcing tendensy of rhogvz  ( G^{1/2} X gamma2 )
    real(8) :: frhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- forcing tendensy of rhogw  ( G^{1/2} X gamma2 )
    real(8) :: frhogw   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- forcing tendensy of rhoge  ( G^{1/2} X gamma2 )
    real(8) :: frhoge   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhoge_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendensy of rhogetot  ( G^{1/2} X gamma2 )
    real(8) :: frhogetot   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogetot_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- forcing tendensy of rhogq  ( G^{1/2} X gamma2 )
    real(8) :: frhogq   (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)
    real(8) :: frhogq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !--- rho X ( G^{1/2} X gamma2 )
    real(8) :: rhog   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- rho X ( G^{1/2} X gamma2 ) X vx
    real(8) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- rho X ( G^{1/2} X gamma2 ) X vy
    real(8) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- rho X ( G^{1/2} X gamma2 ) X vz
    real(8) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- rho X ( G^{1/2} X gamma2 ) X w
    real(8) :: rhogw   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- rho X ( G^{1/2} X gamma2 ) X ein
    real(8) :: rhoge   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhoge_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- rho X ( G^{1/2} X gamma2 ) X q
    real(8) :: rhogq   (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)
    real(8) :: rhogq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !--- old rho X ( G^{1/2} X gamma2 )
    real(8) :: rhog_old   (ADM_gall,   ADM_kall,ADM_lall   )          
    real(8) :: rhog_old_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) 

    !--- rho X ( G^{1/2} X gamma2 ) X q_old
    real(8) :: rhogq_old   (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)
    real(8) :: rhogq_old_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

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

    !--- tendensy of rhog  ( G^{1/2} X gamma2 )
    real(8) :: grhog   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendensy of rhogvx  ( G^{1/2} X gamma2 )
    real(8) :: grhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendensy of rhogvy  ( G^{1/2} X gamma2 )
    real(8) :: grhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendensy of rhogvz  ( G^{1/2} X gamma2 )
    real(8) :: grhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendensy of rhogw  ( G^{1/2} X gamma2 )
    real(8) :: grhogw   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendensy of rhoge  ( G^{1/2} X gamma2 )
    real(8) :: grhoge   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhoge_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendensy of rhogetot  ( G^{1/2} X gamma2 )
    real(8) :: grhogetot   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogetot_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendensy of rhogq  ( G^{1/2} X gamma2 )
    real(8) :: grhogq   (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)
    real(8) :: grhogq_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !--- deviation of rhog from the current large step ( G^{1/2} X gamma2 )
    real(8) :: rhog_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- deviation of rhogvx from the current large step ( G^{1/2} X gamma2 )
    real(8) :: rhogvx_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvx_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- deviation of rhogvy from the current large step ( G^{1/2} X gamma2 )
    real(8) :: rhogvy_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvy_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- deviation of rhogvz from the current large step ( G^{1/2} X gamma2 )
    real(8) :: rhogvz_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvz_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- deviation of rhogw from the current large step ( G^{1/2} X gamma2 )
    real(8) :: rhogw_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogw_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- deviation of rhoge from the current large step ( G^{1/2} X gamma2 )
    real(8) :: rhoge_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhoge_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer, parameter :: vmax_mean_c = 5

    integer, parameter :: I_rhogvx = 1 ! mean rhogvx for tracer advection  ( G^{1/2} X gamma2 )
    integer, parameter :: I_rhogvy = 2 ! mean rhogvy for tracer advection  ( G^{1/2} X gamma2 )
    integer, parameter :: I_rhogvz = 3 ! mean rhogvz for tracer advection  ( G^{1/2} X gamma2 )
    integer, parameter :: I_rhog   = 4 ! mean rhog for tracer advection  ( G^{1/2} X gamma2 )
    integer, parameter :: I_rhogw  = 5 ! mean rhogw for tracer advection  ( G^{1/2} X gamma2 )

    real(8) :: v_mean_c   (ADM_gall,   ADM_kall,ADM_lall   ,vmax_mean_c)
    real(8) :: v_mean_c_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax_mean_c)

    !--- temporary variables
    real(8) :: qd      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: qd_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rrhog   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rrhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: cv      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: cv_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), parameter :: TKE_MIN = 0.D0
    real(8)            :: TKEg_corr

    integer :: small_step_ite
    real(8) :: small_step_dt

    logical, save :: iflag = .true.
    integer, save :: num_of_iteration_lstep    ! number of large steps ( 2-4 )
    integer, save :: num_of_iteration_sstep(4) ! number of small steps in each of large steps

    integer :: ij, k ,l, nq, nl
    !---------------------------------------------------------------------------

    if ( iflag ) then
       iflag = .false.

       select case(trim(TIME_INTEG_TYPE))
       case('RK2')
          num_of_iteration_lstep = 2
          num_of_iteration_sstep(1) = TIME_SSTEP_MAX/2
          num_of_iteration_sstep(2) = TIME_SSTEP_MAX
       case('RK3')
          num_of_iteration_lstep = 3
          num_of_iteration_sstep(1) = TIME_SSTEP_MAX/3
          num_of_iteration_sstep(2) = TIME_SSTEP_MAX/2
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

#ifdef _FJTIMER_
    call timer_sta(2000)
    call timer_sta(2001)
#endif
    !--- get from prg0
    call prgvar_get( rhog,   rhog_pl,   & !--- [OUT]
                     rhogvx, rhogvx_pl, & !--- [OUT]
                     rhogvy, rhogvy_pl, & !--- [OUT]
                     rhogvz, rhogvz_pl, & !--- [OUT]
                     rhogw,  rhogw_pl,  & !--- [OUT]
                     rhoge,  rhoge_pl,  & !--- [OUT]
                     rhogq,  rhogq_pl,  & !--- [OUT]
                     0                  ) !--- [IN]

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
          do ij = 1, ADM_gall
             rhog_old(ij,k,l) = rhog(ij,k,l)
          enddo
       enddo
    enddo
!OCL SERIAL
    do nq = 1, TRC_VMAX
    do l  = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
          do ij = 1, ADM_gall
             rhogq_old(ij,k,l,nq) = rhogq(ij,k,l,nq)
          enddo
       enddo
    enddo
    enddo
!OCL SERIAL
    do l = 1, ADM_lall_pl
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
          do ij = 1, ADM_gall_pl
             rhog_old_pl(ij,k,l) = rhog_pl(ij,k,l)
          enddo
       enddo
    enddo
!OCL SERIAL
    do nq = 1,TRC_VMAX
    do l  = 1, ADM_lall_pl
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
          do ij = 1, ADM_gall_pl
             rhogq_old_pl(ij,k,l,nq) = rhogq_pl(ij,k,l,nq)
          enddo
       enddo
    enddo
    enddo

    !--- initialization of split values
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhog_split  (ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhogvx_split(ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhogvy_split(ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhogvz_split(ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhogw_split (ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhoge_split (ij,k,l) = 0.D0
          enddo
       enddo
    enddo

!OCL SERIAL
    do l = 1, ADM_lall_pl
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
          do ij = 1, ADM_gall_pl
             rhog_split_pl  (ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
          do ij = 1, ADM_gall_pl
             rhogvx_split_pl(ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
          do ij = 1, ADM_gall_pl
             rhogvy_split_pl(ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
          do ij = 1, ADM_gall_pl
             rhogvz_split_pl(ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
          do ij = 1, ADM_gall_pl
             rhogw_split_pl (ij,k,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
          do ij = 1, ADM_gall_pl
             rhoge_split_pl (ij,k,l) = 0.D0
          enddo
       enddo
    enddo

#ifdef _FJTIMER_
    call timer_end(2001)
    call timer_sta(2002)
#endif

    !---------------------------------------------------------------------------
    !
    !> Start large time step integration
    !
    !---------------------------------------------------------------------------
    do nl = 1, num_of_iteration_lstep

       !------------------------------------------------------------------------
       !> LARGE step
       !------------------------------------------------------------------------
#ifdef _FJTIMER_
       call timer_sta(2003)
#endif
       if ( nl /= 1 ) then
          !------ get from prg1
          call prgvar_get( rhog,   rhog_pl,   & !--- [OUT]
                           rhogvx, rhogvx_pl, & !--- [OUT]
                           rhogvy, rhogvy_pl, & !--- [OUT]
                           rhogvz, rhogvz_pl, & !--- [OUT]
                           rhogw,  rhogw_pl,  & !--- [OUT]
                           rhoge,  rhoge_pl,  & !--- [OUT]
                           rhogq,  rhogq_pl,  & !--- [OUT]
                           1                  ) !--- [IN]
       endif

       !---< Generate diagnostic values and set the boudary conditions
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rrhog(ij,k,l) = 1.D0 / rhog(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
               rho(ij,k,l) = rhog(ij,k,l) / VMTR_GSGAM2(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                ein(ij,k,l) = rhoge(ij,k,l) * rrhog(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                vx(ij,k,l) = rhogvx(ij,k,l) * rrhog(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                vy(ij,k,l) = rhogvy(ij,k,l) * rrhog(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                vz(ij,k,l) = rhogvz(ij,k,l) * rrhog(ij,k,l)
             enddo
          enddo
!OCL SERIAL
          do nq = 1, TRC_VMAX
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
                do ij = 1, ADM_gall
                   q(ij,k,l,nq) = rhogq(ij,k,l,nq) * rrhog(ij,k,l)
                enddo
             enddo
          enddo
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
             do ij = 1, ADM_gall
                qd(ij,k,l) = 1.D0
             enddo
          enddo
!OCL SERIAL
          do nq = NQW_STR, NQW_END
!OCL PARALLEL
             do k  = 1, ADM_kall
             do ij = 1, ADM_gall
                qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
             enddo
             enddo
          enddo
!OCL PARALLEL
          do k  = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
             do ij = 1, ADM_gall
                cv(ij,k,l) = qd(ij,k,l) * CNST_CV
             enddo
          enddo
!OCL SERIAL
          do nq = NQW_STR, NQW_END
!OCL PARALLEL
             do k  = 1, ADM_kall
             do ij = 1, ADM_gall
                cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
             enddo
             enddo
          enddo
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                tem(ij,k,l) = ein(ij,k,l) / cv(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                pre(ij,k,l) = rho(ij,k,l) * tem(ij,k,l) * ( qd(ij,k,l)*CNST_RAIR + q(ij,k,l,I_QV)*CNST_RVAP )
             enddo
          enddo
!OCL PARALLEL
          do k = ADM_kmin+1, ADM_kmax
             do ij = 1, ADM_gall
                w(ij,k,l) = rhogw(ij,k,l) &
                          / ( 0.5D0 * ( GRD_afac(k) * rho(ij,k  ,l) &
                                      + GRD_bfac(k) * rho(ij,k-1,l) ) * VMTR_GSGAM2H(ij,k,l) )
             enddo
          enddo

          !--- boundary conditions
          call bndcnd_all( ADM_gall,            & !--- [IN]
                           vx    (:,:,l),       & !--- [INOUT]
                           vy    (:,:,l),       & !--- [INOUT]
                           vz    (:,:,l),       & !--- [INOUT]
                           w     (:,:,l),       & !--- [INOUT]
                           tem   (:,:,l),       & !--- [INOUT]
                           rho   (:,:,l),       & !--- [INOUT]
                           pre   (:,:,l),       & !--- [INOUT]
                           ein   (:,:,l),       & !--- [INOUT]
                           rhog  (:,:,l),       & !--- [INOUT]
                           rhogvx(:,:,l),       & !--- [INOUT]
                           rhogvy(:,:,l),       & !--- [INOUT]
                           rhogvz(:,:,l),       & !--- [INOUT]
                           rhogw (:,:,l),       & !--- [INOUT]
                           rhoge (:,:,l),       & !--- [INOUT]
                           phi         (:,:,l), & !--- [IN]
                           VMTR_GSGAM2 (:,:,l), & !--- [IN]
                           VMTR_GSGAM2H(:,:,l), & !--- [IN]
                           VMTR_GZXH   (:,:,l), & !--- [IN]
                           VMTR_GZYH   (:,:,l), & !--- [IN]
                           VMTR_GZZH   (:,:,l)  ) !--- [IN]

             !--- perturbations ( pred, rhod, temd )
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                pregd(ij,k,l) = ( pre(ij,k,l) - pre_bs(ij,k,l) ) * VMTR_GSGAM2(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rhogd(ij,k,l) = ( rho(ij,k,l) - rho_bs(ij,k,l) ) * VMTR_GSGAM2(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                temd(ij,k,l) = tem(ij,k,l) - tem_bs(ij,k,l)
             enddo
          enddo

          call thrmdyn_th( ADM_gall, th(:,:,l), tem(:,:,l), pre(:,:,l) )

          call thrmdyn_eth( ADM_gall, eth(:,:,l), ein(:,:,l), pre(:,:,l), rho(:,:,l) )       

       enddo ! region LOOP

!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
             do ij = 1, ADM_gall_pl
                rho_pl(ij,k,l) = 0.D0
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then

          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   rrhog_pl(ij,k,l) = 1.D0 / rhog_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                  rho_pl(ij,k,l)  = rhog_pl(ij,k,l) / VMTR_GSGAM2_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   ein_pl(ij,k,l) = rhoge_pl(ij,k,l)  * rrhog_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   vx_pl(ij,k,l)  = rhogvx_pl(ij,k,l) * rrhog_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   vy_pl(ij,k,l)  = rhogvy_pl(ij,k,l) * rrhog_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   vz_pl(ij,k,l)  = rhogvz_pl(ij,k,l) * rrhog_pl(ij,k,l)
                enddo
             enddo
!OCL SERIAL
             do nq = 1, TRC_VMAX
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      q_pl(ij,k,l,nq) = rhogq_pl(ij,k,l,nq) * rrhog_pl(ij,k,l)
                   enddo
                enddo
             enddo
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   qd_pl(ij,k,l) = 1.D0
                enddo
             enddo
!OCL SERIAL
             do nq = NQW_STR, NQW_END
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      qd_pl(ij,k,l) = qd_pl(ij,k,l) - q_pl(ij,k,l,nq)
                   enddo
                enddo
             enddo
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   cv_pl(ij,k,l) = qd_pl(ij,k,l) * CNST_CV
                enddo
             enddo
!OCL SERIAL
             do nq = NQW_STR, NQW_END
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      cv_pl(ij,k,l) = cv_pl(ij,k,l) + q_pl(ij,k,l,nq) * CVW(nq)
                   enddo
                enddo
             enddo
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   tem_pl(ij,k,l) = ein_pl(ij,k,l) / cv_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   pre_pl(ij,k,l) = rho_pl(ij,k,l) * tem_pl(ij,k,l) &
                        * ( qd_pl(ij,k,l)*CNST_RAIR+q_pl(ij,k,l,I_QV)*CNST_RVAP )
                enddo
             enddo
!OCL PARALLEL
             do k = ADM_kmin+1, ADM_kmax
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   w_pl(ij,k,l) = rhogw_pl(ij,k,l) &
                                / ( 0.5D0 * ( GRD_afac(k) * rho_pl(ij,k  ,l) &
                                            + GRD_bfac(k) * rho_pl(ij,k-1,l) ) * VMTR_GSGAM2H_pl(ij,k,l) )
                enddo
             enddo

             call bndcnd_all( ADM_gall_pl,            & !--- [IN]
                              vx_pl    (:,:,l),       & !--- [INOUT]
                              vy_pl    (:,:,l),       & !--- [INOUT]
                              vz_pl    (:,:,l),       & !--- [INOUT]
                              w_pl     (:,:,l),       & !--- [INOUT]
                              tem_pl   (:,:,l),       & !--- [INOUT]
                              rho_pl   (:,:,l),       & !--- [INOUT]
                              pre_pl   (:,:,l),       & !--- [INOUT]
                              ein_pl   (:,:,l),       & !--- [INOUT]
                              rhog_pl  (:,:,l),       & !--- [INOUT]
                              rhogvx_pl(:,:,l),       & !--- [INOUT]
                              rhogvy_pl(:,:,l),       & !--- [INOUT]
                              rhogvz_pl(:,:,l),       & !--- [INOUT]
                              rhogw_pl (:,:,l),       & !--- [INOUT]
                              rhoge_pl (:,:,l),       & !--- [INOUT]
                              phi_pl   (:,:,l),       & !--- [IN]
                              VMTR_GSGAM2_pl (:,:,l), & !--- [IN]
                              VMTR_GSGAM2H_pl(:,:,l), & !--- [IN]
                              VMTR_GZXH_pl   (:,:,l), & !--- [IN]
                              VMTR_GZYH_pl   (:,:,l), & !--- [IN]
                              VMTR_GZZH_pl   (:,:,l)  ) !--- [IN]

!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   pregd_pl(ij,k,l) = ( pre_pl(ij,k,l) - pre_bs_pl(ij,k,l) ) * VMTR_GSGAM2_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   rhogd_pl(ij,k,l) = ( rho_pl(ij,k,l) - rho_bs_pl(ij,k,l) ) * VMTR_GSGAM2_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   temd_pl(ij,k,l) = tem_pl(ij,k,l) - tem_bs_pl(ij,k,l)
                enddo
             enddo

             call thrmdyn_th( ADM_gall_pl, th_pl(:,:,l), tem_pl(:,:,l), pre_pl(:,:,l) )

             call thrmdyn_eth( ADM_gall_pl, eth_pl(:,:,l), ein_pl(:,:,l), pre_pl(:,:,l), rho_pl(:,:,l) )       

          enddo
       endif

#ifdef _FJTIMER_
       call timer_end(2003)
       call timer_sta(2004)
#endif

       !---< numerical diffusion term

       if ( NDIFF_LOCATION == 'IN_LARGE_STEP' ) then

          if ( nl == 1 ) then ! only first step
!OCL SERIAL
             do l = 1, ADM_lall
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhog  (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogvx(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogvy(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogvz(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogw (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhoge (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogetot(ij,k,l) = 0.D0
                   enddo
                enddo
             enddo
!OCL SERIAL
             do nq = 1,TRC_VMAX
!OCL SERIAL
             do l = 1, ADM_lall
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
                   do ij = 1, ADM_gall
                      frhogq(ij,k,l,nq) = 0.D0
                   enddo
                enddo
             enddo
             enddo

!OCL SERIAL
             do l = 1, ADM_lall_pl
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhog_pl  (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogvx_pl(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogvy_pl(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogvz_pl(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogw_pl (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhoge_pl (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogetot_pl(ij,k,l) = 0.D0
                   enddo
                enddo
             enddo
!OCL SERIAL
             do nq = 1, TRC_VMAX
!OCL SERIAL
             do l  = 1, ADM_lall_pl
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogq_pl(ij,k,l,nq) = 0.D0
                   enddo
                enddo
             enddo
             enddo

             if ( NUMFILTER_DOrayleigh ) then ! [add] H.Yashiro 20120530

             !------ rayleigh damping
             call numfilter_rayleigh_damping( rho,     rho_pl,     & !--- [IN]
                                              vx,      vx_pl,      & !--- [IN]
                                              vy,      vy_pl,      & !--- [IN]
                                              vz,      vz_pl,      & !--- [IN]
                                              w,       w_pl,       & !--- [IN]
                                              frhogvx, frhogvx_pl, & !--- [INOUT]
                                              frhogvy, frhogvy_pl, & !--- [INOUT]
                                              frhogvz, frhogvz_pl, & !--- [INOUT]
                                              frhogw,  frhogw_pl   ) !--- [INOUT]

             endif

             !------ numerical diffusion
             call numfilter_numerical_hdiff( rho,       rho_pl,       & !--- [IN]
                                             vx,        vx_pl,        & !--- [IN]
                                             vy,        vy_pl,        & !--- [IN]
                                             vz,        vz_pl,        & !--- [IN]
                                             w,         w_pl,         & !--- [IN]
                                             temd,      temd_pl,      & !--- [IN]
                                             q,         q_pl,         & !--- [IN]
                                             frhog,     frhog_pl,     & !--- [INOUT]
                                             frhogvx,   frhogvx_pl,   & !--- [INOUT]
                                             frhogvy,   frhogvy_pl,   & !--- [INOUT]
                                             frhogvz,   frhogvz_pl,   & !--- [INOUT]
                                             frhogw,    frhogw_pl,    & !--- [INOUT]
                                             frhoge,    frhoge_pl,    & !--- [INOUT]
                                             frhogetot, frhogetot_pl, & !--- [INOUT]
                                             frhogq,    frhogq_pl     ) !--- [INOUT]

             if ( NUMFILTER_DOverticaldiff ) then ! [add] H.Yashiro 20120530

             call numfilter_numerical_vdiff( rho,       rho_pl,       & !--- [IN]
                                             vx,        vx_pl,        & !--- [IN]
                                             vy,        vy_pl,        & !--- [IN]
                                             vz,        vz_pl,        & !--- [IN]
                                             w,         w_pl,         & !--- [IN]
                                             temd,      temd_pl,      & !--- [IN]
                                             q,         q_pl,         & !--- [IN]
                                             frhog,     frhog_pl,     & !--- [INOUT]
                                             frhogvx,   frhogvx_pl,   & !--- [INOUT]
                                             frhogvy,   frhogvy_pl,   & !--- [INOUT]
                                             frhogvz,   frhogvz_pl,   & !--- [INOUT]
                                             frhogw,    frhogw_pl,    & !--- [INOUT]
                                             frhoge,    frhoge_pl,    & !--- [INOUT]
                                             frhogetot, frhogetot_pl, & !--- [INOUT]
                                             frhogq,    frhogq_pl     ) !--- [INOUT]

             endif

          endif

       elseif( NDIFF_LOCATION == 'IN_LARGE_STEP2' ) then

!OCL SERIAL
          do l = 1, ADM_lall
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhog  (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogvx(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogvy(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogvz(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogw (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhoge (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogetot(ij,k,l) = 0.D0
                enddo
             enddo
          enddo
!OCL SERIAL
          do nq = 1, TRC_VMAX
!OCL SERIAL
          do l  = 1, ADM_lall
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
                do ij = 1, ADM_gall
                   frhogq(ij,k,l,nq) = 0.D0
                enddo
             enddo
          enddo
          enddo

!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhog_pl  (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogvx_pl(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogvy_pl(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogvz_pl(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogw_pl (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhoge_pl (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogetot_pl(ij,k,l) = 0.D0
                enddo
             enddo
          enddo
!OCL SERIAL
          do nq = 1, TRC_VMAX
!OCL SERIAL
          do l  = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogq_pl(ij,k,l,nq) = 0.D0
                enddo
             enddo
          enddo
          enddo

          if ( NUMFILTER_DOrayleigh ) then ! [add] H.Yashiro 20120530

             !------ rayleigh damping
             call numfilter_rayleigh_damping( rho,     rho_pl,     & !--- [IN]
                                              vx,      vx_pl,      & !--- [IN]
                                              vy,      vy_pl,      & !--- [IN]
                                              vz,      vz_pl,      & !--- [IN]
                                              w,       w_pl,       & !--- [IN]
                                              frhogvx, frhogvx_pl, & !--- [INOUT]
                                              frhogvy, frhogvy_pl, & !--- [INOUT]
                                              frhogvz, frhogvz_pl, & !--- [INOUT]
                                              frhogw,  frhogw_pl   ) !--- [INOUT]

          endif

          !------ numerical diffusion
          call numfilter_numerical_hdiff( rho,       rho_pl,       & !--- [IN]
                                          vx,        vx_pl,        & !--- [IN]
                                          vy,        vy_pl,        & !--- [IN]
                                          vz,        vz_pl,        & !--- [IN]
                                          w,         w_pl,         & !--- [IN]
                                          temd,      temd_pl,      & !--- [IN]
                                          q,         q_pl,         & !--- [IN]
                                          frhog,     frhog_pl,     & !--- [INOUT]
                                          frhogvx,   frhogvx_pl,   & !--- [INOUT]
                                          frhogvy,   frhogvy_pl,   & !--- [INOUT]
                                          frhogvz,   frhogvz_pl,   & !--- [INOUT]
                                          frhogw,    frhogw_pl,    & !--- [INOUT]
                                          frhoge,    frhoge_pl,    & !--- [INOUT]
                                          frhogetot, frhogetot_pl, & !--- [INOUT]
                                          frhogq,    frhogq_pl     ) !--- [INOUT]

          if ( NUMFILTER_DOverticaldiff ) then ! [add] H.Yashiro 20120530

             call numfilter_numerical_vdiff( rho,       rho_pl,       & !--- [IN]
                                             vx,        vx_pl,        & !--- [IN]
                                             vy,        vy_pl,        & !--- [IN]
                                             vz,        vz_pl,        & !--- [IN]
                                             w,         w_pl,         & !--- [IN]
                                             temd,      temd_pl,      & !--- [IN]
                                             q,         q_pl,         & !--- [IN]
                                             frhog,     frhog_pl,     & !--- [INOUT]
                                             frhogvx,   frhogvx_pl,   & !--- [INOUT]
                                             frhogvy,   frhogvy_pl,   & !--- [INOUT]
                                             frhogvz,   frhogvz_pl,   & !--- [INOUT]
                                             frhogw,    frhogw_pl,    & !--- [INOUT]
                                             frhoge,    frhoge_pl,    & !--- [INOUT]
                                             frhogetot, frhogetot_pl, & !--- [INOUT]
                                             frhogq,    frhogq_pl     ) !--- [INOUT]

          endif

       else !--- OUT_LARGE

!OCL SERIAL
          do l = 1, ADM_lall
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhog  (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogvx(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogvy(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogvz(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogw (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhoge (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do ij = 1, ADM_gall
                   frhogetot(ij,k,l) = 0.D0
                enddo
             enddo
          enddo
!OCL SERIAL
          do nq = 1, TRC_VMAX
!OCL SERIAL
          do l  = 1, ADM_lall
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
                do ij = 1, ADM_gall
                   frhogq(ij,k,l,nq) = 0.D0
                enddo
             enddo
          enddo
          enddo

!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhog_pl  (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogvx_pl(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogvy_pl(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogvz_pl(ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogw_pl (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhoge_pl (ij,k,l) = 0.D0
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogetot_pl(ij,k,l) = 0.D0
                enddo
             enddo
          enddo
!OCL SERIAL
          do nq = 1, TRC_VMAX
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   frhogq_pl(ij,k,l,nq) = 0.D0
                enddo
             enddo
          enddo
          enddo

       endif

#ifdef _FJTIMER_
       call timer_end(2004)
       call timer_sta(2005)
#endif

       ! Smagorinksy-type SGS model [add] A.Noda 10/11/29
       if ( trim(TB_TYPE ) == 'SMG') then

         if ( ADM_prc_me /= ADM_prc_pl ) then ! for safety
            th_pl(:,:,:) = 0.D0
            vx_pl(:,:,:) = 0.D0
            vy_pl(:,:,:) = 0.D0
            vz_pl(:,:,:) = 0.D0
            w_pl (:,:,:) = 0.D0
         endif

         call tb_smg_driver( nl,              &
                             rho,   rho_pl,           & !--- [IN] : density
                             rhog,  rhog_pl,          & !--- [IN] : density
                             rhogq, rhogq_pl,         & !--- [IN] : tracers
                             vx,    vx_pl,            & !--- [IN] : Vx
                             vy,    vy_pl,            & !--- [IN] : Vy
                             vz,    vz_pl,            & !--- [IN] : Vz
                             w,     w_pl,             & !--- [IN] : w
                             !temd, temd_pl,          & !--- [IN] : temperature (deviation from tem_bs)
                             tem,   tem_pl,           & !--- [IN] : temperature
                             q,     q_pl,             & !--- [IN] : q
                             th,    th_pl,            & !--- [IN] : potential temperature 
                             frhog,     frhog_pl,     & !--- [INOUT] : tend. of rhog
                             frhogvx,   frhogvx_pl,   & !--- [INOUT] : tend. of rhogvx
                             frhogvy,   frhogvy_pl,   & !--- [INOUT] : tend. of rhogvy
                             frhogvz,   frhogvz_pl,   & !--- [INOUT] : tend. of rhogvz
                             frhogw,    frhogw_pl,    & !--- [INOUT] : tend. of rhogw
                             frhoge,    frhoge_pl,    & !--- [INOUT] : tend. of rhoge
                             frhogetot, frhogetot_pl, & !--- [INOUT] : tend. of rhoge
                             frhogq,    frhogq_pl     ) !--- [INOUT] : tend. of rhogq
       endif

       !--- Nudging routines [add] Y.Niwa 08/09/09
       if ( FLAG_NUDGING ) then

          if ( nl == 1 ) then
             call ndg_update_var
          endif

          if ( nl == num_of_iteration_lstep ) then
             call ndg_nudging_uvtp(  &
                  rho,       rho_pl,       & !--- [IN]
                  vx,        vx_pl,        & !--- [IN]
                  vy,        vy_pl,        & !--- [IN]
                  vz,        vz_pl,        & !--- [IN]
                  w,         w_pl,         & !--- [IN]
                  tem,       tem_pl,       & !--- [IN]
                  pre,       pre_pl,       & !--- [IN]
!                  frhog,     frhog_pl,     & !--- [INOUT] !S.Iga del 090908
                  frhogvx,   frhogvx_pl,   & !--- [INOUT]
                  frhogvy,   frhogvy_pl,   & !--- [INOUT]
                  frhogvz,   frhogvz_pl,   & !--- [INOUT]
                  frhogw,    frhogw_pl,    & !--- [INOUT]
                  frhoge,    frhoge_pl,    & !--- [INOUT]
                  frhogetot, frhogetot_pl, & !--- [INOUT]
                  .true.                   ) !--- [IN] ( tendency out )
          else
             call ndg_nudging_uvtp(  &
                  rho,       rho_pl,       & !--- [IN]
                  vx,        vx_pl,        & !--- [IN]
                  vy,        vy_pl,        & !--- [IN]
                  vz,        vz_pl,        & !--- [IN]
                  w,         w_pl,         & !--- [IN]
                  tem,       tem_pl,       & !--- [IN]
                  pre,       pre_pl,       & !--- [IN]
!                  frhog,     frhog_pl,     & !--- [INOUT] !S.Iga del 090908
                  frhogvx,   frhogvx_pl,   & !--- [INOUT]
                  frhogvy,   frhogvy_pl,   & !--- [INOUT]
                  frhogvz,   frhogvz_pl,   & !--- [INOUT]
                  frhogw,    frhogw_pl,    & !--- [INOUT]
                  frhoge,    frhoge_pl,    & !--- [INOUT]
                  frhogetot, frhogetot_pl, & !--- [INOUT]
                  .false.                  ) !--- [IN] ( tendency out )
          endif
       endif

       !--- calculation of advection tendency including Coriolis force
       call src_advection_convergence_v(  &  
            vx, vx_pl,                    & !--- [IN]
            vy, vy_pl,                    & !--- [IN]
            vz, vz_pl,                    & !--- [IN]
            w,  w_pl,                     & !--- [IN]
            rhog, rhog_pl,                & !--- [IN]
            rhogvx, rhogvx_pl,            & !--- [IN]
            rhogvy, rhogvy_pl,            & !--- [IN]
            rhogvz, rhogvz_pl,            & !--- [IN]
            rhogw,  rhogw_pl,             & !--- [IN]
            grhogvx, grhogvx_pl,          & !--- [OUT]
            grhogvy, grhogvy_pl,          & !--- [OUT]
            grhogvz, grhogvz_pl,          & !--- [OUT]
            grhogw,  grhogw_pl            ) !--- [OUT]

       !--- sum the large step tendency ( advection + coriolis + num.diff.,SGS,nudge )
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k  = 1, ADM_kall
          do ij = 1, ADM_gall
             grhog  (ij,k,l)   =                   frhog  (ij,k,l)
             grhogvx(ij,k,l)   = grhogvx(ij,k,l) + frhogvx(ij,k,l)
             grhogvy(ij,k,l)   = grhogvy(ij,k,l) + frhogvy(ij,k,l)
             grhogvz(ij,k,l)   = grhogvz(ij,k,l) + frhogvz(ij,k,l)
             grhogw (ij,k,l)   = grhogw (ij,k,l) + frhogw (ij,k,l)
             grhoge (ij,k,l)   =                   frhoge (ij,k,l)
             grhogetot(ij,k,l) =                   frhogetot(ij,k,l)
          enddo
          enddo
       enddo

       if (ADM_prc_pl==ADM_prc_me) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k  = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                grhog_pl  (ij,k,l)   =                      frhog_pl  (ij,k,l)
                grhogvx_pl(ij,k,l)   = grhogvx_pl(ij,k,l) + frhogvx_pl(ij,k,l)
                grhogvy_pl(ij,k,l)   = grhogvy_pl(ij,k,l) + frhogvy_pl(ij,k,l)
                grhogvz_pl(ij,k,l)   = grhogvz_pl(ij,k,l) + frhogvz_pl(ij,k,l)
                grhogw_pl (ij,k,l)   = grhogw_pl (ij,k,l) + frhogw_pl (ij,k,l)
                grhoge_pl (ij,k,l)   =                      frhoge_pl (ij,k,l)
                grhogetot_pl(ij,k,l) =                      frhogetot_pl(ij,k,l)
             enddo
             enddo
          enddo
       endif

#ifdef _FJTIMER_
    call timer_end(2005)
    call timer_sta(2006)
#endif
       !------------------------------------------------------------------------
       !> SMALL step
       !------------------------------------------------------------------------
!--->[start] for less diff
if ( nl == 1 ) then
    !--- get from prg0
    call prgvar_get( rhog,   rhog_pl,   & !--- [OUT]
                     rhogvx, rhogvx_pl, & !--- [OUT]
                     rhogvy, rhogvy_pl, & !--- [OUT]
                     rhogvz, rhogvz_pl, & !--- [OUT]
                     rhogw,  rhogw_pl,  & !--- [OUT]
                     rhoge,  rhoge_pl,  & !--- [OUT]
                     rhogq,  rhogq_pl,  & !--- [OUT]
                     0                  ) !--- [IN]
else
    !--- get from prg0
    call prgvar_get( rhog,   rhog_pl,   & !--- [OUT]
                     rhogvx, rhogvx_pl, & !--- [OUT]
                     rhogvy, rhogvy_pl, & !--- [OUT]
                     rhogvz, rhogvz_pl, & !--- [OUT]
                     rhogw,  rhogw_pl,  & !--- [OUT]
                     rhoge,  rhoge_pl,  & !--- [OUT]
                     rhogq,  rhogq_pl,  & !--- [OUT]
                     1                  ) !--- [IN]
endif
!<---[end] for less diff

       if ( nl /= 1 ) then ! update split values

          call prgvar_get_noq( rhog_split,   rhog_split_pl,   & !--- [OUT]
                               rhogvx_split, rhogvx_split_pl, & !--- [OUT]
                               rhogvy_split, rhogvy_split_pl, & !--- [OUT]
                               rhogvz_split, rhogvz_split_pl, & !--- [OUT]
                               rhogw_split,  rhogw_split_pl,  & !--- [OUT]
                               rhoge_split,  rhoge_split_pl   ) !--- [IN]

!OCL SERIAL
          do l = 1, ADM_lall
!OCL PARALLEL
             do k  = 1, ADM_kall
             do ij = 1, ADM_gall
                rhog_split  (ij,k,l) = rhog_split  (ij,k,l) - rhog  (ij,k,l)
                rhogvx_split(ij,k,l) = rhogvx_split(ij,k,l) - rhogvx(ij,k,l)
                rhogvy_split(ij,k,l) = rhogvy_split(ij,k,l) - rhogvy(ij,k,l)
                rhogvz_split(ij,k,l) = rhogvz_split(ij,k,l) - rhogvz(ij,k,l)
                rhogw_split (ij,k,l) = rhogw_split (ij,k,l) - rhogw (ij,k,l)
                rhoge_split (ij,k,l) = rhoge_split (ij,k,l) - rhoge (ij,k,l)
             enddo
             enddo
          enddo
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                rhog_split_pl  (ij,k,l) = rhog_split_pl  (ij,k,l) - rhog_pl  (ij,k,l)
                rhogvx_split_pl(ij,k,l) = rhogvx_split_pl(ij,k,l) - rhogvx_pl(ij,k,l)
                rhogvy_split_pl(ij,k,l) = rhogvy_split_pl(ij,k,l) - rhogvy_pl(ij,k,l)
                rhogvz_split_pl(ij,k,l) = rhogvz_split_pl(ij,k,l) - rhogvz_pl(ij,k,l)
                rhogw_split_pl (ij,k,l) = rhogw_split_pl (ij,k,l) - rhogw_pl (ij,k,l)
                rhoge_split_pl (ij,k,l) = rhoge_split_pl (ij,k,l) - rhoge_pl (ij,k,l)
             enddo
             enddo
          enddo

       endif

!--->[start] for less diff
do l = 1, ADM_lall
             call bndcnd_all(          &
                  ADM_gall,            &  !--- in
                  vx(:,:,l),           &  !--- inout
                  vy(:,:,l),           &  !--- inout
                  vz(:,:,l),           &  !--- inout
                  w(:,:,l),            &  !--- inout
                  tem(:,:,l),          &  !--- inout
                  rho(:,:,l),          &  !--- inout
                  pre(:,:,l),          &  !--- inout
                  ein(:,:,l),          &  !--- inout
                  rhog(:,:,l),         &  !--- inout
                  rhogvx(:,:,l),       &  !--- inout
                  rhogvy(:,:,l),       &  !--- inout
                  rhogvz(:,:,l),       &  !--- inout
                  rhogw(:,:,l),        &  !--- inout
                  rhoge(:,:,l),        &  !--- inout
                  phi(:,:,l),          &  !--- in
                  VMTR_GSGAM2(:,:,l),  &  !--- in
                  VMTR_GSGAM2H(:,:,l), &  !--- in
                  VMTR_GZXH(:,:,l),    &  !--- in
                  VMTR_GZYH(:,:,l),    &  !--- in
                  VMTR_GZZH(:,:,l)     )  !--- in
enddo

if ( ADM_prc_pl == ADM_prc_me ) then
do l = 1, ADM_lall_pl
                call bndcnd_all(          &
                     ADM_gall_pl,            &  !--- in
                     vx_pl(:,:,l),           &  !--- inout
                     vy_pl(:,:,l),           &  !--- inout
                     vz_pl(:,:,l),           &  !--- inout
                     w_pl(:,:,l),            &  !--- inout
                     tem_pl(:,:,l),          &  !--- inout
                     rho_pl(:,:,l),          &  !--- inout
                     pre_pl(:,:,l),          &  !--- inout
                     ein_pl(:,:,l),          &  !--- inout
                     rhog_pl(:,:,l),         &  !--- inout
                     rhogvx_pl(:,:,l),       &  !--- inout
                     rhogvy_pl(:,:,l),       &  !--- inout
                     rhogvz_pl(:,:,l),       &  !--- inout
                     rhogw_pl(:,:,l),        &  !--- inout
                     rhoge_pl(:,:,l),        &  !--- inout
                     phi_pl(:,:,l),          &  !--- in
                     VMTR_GSGAM2_pl(:,:,l),  &  !--- in
                     VMTR_GSGAM2H_pl(:,:,l), &  !--- in
                     VMTR_GZXH_pl(:,:,l),    &  !--- in
                     VMTR_GZYH_pl(:,:,l),    &  !--- in
                     VMTR_GZZH_pl(:,:,l)     )  !--- in
enddo
endif
!--->[end] for less diff

#ifdef _FJTIMER_
       call timer_end(2006)
       call timer_sta(2007)
#endif

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

       call vi_small_step( rhog,         rhog_pl,         & !--- [INOUT] prognostic value
                           rhogvx,       rhogvx_pl,       & !--- [INOUT]
                           rhogvy,       rhogvy_pl,       & !--- [INOUT]
                           rhogvz,       rhogvz_pl,       & !--- [INOUT]
                           rhogw,        rhogw_pl,        & !--- [INOUT]
                           rhoge,        rhoge_pl,        & !--- [INOUT]
                           vx,           vx_pl,           & !--- [IN] diagnostic value
                           vy,           vy_pl,           & !--- [IN]
                           vz,           vz_pl,           & !--- [IN]
                           eth,          eth_pl,          & !--- [IN]
                           rhogd,        rhogd_pl,        & !--- [IN]
                           pregd,        pregd_pl,        & !--- [IN]
                           grhog,        grhog_pl,        & !--- [IN] large step tendency
                           grhogvx,      grhogvx_pl,      & !--- [IN]
                           grhogvy,      grhogvy_pl,      & !--- [IN]
                           grhogvz,      grhogvz_pl,      & !--- [IN]
                           grhogw,       grhogw_pl,       & !--- [IN]
                           grhoge,       grhoge_pl,       & !--- [IN]
                           grhogetot,    grhogetot_pl,    & !--- [IN]
                           rhog_split,   rhog_split_pl,   & !--- [INOUT] split value
                           rhogvx_split, rhogvx_split_pl, & !--- [INOUT]
                           rhogvy_split, rhogvy_split_pl, & !--- [INOUT]
                           rhogvz_split, rhogvz_split_pl, & !--- [INOUT]
                           rhogw_split,  rhogw_split_pl,  & !--- [INOUT]
                           rhoge_split,  rhoge_split_pl,  & !--- [INOUT]
                           v_mean_c,     v_mean_c_pl,     & !--- [OUT] mean value
                           small_step_ite,                & !--- [IN]
                           small_step_dt                  ) !--- [IN]

#ifdef _FJTIMER_
       call timer_end(2007)
       call timer_sta(2008)
#endif
       !------------------------------------------------------------------------
       !>  Tracer advection
       !------------------------------------------------------------------------
       if ( TRC_ADV_TYPE == 'MIURA2004' ) then

          if ( nl == num_of_iteration_lstep ) then

             call src_update_tracer_rev( TRC_VMAX,                                              & !--- [IN]
                                         rhogq,                    rhogq_pl,                    & !--- [INOUT]
                                         rhog_old,                 rhog_old_pl,                 & !--- [IN]
                                         v_mean_c(:,:,:,I_rhog),   v_mean_c_pl(:,:,:,I_rhog),   & !--- [IN]
                                         v_mean_c(:,:,:,I_rhogvx), v_mean_c_pl(:,:,:,I_rhogvx), & !--- [IN]
                                         v_mean_c(:,:,:,I_rhogvy), v_mean_c_pl(:,:,:,I_rhogvy), & !--- [IN]
                                         v_mean_c(:,:,:,I_rhogvz), v_mean_c_pl(:,:,:,I_rhogvz), & !--- [IN]
                                         v_mean_c(:,:,:,I_rhogw),  v_mean_c_pl(:,:,:,I_rhogw),  & !--- [IN]
                                         frhog,                    frhog_pl,                    & !--- [IN]
                                         TIME_DTL                                               ) !--- [IN]

!OCL SERIAL
             do nq = 1,TRC_VMAX
!OCL SERIAL
             do l  = 1, ADM_lall
!OCL PARALLEL
                do k  = 1, ADM_kall
                do ij = 1, ADM_gall
                   rhogq(ij,k,l,nq) = rhogq(ij,k,l,nq) + TIME_DTL * frhogq(ij,k,l,nq) ! update rhogq by viscosity
                enddo
                enddo
!OCL PARALLEL
                do ij = 1, ADM_gall
                   rhogq(ij,ADM_kmin-1,l,nq) = 0.D0
                   rhogq(ij,ADM_kmax+1,l,nq) = 0.D0
                enddo
             enddo
             enddo

             if ( ADM_prc_pl == ADM_prc_me ) then
!OCL SERIAL
                do nq = 1, TRC_VMAX
!OCL SERIAL
                do l  = 1, ADM_lall_pl
!OCL PARALLEL
                   do k  = 1, ADM_kall
                   do ij = 1, ADM_gall_pl
                      rhogq_pl(ij,k,l,nq) = rhogq_pl(ij,k,l,nq) + TIME_DTL * frhogq_pl(ij,k,l,nq)
                   enddo
                   enddo
!OCL PARALLEL
                   do ij = 1, ADM_gall_pl
                      rhogq_pl(ij,ADM_kmin-1,l,nq) = 0.D0
                      rhogq_pl(ij,ADM_kmax+1,l,nq) = 0.D0
                   enddo
                enddo
                enddo
             endif

             !<---- I don't recommend adding the hyperviscosity term
             !<---- because of numerical instability in this case. ( H.Tomita )

          endif ! Last large step only

       elseif( TRC_ADV_TYPE == 'MIURA2004OLD' ) then

          if ( nl == num_of_iteration_lstep ) then

             call src_update_tracer( rhogq,                    rhogq_pl,                    & !--- [INOUT]
                                     v_mean_c(:,:,:,I_rhog),   v_mean_c_pl(:,:,:,I_rhog),   & !--- [IN]
                                     v_mean_c(:,:,:,I_rhogvx), v_mean_c_pl(:,:,:,I_rhogvx), & !--- [IN]
                                     v_mean_c(:,:,:,I_rhogvy), v_mean_c_pl(:,:,:,I_rhogvy), & !--- [IN]
                                     v_mean_c(:,:,:,I_rhogvz), v_mean_c_pl(:,:,:,I_rhogvz), & !--- [IN]
                                     v_mean_c(:,:,:,I_rhogw),  v_mean_c_pl(:,:,:,I_rhogw),  & !--- [IN]
                                     TIME_DTL                                               ) !--- [IN]

!OCL SERIAL
             do nq = 1,TRC_VMAX
!OCL SERIAL
             do l  = 1, ADM_lall
!OCL PARALLEL
                do k  = 1, ADM_kall
                do ij = 1, ADM_gall
                   rhogq(ij,k,l,nq) = rhogq(ij,k,l,nq) + TIME_DTL * frhogq(ij,k,l,nq) ! update rhogq by viscosity
                enddo
                enddo
!OCL PARALLEL
                do ij = 1, ADM_gall
                   rhogq(ij,ADM_kmin-1,l,nq) = 0.D0
                   rhogq(ij,ADM_kmax+1,l,nq) = 0.D0
                enddo
             enddo
             enddo

             if ( ADM_prc_pl == ADM_prc_me ) then
!OCL SERIAL
                do nq = 1, TRC_VMAX
!OCL SERIAL
                do l  = 1, ADM_lall_pl
!OCL PARALLEL
                   do k  = 1, ADM_kall
                   do ij = 1, ADM_gall_pl
                      rhogq_pl(ij,k,l,nq) = rhogq_pl(ij,k,l,nq) + TIME_DTL * frhogq_pl(ij,k,l,nq)
                   enddo
                   enddo
!OCL PARALLEL
                   do ij = 1, ADM_gall_pl
                      rhogq_pl(ij,ADM_kmin-1,l,nq) = 0.D0
                      rhogq_pl(ij,ADM_kmax+1,l,nq) = 0.D0
                   enddo
                enddo
                enddo
             endif

             !<---- I don't recommend adding the hyperviscosity term
             !<---- because of numerical instability in this case. ( H.Tomita )

          endif ! Last large step only

       elseif( TRC_ADV_TYPE == 'DEFAULT' ) then

!OCL SERIAL
          do nq = 1, TRC_VMAX

             call src_advection_convergence( v_mean_c(:,:,:,I_rhogvx), v_mean_c_pl(:,:,:,I_rhogvx), & !--- [IN]
                                             v_mean_c(:,:,:,I_rhogvy), v_mean_c_pl(:,:,:,I_rhogvy), & !--- [IN]
                                             v_mean_c(:,:,:,I_rhogvz), v_mean_c_pl(:,:,:,I_rhogvz), & !--- [IN]
                                             v_mean_c(:,:,:,I_rhogw),  v_mean_c_pl(:,:,:,I_rhogw),  & !--- [IN]
                                             q(:,:,:,nq),      q_pl(:,:,:,nq),                      & !--- [IN]
                                             grhogq(:,:,:,nq), grhogq_pl(:,:,:,nq),                 & !--- [OUT]
                                             I_SRC_default                                          ) !--- [IN]  [mod] H.Yashiro 20120530

!OCL SERIAL
             do l  = 1, ADM_lall
!OCL PARALLEL
                do k  = 1, ADM_kall
                do ij = 1, ADM_gall
                   rhogq(ij,k,l,nq) = rhogq_old(ij,k,l,nq)                      &
                                    + ( num_of_iteration_sstep(nl) * TIME_DTS ) &
                                    * ( grhogq(ij,k,l,nq) + frhogq(ij,k,l,nq) )
                enddo
                enddo
!OCL PARALLEL
                do ij = 1, ADM_gall
                   rhogq(ij,ADM_kmin-1,l,nq) = 0.D0
                   rhogq(ij,ADM_kmax+1,l,nq) = 0.D0
                enddo
             enddo

             if ( ADM_prc_pl == ADM_prc_me ) then
!OCL SERIAL
                do l  = 1, ADM_lall_pl
!OCL PARALLEL
                   do k = 1, ADM_kall
                   do ij = 1, ADM_gall_pl
                      rhogq_pl(ij,k,l,nq) = rhogq_old_pl(ij,k,l,nq)                         &
                                          + ( num_of_iteration_sstep(nl) * TIME_DTS )       &
                                          * ( grhogq_pl(ij,k,l,nq) + frhogq_pl(ij,k,l,nq) )
                   enddo
                   enddo
!OCL PARALLEL
                   do ij = 1, ADM_gall_pl
                      rhogq_pl(ij,ADM_kmin-1,l,nq) = 0.D0
                      rhogq_pl(ij,ADM_kmax+1,l,nq) = 0.D0
                   enddo
                enddo
             endif

          enddo ! tracer LOOP

       endif

       !--- TKE fixer ( TKE >= 0.D0 )
       ! 2011/08/16 M.Satoh [comment] need this fixer for every small time steps
       if ( I_TKE >= 0 ) then
          if ( TRC_ADV_TYPE == 'DEFAULT' .OR. nl == num_of_iteration_lstep ) then
!OCL SERIAL
             do l  = 1, ADM_lall
!OCL PARALLEL
                do k  = 1, ADM_kall
                do ij = 1, ADM_gall
                   TKEg_corr = TKE_MIN * VMTR_GSGAM2(ij,k,l) - rhogq(ij,k,l,I_TKE)

                   if ( TKEg_corr >= 0.D0 ) then
                      rhoge(ij,k,l)       = rhoge(ij,k,l)       - TKEg_corr
                      rhogq(ij,k,l,I_TKE) = rhogq(ij,k,l,I_TKE) + TKEg_corr
                   endif
                enddo
                enddo
             enddo

             if ( ADM_prc_pl == ADM_prc_me ) then
!OCL SERIAL
                do l  = 1, ADM_lall_pl
!OCL PARALLEL
                   do k  = 1, ADM_kall
                   do ij = 1, ADM_gall_pl
                      TKEg_corr = TKE_MIN * VMTR_GSGAM2_pl(ij,k,l) - rhogq_pl(ij,k,l,I_TKE)

                      if ( TKEg_corr >= 0.D0 ) then
                         rhoge_pl(ij,k,l)       = rhoge_pl(ij,k,l)       - TKEg_corr
                         rhogq_pl(ij,k,l,I_TKE) = rhogq_pl(ij,k,l,I_TKE) + TKEg_corr
                      endif
                   enddo
                   enddo
                enddo
             endif

          endif
       endif

#ifdef _FJTIMER_
       call timer_end(2008)
#endif
       !------ Update
       if ( nl /= num_of_iteration_lstep ) then
          !------ Temporaly process --> prg1
          call prgvar_set( rhog,   rhog_pl,   & !--- [IN]
                           rhogvx, rhogvx_pl, & !--- [IN]
                           rhogvy, rhogvy_pl, & !--- [IN]
                           rhogvz, rhogvz_pl, & !--- [IN]
                           rhogw,  rhogw_pl,  & !--- [IN]
                           rhoge,  rhoge_pl,  & !--- [IN]
                           rhogq,  rhogq_pl,  & !--- [IN]
                           1                  ) !--- [IN]
       endif

    enddo !--- large step 

#ifdef _FJTIMER_
    call timer_end(2002)
    call timer_sta(2009)
#endif

    ! 2011/08/16 M.Satoh move from above =>
    ! 2010/05/06 M.Satoh (comment on TDK)
    ! 2011/08/16 M.Satoh bug fix for TDK: conv => tendency
    ! 2011/08/16 M.Satoh
    !    use rhogq-rhogq_old to calculate qconv & qtendency
    ! Kuo param.    : total convergence of water vapor [integ over time]
    ! Tiedtke scheme: tendency of water vapor [per time]
    if ( CP_TYPE == 'KUO' ) then
       diagvar(:,:,:,I_RHOGQV_CONV) = rhogq(:,:,:,I_QV) - rhogq_old(:,:,:,I_QV)
       if ( ADM_prc_pl == ADM_prc_me ) then
          diagvar_pl(:,:,:,I_RHOGQV_CONV) = rhogq_pl(:,:,:,I_QV) - rhogq_old_pl(:,:,:,I_QV)
       endif
    endif

    if ( CP_TYPE == 'TDK' ) then
       ! 2011.08.16 M.Satoh,bug fix: conv => tendency
       ! qv_dyn_tend = v grad q
       !             = ( div(rho v q) - div(rho v)*q )/rho
       diagvar(:,:,:,I_QV_DYN_TEND) = ( ( rhogq(:,:,:,I_QV)-rhogq_old(:,:,:,I_QV) ) &
                                      - ( rhog(:,:,:)-rhog_old(:,:,:) )             &
                                        * rhogq(:,:,:,I_QV) / rhog(:,:,:)           & ! = q(:,:,:,I_QV)
                                      ) / TIME_DTL / rhog(:,:,:)
       if ( ADM_prc_pl == ADM_prc_me ) then
          diagvar_pl(:,:,:,I_QV_DYN_TEND) = ( ( rhogq_pl(:,:,:,I_QV)-rhogq_old_pl(:,:,:,I_QV) ) &
                                            - ( rhog_pl(:,:,:)-rhog_old_pl(:,:,:) )             &
                                              * rhogq_pl(:,:,:,I_QV) / rhog_pl(:,:,:)           & ! = q_pl(:,:,:,I_QV)
                                            ) / TIME_DTL / rhog_pl(:,:,:)
       endif
    endif
    ! <= 2011/08/16 M.Satoh move from above

    call prgvar_set( rhog,   rhog_pl,   & !--- [IN]
                     rhogvx, rhogvx_pl, & !--- [IN]
                     rhogvy, rhogvy_pl, & !--- [IN]
                     rhogvz, rhogvz_pl, & !--- [IN]
                     rhogw,  rhogw_pl,  & !--- [IN]
                     rhoge,  rhoge_pl,  & !--- [IN]
                     rhogq,  rhogq_pl,  & !--- [IN]
                     0                  )

#ifdef _FJTIMER_
    call timer_end(2009)
    call timer_end(2000)
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
       VMTR_GZZH_pl
    use mod_runconf, only: &
       TRC_VMAX,       &
       I_QV,           &
       NQW_STR,        &
       NQW_END,        &
       CVW,            &
       NDIFF_LOCATION
    use mod_bsstate, only: &
       pre_bs, pre_bs_pl, &
       tem_bs, tem_bs_pl, &
       rho_bs, rho_bs_pl, &
       phi, phi_pl
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
    real(8) :: frhog(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhogvx  ( G^{1/2} X gamma2 )
    real(8) :: frhogvx(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhogvy  ( G^{1/2} X gamma2 )
    real(8) :: frhogvy(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: frhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- forcing tendensy of rhogvz  ( G^{1/2} X gamma2 )
    real(8) :: frhogvz(ADM_gall,   ADM_kall,ADM_lall   )
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
    !--- density deviation from the base state ( G^{1/2} X gamma2 )
    real(8) :: rhogd(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    !--- pressure deviation from the base state ( G^{1/2} X gamma2 )
    real(8) :: pregd(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: pregd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
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

    if ( NDIFF_LOCATION == 'OUT_LARGE_STEP' ) then

       !--- get from prg0
       call prgvar_get( rhog,   rhog_pl,   & !--- [OUT]
                        rhogvx, rhogvx_pl, & !--- [OUT]
                        rhogvy, rhogvy_pl, & !--- [OUT]
                        rhogvz, rhogvz_pl, & !--- [OUT]
                        rhogw,  rhogw_pl,  & !--- [OUT]
                        rhoge,  rhoge_pl,  & !--- [OUT]
                        rhogq,  rhogq_pl,  & !--- [OUT]
                        0                  ) !--- [IN]

       !---< Generate diagnostic values and set the boudary conditions
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rrhog(ij,k,l) = 1.D0 / rhog(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
               rho(ij,k,l) = rhog(ij,k,l) / VMTR_GSGAM2(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                ein(ij,k,l) = rhoge(ij,k,l) * rrhog(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                vx(ij,k,l) = rhogvx(ij,k,l) * rrhog(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                vy(ij,k,l) = rhogvy(ij,k,l) * rrhog(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                vz(ij,k,l) = rhogvz(ij,k,l) * rrhog(ij,k,l)
             enddo
          enddo
!OCL SERIAL
          do nq = 1, TRC_VMAX
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
                do ij = 1, ADM_gall
                   q(ij,k,l,nq) = rhogq(ij,k,l,nq) * rrhog(ij,k,l)
                enddo
             enddo
          enddo
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
             do ij = 1, ADM_gall
                qd(ij,k,l) = 1.D0
             enddo
          enddo
!OCL SERIAL
          do nq = NQW_STR, NQW_END
!OCL PARALLEL
             do k  = 1, ADM_kall
             do ij = 1, ADM_gall
                qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
             enddo
             enddo
          enddo
!OCL PARALLEL
          do k  = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
             do ij = 1, ADM_gall
                cv(ij,k,l) = qd(ij,k,l) * CNST_CV
             enddo
          enddo
!OCL SERIAL
          do nq = NQW_STR, NQW_END
!OCL PARALLEL
             do k  = 1, ADM_kall
             do ij = 1, ADM_gall
                cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
             enddo
             enddo
          enddo
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                tem(ij,k,l) = ein(ij,k,l) / cv(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                pre(ij,k,l) = rho(ij,k,l) * tem(ij,k,l) * ( qd(ij,k,l)*CNST_RAIR + q(ij,k,l,I_QV)*CNST_RVAP )
             enddo
          enddo
!OCL PARALLEL
          do k = ADM_kmin+1, ADM_kmax
             do ij = 1, ADM_gall
                w(ij,k,l) = rhogw(ij,k,l) &
                          / ( 0.5D0 * ( GRD_afac(k) * rho(ij,k  ,l) &
                                      + GRD_bfac(k) * rho(ij,k-1,l) ) * VMTR_GSGAM2H(ij,k,l) )
             enddo
          enddo

          !--- boundary conditions
          call bndcnd_all( ADM_gall,            & !--- [IN]
                           vx    (:,:,l),       & !--- [INOUT]
                           vy    (:,:,l),       & !--- [INOUT]
                           vz    (:,:,l),       & !--- [INOUT]
                           w     (:,:,l),       & !--- [INOUT]
                           tem   (:,:,l),       & !--- [INOUT]
                           rho   (:,:,l),       & !--- [INOUT]
                           pre   (:,:,l),       & !--- [INOUT]
                           ein   (:,:,l),       & !--- [INOUT]
                           rhog  (:,:,l),       & !--- [INOUT]
                           rhogvx(:,:,l),       & !--- [INOUT]
                           rhogvy(:,:,l),       & !--- [INOUT]
                           rhogvz(:,:,l),       & !--- [INOUT]
                           rhogw (:,:,l),       & !--- [INOUT]
                           rhoge (:,:,l),       & !--- [INOUT]
                           phi         (:,:,l), & !--- [IN]
                           VMTR_GSGAM2 (:,:,l), & !--- [IN]
                           VMTR_GSGAM2H(:,:,l), & !--- [IN]
                           VMTR_GZXH   (:,:,l), & !--- [IN]
                           VMTR_GZYH   (:,:,l), & !--- [IN]
                           VMTR_GZZH   (:,:,l)  ) !--- [IN]

             !--- perturbations ( pred, rhod, temd )
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                temd(ij,k,l) = tem(ij,k,l) - tem_bs(ij,k,l)
             enddo
          enddo

       enddo ! region LOOP

       if ( ADM_prc_me == ADM_prc_pl ) then

          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   rrhog_pl(ij,k,l) = 1.D0 / rhog_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                  rho_pl(ij,k,l)  = rhog_pl(ij,k,l) / VMTR_GSGAM2_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   ein_pl(ij,k,l) = rhoge_pl(ij,k,l)  * rrhog_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   vx_pl(ij,k,l)  = rhogvx_pl(ij,k,l) * rrhog_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   vy_pl(ij,k,l)  = rhogvy_pl(ij,k,l) * rrhog_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   vz_pl(ij,k,l)  = rhogvz_pl(ij,k,l) * rrhog_pl(ij,k,l)
                enddo
             enddo
!OCL SERIAL
             do nq = 1, TRC_VMAX
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      q_pl(ij,k,l,nq) = rhogq_pl(ij,k,l,nq) * rrhog_pl(ij,k,l)
                   enddo
                enddo
             enddo
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   qd_pl(ij,k,l) = 1.D0
                enddo
             enddo
!OCL SERIAL
             do nq = NQW_STR, NQW_END
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      qd_pl(ij,k,l) = qd_pl(ij,k,l) - q_pl(ij,k,l,nq)
                   enddo
                enddo
             enddo
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   cv_pl(ij,k,l) = qd_pl(ij,k,l) * CNST_CV
                enddo
             enddo
!OCL SERIAL
             do nq = NQW_STR, NQW_END
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      cv_pl(ij,k,l) = cv_pl(ij,k,l) + q_pl(ij,k,l,nq) * CVW(nq)
                   enddo
                enddo
             enddo
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   tem_pl(ij,k,l) = ein_pl(ij,k,l) / cv_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   pre_pl(ij,k,l) = rho_pl(ij,k,l) * tem_pl(ij,k,l) &
                        * ( qd_pl(ij,k,l)*CNST_RAIR+q_pl(ij,k,l,I_QV)*CNST_RVAP )
                enddo
             enddo
!OCL PARALLEL
             do k = ADM_kmin+1, ADM_kmax
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   w_pl(ij,k,l) = rhogw_pl(ij,k,l) &
                                / ( 0.5D0 * ( GRD_afac(k) * rho_pl(ij,k  ,l) &
                                            + GRD_bfac(k) * rho_pl(ij,k-1,l) ) * VMTR_GSGAM2H_pl(ij,k,l) )
                enddo
             enddo

             call bndcnd_all( ADM_gall_pl,            & !--- [IN]
                              vx_pl    (:,:,l),       & !--- [INOUT]
                              vy_pl    (:,:,l),       & !--- [INOUT]
                              vz_pl    (:,:,l),       & !--- [INOUT]
                              w_pl     (:,:,l),       & !--- [INOUT]
                              tem_pl   (:,:,l),       & !--- [INOUT]
                              rho_pl   (:,:,l),       & !--- [INOUT]
                              pre_pl   (:,:,l),       & !--- [INOUT]
                              ein_pl   (:,:,l),       & !--- [INOUT]
                              rhog_pl  (:,:,l),       & !--- [INOUT]
                              rhogvx_pl(:,:,l),       & !--- [INOUT]
                              rhogvy_pl(:,:,l),       & !--- [INOUT]
                              rhogvz_pl(:,:,l),       & !--- [INOUT]
                              rhogw_pl (:,:,l),       & !--- [INOUT]
                              rhoge_pl (:,:,l),       & !--- [INOUT]
                              phi_pl   (:,:,l),       & !--- [IN]
                              VMTR_GSGAM2_pl (:,:,l), & !--- [IN]
                              VMTR_GSGAM2H_pl(:,:,l), & !--- [IN]
                              VMTR_GZXH_pl   (:,:,l), & !--- [IN]
                              VMTR_GZYH_pl   (:,:,l), & !--- [IN]
                              VMTR_GZZH_pl   (:,:,l)  ) !--- [IN]

!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   temd_pl(ij,k,l) = tem_pl(ij,k,l) - tem_bs_pl(ij,k,l)
                enddo
             enddo

          enddo
       endif

!OCL SERIAL
             do l = 1, ADM_lall
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhog  (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogvx(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogvy(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogvz(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogw (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhoge (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                   do ij = 1, ADM_gall
                      frhogetot(ij,k,l) = 0.D0
                   enddo
                enddo
             enddo
!OCL SERIAL
             do nq = 1,TRC_VMAX
             do l = 1, ADM_lall
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
                   do ij = 1, ADM_gall
                      frhogq(ij,k,l,nq) = 0.D0
                   enddo
                enddo
             enddo
             enddo

!OCL SERIAL
             do l = 1, ADM_lall_pl
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhog_pl  (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogvx_pl(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogvy_pl(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogvz_pl(ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogw_pl (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhoge_pl (ij,k,l) = 0.D0
                   enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogetot_pl(ij,k,l) = 0.D0
                   enddo
                enddo
             enddo
!OCL SERIAL
             do nq = 1, TRC_VMAX
             do l  = 1, ADM_lall_pl
!OCL PARALLEL
                do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                   do ij = 1, ADM_gall_pl
                      frhogq_pl(ij,k,l,nq) = 0.D0
                   enddo
                enddo
             enddo
             enddo

             if ( NUMFILTER_DOrayleigh ) then ! [add] H.Yashiro 20120530

             !------ rayleigh damping
             call numfilter_rayleigh_damping( &
                  rho, rho_pl,                & !--- [IN]
                  vx, vx_pl,                  & !--- [IN]
                  vy, vy_pl,                  & !--- [IN]
                  vz, vz_pl,                  & !--- [IN]
                  w, w_pl,                    & !--- [IN]
                  frhogvx, frhogvx_pl,        & !--- [INOUT]
                  frhogvy, frhogvy_pl,        & !--- [INOUT]
                  frhogvz, frhogvz_pl,        & !--- [INOUT]
                  frhogw, frhogw_pl           ) !--- [INOUT]

             endif

             !------ numerical diffusion
             call numfilter_numerical_hdiff(  &
                  rho, rho_pl,                & !--- [IN]
                  vx, vx_pl,                  & !--- [IN]
                  vy, vy_pl,                  & !--- [IN]
                  vz, vz_pl,                  & !--- [IN]
                  w , w_pl,                   & !--- [IN]
                  temd, temd_pl,              & !--- [IN]
                  q , q_pl,                   & !--- [IN]
                  frhog, frhog_pl,            & !--- [INOUT]
                  frhogvx, frhogvx_pl,        & !--- [INOUT]
                  frhogvy, frhogvy_pl,        & !--- [INOUT]
                  frhogvz, frhogvz_pl,        & !--- [INOUT]
                  frhogw , frhogw_pl,         & !--- [INOUT]
                  frhoge , frhoge_pl,         & !--- [INOUT]
                  frhogetot , frhogetot_pl,         & !--- [INOUT]
                  frhogq , frhogq_pl )          !--- [INOUT]

             if ( NUMFILTER_DOverticaldiff ) then ! [add] H.Yashiro 20120530

             call numfilter_numerical_vdiff(  &
                  rho,   rho_pl,              &  !--- [IN]
                  vx, vx_pl,                  &  !--- [IN]
                  vy, vy_pl,                  &  !--- [IN]
                  vz, vz_pl,                  &  !--- [IN]
                  w , w_pl,                   &  !--- [IN]
                  temd,temd_pl,               &  !--- [IN]
                  q,q_pl,                     &  !--- [IN]
                  frhog, frhog_pl,            &  !--- [INOUT]
                  frhogvx, frhogvx_pl,        &  !--- [INOUT]
                  frhogvy, frhogvy_pl,        &  !--- [INOUT]
                  frhogvz, frhogvz_pl,        &  !--- [INOUT]
                  frhogw , frhogw_pl,         &  !--- [INOUT]
                  frhoge , frhoge_pl,         &  !--- [INOUT]
                  frhogetot , frhogetot_pl,         & !--- [INOUT]
                  frhogq , frhogq_pl          )  !--- [INOUT]

             endif

       call OPRT_horizontalize_vec( frhogvx, frhogvx_pl, & !--- [INOUT]
                                    frhogvy, frhogvy_pl, & !--- [INOUT]
                                    frhogvz, frhogvz_pl  ) !--- [INOUT]

       rhog   = rhog   + TIME_DTL * frhog
       rhogvx = rhogvx + TIME_DTL * frhogvx
       rhogvy = rhogvy + TIME_DTL * frhogvy
       rhogvz = rhogvz + TIME_DTL * frhogvz
       rhogw  = rhogw  + TIME_DTL * frhogw
       rhoge  = rhoge  + TIME_DTL * frhoge
       rhogq  = rhogq  + TIME_DTL * frhogq

       if ( ADM_prc_me == ADM_prc_pl ) then
          rhog_pl   = rhog_pl   + TIME_DTL * frhog_pl
          rhogvx_pl = rhogvx_pl + TIME_DTL * frhogvx_pl
          rhogvy_pl = rhogvy_pl + TIME_DTL * frhogvy_pl
          rhogvz_pl = rhogvz_pl + TIME_DTL * frhogvz_pl
          rhogw_pl  = rhogw_pl  + TIME_DTL * frhogw_pl
          rhoge_pl  = rhoge_pl  + TIME_DTL * frhoge_pl
          rhogq_pl  = rhogq_pl  + TIME_DTL * frhogq_pl
       endif

       call prgvar_set( rhog,   rhog_pl,   & !--- [IN]
                        rhogvx, rhogvx_pl, & !--- [IN]
                        rhogvy, rhogvy_pl, & !--- [IN]
                        rhogvz, rhogvz_pl, & !--- [IN]
                        rhogw,  rhogw_pl,  & !--- [IN]
                        rhoge,  rhoge_pl,  & !--- [IN]
                        rhogq,  rhogq_pl,  & !--- [IN]
                        0                  )

    endif

    return
  end subroutine dynstep2

end module mod_dynstep
