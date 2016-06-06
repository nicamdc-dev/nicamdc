!-------------------------------------------------------------------------------
!>
!! Basic state module
!!
!! @par Description
!!         This module provides the subroutines for basic state.
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita) Imported from igdc-4.33
!! @li      2011-08-13 (A.Noda)   Add twp-ice exp.
!!
!<
module mod_bsstate
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_MAXFNAME, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: bsstate_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: rho_bs   (:,:,:)
  real(RP), public, allocatable :: rho_bs_pl(:,:,:)
  real(RP), public, allocatable :: pre_bs   (:,:,:)
  real(RP), public, allocatable :: pre_bs_pl(:,:,:)
  real(RP), public, allocatable :: tem_bs   (:,:,:)
  real(RP), public, allocatable :: tem_bs_pl(:,:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: bsstate_input_ref
  private :: bsstate_output_ref
  private :: bsstate_generate
  private :: set_basicstate
  private :: output_info

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=ADM_NSYS), private :: ref_type = 'NOBASE'  !--- Basic state type
  !                                            = 'NOBASE' : no basic state
  !                                            = 'INPUT'  : input
  !                                            = 'TEM'    : temperature is given.
  !                                            = 'TH'     : potential temperature is given.

  character(len=ADM_MAXFNAME), private :: ref_fname = 'ref.dat'

  real(RP), private, allocatable :: phi_ref(:) ! reference phi
  real(RP), private, allocatable :: rho_ref(:) ! reference density
  real(RP), private, allocatable :: pre_ref(:) ! reference pressure
  real(RP), private, allocatable :: tem_ref(:) ! reference temperature
  real(RP), private, allocatable :: qv_ref (:) ! water vapor
  real(RP), private, allocatable :: th_ref (:) ! reference potential temperature

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine bsstate_setup
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_kall,      &
       ADM_gall_pl,   &
       ADM_gall
    implicit none

    namelist / BSSTATEPARAM / &
       ref_type,  &
       ref_fname

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[basic state]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=BSSTATEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** BSSTATEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist BSSTATEPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist BSSTATEPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=BSSTATEPARAM)

    !--- allocation of reference variables
    allocate( phi_ref(ADM_kall) )
    allocate( rho_ref(ADM_kall) )
    allocate( pre_ref(ADM_kall) )
    allocate( tem_ref(ADM_kall) )
    allocate( th_ref (ADM_kall) )
    allocate( qv_ref (ADM_kall) )
    phi_ref(:) = 0.0_RP
    rho_ref(:) = 0.0_RP
    pre_ref(:) = 0.0_RP
    tem_ref(:) = 0.0_RP
    th_ref (:) = 0.0_RP
    qv_ref (:) = 0.0_RP

    allocate( rho_bs   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( rho_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    rho_bs   (:,:,:) = 0.0_RP
    rho_bs_pl(:,:,:) = 0.0_RP

    allocate( pre_bs   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( pre_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    pre_bs   (:,:,:) = 0.0_RP
    pre_bs_pl(:,:,:) = 0.0_RP

    allocate( tem_bs   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( tem_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    tem_bs   (:,:,:) = 0.0_RP
    tem_bs_pl(:,:,:) = 0.0_RP

    if    ( ref_type == 'INPUT' ) then
       call bsstate_input_ref(ref_fname)
    elseif( ref_type == 'INIT' ) then
       call bsstate_generate
       call bsstate_output_ref(ref_fname)
    endif

    !--- calculation of basic state
    call set_basicstate

    !--- output the information
    call output_info

    return
  end subroutine bsstate_setup

  !-----------------------------------------------------------------------------
  subroutine bsstate_input_ref( basename )
    use mod_misc, only: &
       MISC_get_available_fid
    use mod_adm, only: &
       ADM_kall
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV, &
       Rdry  => CNST_RAIR,  &
       CPdry => CNST_CP,    &
       CVdry => CNST_CV,    &
       Rvap  => CNST_RVAP,  &
       PRE00 => CNST_PRE00
    use mod_grd, only: &
       GRD_gz
    implicit none

    Character(*), Intent(in) :: basename

    real(DP) :: pre_ref_DP(ADM_kall)
    real(DP) :: tem_ref_DP(ADM_kall)
    real(DP) :: qv_ref_DP (ADM_kall)

    real(RP) :: kappa
    integer  :: fid
    integer  :: k
    !---------------------------------------------------------------------------

    kappa = Rdry / CPdry

    fid = MISC_get_available_fid()
    open( unit   = fid,            &
          file   = trim(basename), &
          status = 'old',          &
          form   = 'unformatted'   )
       read(fid) pre_ref_DP(:)
       read(fid) tem_ref_DP(:)
       read(fid) qv_ref_DP (:)
    close(fid)

    pre_ref(:) = real(pre_ref_DP(:),kind=RP)
    tem_ref(:) = real(tem_ref_DP(:),kind=RP)
    qv_ref (:) = real(qv_ref_DP (:),kind=RP)

    !--- additional reference state.
    do k = 1, ADM_kall
       th_ref (k) = tem_ref(k) * ( PRE00 / pre_ref(k) )**kappa
       phi_ref(k) = GRAV * GRD_gz(k)
       rho_ref(k) = pre_ref(k) / tem_ref(k) / ( ( 1.0_RP - qv_ref(k) ) * Rdry &
                                              + (          qv_ref(k) ) * Rvap )
    enddo

    return
  end subroutine bsstate_input_ref

  !-----------------------------------------------------------------------------
  subroutine bsstate_output_ref( basename )
    use mod_misc, only: &
       MISC_get_available_fid
    implicit none

    Character(*), Intent(in) :: basename

    integer :: fid
    !---------------------------------------------------------------------------

    fid = MISC_get_available_fid()
    open( unit   = fid,            &
          file   = trim(basename), &
          status = 'replace',      &
          form   = 'unformatted'   )
       write(fid) pre_ref(:)
       write(fid) tem_ref(:)
       write(fid) qv_ref (:)
    close(fid)

    return
  end subroutine bsstate_output_ref

  !-----------------------------------------------------------------------------
  subroutine bsstate_generate
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV, &
       Rdry  => CNST_RAIR,  &
       CPdry => CNST_CP,    &
       CVdry => CNST_CV,    &
       Rvap  => CNST_RVAP,  &
       PRE00 => CNST_PRE00
    use mod_grd, only: &
       GRD_gz
    use mod_gtl, only: &
       GTL_global_mean_eachlayer
    use mod_runconf, only: &
       TRC_vmax, &
       I_QV
    use mod_prgvar, only: &
       prgvar_get_withdiag
    implicit none

    real(RP) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogq    (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(RP) :: rho      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: pre      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tem      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vx       (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vy       (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vz       (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: w        (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q        (ADM_gall,   ADM_kall,ADM_lall,   TRC_vmax)
    real(RP) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(RP) :: kappa

    integer :: k
    !---------------------------------------------------------------------------

    kappa = Rdry / CPdry

    call prgvar_get_withdiag( rhog,   rhog_pl,   & ! [OUT]
                              rhogvx, rhogvx_pl, & ! [OUT]
                              rhogvy, rhogvy_pl, & ! [OUT]
                              rhogvz, rhogvz_pl, & ! [OUT]
                              rhogw,  rhogw_pl,  & ! [OUT]
                              rhoge,  rhoge_pl,  & ! [OUT]
                              rhogq,  rhogq_pl,  & ! [OUT]
                              rho,    rho_pl,    & ! [OUT]
                              pre,    pre_pl,    & ! [OUT]
                              tem,    tem_pl,    & ! [OUT]
                              vx,     vx_pl,     & ! [OUT]
                              vy,     vy_pl,     & ! [OUT]
                              vz,     vz_pl,     & ! [OUT]
                              w,      w_pl,      & ! [OUT]
                              q,      q_pl       ) ! [OUT]

    call GTL_global_mean_eachlayer( pre(:,:,:),      pre_pl(:,:,:)     , pre_ref(:) )
    call GTL_global_mean_eachlayer( tem(:,:,:),      tem_pl(:,:,:)     , tem_ref(:) )
    call GTL_global_mean_eachlayer( q  (:,:,:,I_QV), q_pl  (:,:,:,I_QV), qv_ref (:) )

    !--- additional reference state.
    do k = 1, ADM_kall
       th_ref (k) = tem_ref(k) * ( PRE00 / pre_ref(k) )**kappa
       phi_ref(k) = GRAV * GRD_gz(k)
       rho_ref(k) = pre_ref(k) / tem_ref(k) / ( ( 1.0_RP - qv_ref(k) ) * Rdry &
                                              + (          qv_ref(k) ) * Rvap )
    enddo

    return
  end subroutine bsstate_generate

  !-----------------------------------------------------------------------------
  !> generation of basic state from reference state
  subroutine set_basicstate
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV
    use mod_adm, only: &
       ADM_prc_pl,  &
       ADM_prc_me,  &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_grd, only: &
       GRD_Z,     &
       GRD_vz,    &
       GRD_vz_pl
    use mod_vintrpl, only: &
       VINTRPL_zstar_level
    use mod_bndcnd, only: &
       bndcnd_thermo
    use mod_thrmdyn, only: &
       thrmdyn_th,  &
       thrmdyn_rho, &
       thrmdyn_qd
    use mod_runconf, only: &
       TRC_VMAX, &
       I_QV
    implicit none

    real(RP) :: phi     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: phi_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: qd      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qd_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q       (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP) :: q_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(RP) :: th_bs   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: th_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: qv_bs   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qv_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: k, l
    !---------------------------------------------------------------------------

    !--- calculation of geo-potential
    phi(:,:,:) = GRAV * GRD_vz(:,:,:,GRD_Z)
    if ( ADM_prc_me == ADM_prc_pl ) then
       phi_pl(:,:,:) = GRAV * GRD_vz_pl(:,:,:,GRD_Z)
    endif

    if ( ref_type == 'NOBASE' ) return

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       pre_bs(:,k,l) = pre_ref(k)
       tem_bs(:,k,l) = tem_ref(k)
       th_bs (:,k,l) = th_ref (k)
       qv_bs (:,k,l) = qv_ref (k)
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          pre_bs_pl(:,k,l) = pre_ref(k)
          tem_bs_pl(:,k,l) = tem_ref(k)
          th_bs_pl (:,k,l) = th_ref (k)
          qv_bs_pl (:,k,l) = qv_ref (k)
       enddo
       enddo
    endif

    !-- from z-level to zstar-level
    call VINTRPL_zstar_level( pre_bs, pre_bs_pl, .false. )
    call VINTRPL_zstar_level( tem_bs, tem_bs_pl, .false. )
    call VINTRPL_zstar_level( qv_bs,  qv_bs_pl,  .false. )

    do l = 1, ADM_lall
       !--- Setting of mass concentration
       !--- Note :: The basic state is "dry" and TKE=0
       q(:,:,l,:)    = 0.0_RP
       q(:,:,l,I_QV) = qv_bs(:,:,l)

       call thrmdyn_qd( ADM_gall,    & ! [IN]
                        ADM_kall,    & ! [IN]
                        q (:,:,l,:), & ! [IN]
                        qd(:,:,l)    ) ! [OUT]

       !--- calculation of density
       call thrmdyn_rho( ADM_gall,        & ! [IN]
                         ADM_kall,        & ! [IN]
                         pre_bs(:,:,l),   & ! [IN]
                         tem_bs(:,:,l),   & ! [IN]
                         qd    (:,:,l),   & ! [IN]
                         q     (:,:,l,:), & ! [IN]
                         rho_bs(:,:,l)    )! [OUT]

       !--- set boundary conditions of basic state
       call bndcnd_thermo( ADM_gall,      & ! [IN]
                           tem_bs(:,:,l), & ! [INOUT]
                           rho_bs(:,:,l), & ! [INOUT]
                           pre_bs(:,:,l), & ! [INOUT]
                           phi   (:,:,l)  ) ! [IN]

       call thrmdyn_th( ADM_gall,      & ! [IN]
                        ADM_kall,      & ! [IN]
                        tem_bs(:,:,l), & ! [IN]
                        pre_bs(:,:,l), & ! [IN]
                        th_bs (:,:,l)  ) ! [OUT]
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          q_pl(:,:,l,:)    = 0.0_RP
          q_pl(:,:,l,I_QV) = qv_bs_pl(:,:,l)

          call thrmdyn_qd( ADM_gall_pl,    & ! [IN]
                           ADM_kall,       & ! [IN]
                           q_pl (:,:,l,:), & ! [IN]
                           qd_pl(:,:,l)    ) ! [OUT]

          call thrmdyn_rho( ADM_gall_pl,        & ! [IN]
                            ADM_kall,           & ! [IN]
                            tem_bs_pl(:,:,l),   & ! [IN]
                            pre_bs_pl(:,:,l),   & ! [IN]
                            qd_pl    (:,:,l),   & ! [IN]
                            q_pl     (:,:,l,:), & ! [IN]
                            rho_bs_pl(:,:,l)    ) ! [OUT]

          call bndcnd_thermo( ADM_gall_pl,      & ! [IN]
                              tem_bs_pl(:,:,l), & ! [INOUT]
                              rho_bs_pl(:,:,l), & ! [INOUT]
                              pre_bs_pl(:,:,l), & ! [INOUT]
                              phi_pl   (:,:,l)  ) ! [IN]

          call thrmdyn_th( ADM_gall_pl,      & ! [IN]
                           ADM_kall,         & ! [IN]
                           tem_bs_pl(:,:,l), & ! [IN]
                           pre_bs_pl(:,:,l), & ! [IN]
                           th_bs_pl (:,:,l)  ) ! [OUT]
       enddo
    endif

    return
  end subroutine set_basicstate

  !-----------------------------------------------------------------------------
  subroutine output_info
    use mod_adm, only: &
       ADM_LOG_FID, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** Basic state information ***'
    write(ADM_LOG_FID,*) '--- Basic state type         : ', trim(ref_type)

    write(ADM_LOG_FID,*) '-------------------------------------------------------'
    write(ADM_LOG_FID,*) 'Level   Density  Pressure     Temp. Pot. Tem.        qv'
    do k = 1, ADM_kall
       write(ADM_LOG_FID,'(I4,F12.4,F10.2,F10.2,F10.2,F10.7)') k, rho_ref(k), pre_ref(k), tem_ref(k),th_ref(k), qv_ref(k)
       if (  k == ADM_kmin-1 ) write(ADM_LOG_FID,*) '-------------------------------------------------------'
       if (  k == ADM_kmax   ) write(ADM_LOG_FID,*) '-------------------------------------------------------'
    enddo

    return
  end subroutine output_info

end module mod_bsstate
