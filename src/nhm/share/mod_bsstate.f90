!-------------------------------------------------------------------------------
!> Module basic state
!!
!! @par Description
!!          This module is for the set of basic state for non-hydrostatic model
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_bsstate
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

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: ref_type  = 'NOBASE' !--- Basic state type
                                               ! 'NOBASE' : no basic state
                                               ! 'INPUT'  : input
                                               ! 'TEM'    : temperature is given.
                                               ! 'TH'     : potential temperature is given.

  character(len=H_LONG),  private :: ref_fname      = 'ref.dat'
  character(len=H_LONG),  private :: sounding_fname = ''
  real(RP),               private :: tem_sfc                    ! Surface temperature [K]
  real(RP),               private :: pre_sfc                    ! Surface pressure    [Pa]
  real(RP),               private :: BV_freq                    ! Brunt-Vaisala frequency [1/s]
  real(RP),               private :: lapse_rate                 ! temperature lapse rate [K/m]
  real(RP),               private :: Z_tem                      ! constant value (tem/th) is applied below this altitude [m]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine bsstate_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only: &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       PRE00 => CONST_PRE00, &
       Pstd  => CONST_Pstd,  &
       Tstd  => CONST_Tstd
    use mod_adm, only: &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall_pl, &
       ADM_gall,    &
       ADM_kmax,    &
       ADM_kmin
    implicit none

    namelist / BSSTATEPARAM / &
       ref_type,       &
       ref_fname,      &
       sounding_fname, &
       pre_sfc,        &
       tem_sfc,        &
       BV_freq,        &
       lapse_rate,     &
       Z_tem

    real(RP) :: pre_ref(ADM_kall) ! reference pressure
    real(RP) :: tem_ref(ADM_kall) ! reference temperature
    real(RP) :: qv_ref (ADM_kall) ! water vapor
    real(RP) :: rho_ref(ADM_kall) ! density
    real(RP) :: th_ref (ADM_kall) ! potentical temperature (dry)

    integer  :: k
    integer  :: ierr
    !---------------------------------------------------------------------------

    pre_sfc    = Pstd
    tem_sfc    = Tstd
    BV_freq    = 0.0_RP
    lapse_rate = 0.0_RP
    Z_tem      = 0.0_RP

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[basic state]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=BSSTATEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** BSSTATEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist BSSTATEPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=BSSTATEPARAM)

    !--- allocation of reference variables
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Basic state information ***'
    if( IO_L ) write(IO_FID_LOG,*) '--- Basic state type : ', trim(ref_type)

    if    ( ref_type == 'NOBASE' ) then

       ! do nothing

    elseif( ref_type == 'INPUT' ) then

       call bsstate_input_ref ( ref_fname,  & ! [IN]
                                pre_ref(:), & ! [OUT]
                                tem_ref(:), & ! [OUT]
                                qv_ref (:)  ) ! [OUT]

    else

       call bsstate_generate  ( sounding_fname, & ! [IN]
                                pre_ref(:),     & ! [OUT]
                                tem_ref(:),     & ! [OUT]
                                qv_ref (:)      ) ! [OUT]

       call bsstate_output_ref( ref_fname,  & ! [IN]
                                pre_ref(:), & ! [IN]
                                tem_ref(:), & ! [IN]
                                qv_ref (:)  ) ! [IN]

    endif

    if ( ref_type /= 'NOBASE' ) then
       ! set 3-D basic state
       call set_basicstate( pre_ref(:), & ! [IN]
                            tem_ref(:), & ! [IN]
                            qv_ref (:)  ) ! [IN]

       if( IO_L ) write(IO_FID_LOG,*) '-------------------------------------------------------'
       if( IO_L ) write(IO_FID_LOG,*) 'Level   Density  Pressure     Temp. Pot. Tem.        qv'

       do k = ADM_kall, 1, -1
          th_ref (k) = tem_ref(k) * ( PRE00 / pre_ref(k) )**(Rdry/CPdry)
          rho_ref(k) = pre_ref(k) / tem_ref(k) / ( ( 1.0_RP - qv_ref(k) ) * Rdry &
                                                 + (          qv_ref(k) ) * Rvap )

          if ( k == ADM_kmax ) then
             if( IO_L ) write(IO_FID_LOG,*) '-------------------------------------------------------'
          endif
          if( IO_L ) write(IO_FID_LOG,'(I4,F12.4,3F10.2,F10.7)') k,rho_ref(k),pre_ref(k),tem_ref(k),th_ref(k),qv_ref(k)
          if ( k == ADM_kmin ) then
             if( IO_L ) write(IO_FID_LOG,*) '-------------------------------------------------------'
          endif
       enddo
    endif

    return
  end subroutine bsstate_setup

  !-----------------------------------------------------------------------------
  subroutine bsstate_input_ref( &
       fname,   &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall
    implicit none

    character(len=*), intent(in)  :: fname
    real(RP),         intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP),         intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP),         intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(DP) :: pre_ref_DP(ADM_kall)
    real(DP) :: tem_ref_DP(ADM_kall)
    real(DP) :: qv_ref_DP (ADM_kall)

    integer  :: fid
    !---------------------------------------------------------------------------

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'old'          )

       read(fid) pre_ref_DP(:)
       read(fid) tem_ref_DP(:)
       read(fid) qv_ref_DP (:)

    close(fid)

    pre_ref(:) = real(pre_ref_DP(:),kind=RP)
    tem_ref(:) = real(tem_ref_DP(:),kind=RP)
    qv_ref (:) = real(qv_ref_DP (:),kind=RP)

    return
  end subroutine bsstate_input_ref

  !-----------------------------------------------------------------------------
  subroutine bsstate_output_ref( &
       fname,   &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall
    implicit none

    character(len=*), intent(in)  :: fname
    real(RP),         intent(in)  :: pre_ref(ADM_kall) ! reference pressure
    real(RP),         intent(in)  :: tem_ref(ADM_kall) ! reference temperature
    real(RP),         intent(in)  :: qv_ref (ADM_kall) ! water vapor

    real(DP) :: pre_ref_DP(ADM_kall)
    real(DP) :: tem_ref_DP(ADM_kall)
    real(DP) :: qv_ref_DP (ADM_kall)

    integer  :: fid
    !---------------------------------------------------------------------------

    pre_ref_DP(:) = real(pre_ref(:),kind=DP)
    tem_ref_DP(:) = real(tem_ref(:),kind=DP)
    qv_ref_DP (:) = real(qv_ref (:),kind=DP)

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new'          )

       write(fid) pre_ref_DP(:)
       write(fid) tem_ref_DP(:)
       write(fid) qv_ref_DP (:)

    close(fid)

    return
  end subroutine bsstate_output_ref

  !-----------------------------------------------------------------------------
  subroutine bsstate_generate( &
       sounding_fname, &
       pre_ref,        &
       tem_ref,        &
       qv_ref          )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_gtl, only: &
       GTL_global_mean_eachlayer
    use mod_runconf, only: &
       TRC_vmax, &
       I_QV
    use mod_prgvar, only: &
       prgvar_get_withdiag
    implicit none

    character(len=*), intent(in)  :: sounding_fname
    real(RP),         intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP),         intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP),         intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(RP) :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhoge    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogq    (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(RP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: pre      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tem      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vx       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vy       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vz       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: w        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q        (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    !---------------------------------------------------------------------------

    if ( ref_type == 'INIT' ) then

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

    elseif( ref_type == 'TEM' ) then ! based on the reference temperature

       call refgen_tem( pre_ref(:), & ! [OUT]
                        tem_ref(:), & ! [OUT]
                        qv_ref (:)  ) ! [OUT]

    elseif( ref_type == 'RHO' ) then ! based on the reference density

       call refgen_rho( pre_ref(:), & ! [OUT]
                        tem_ref(:), & ! [OUT]
                        qv_ref (:)  ) ! [OUT]

    elseif( ref_type == 'TH' ) then

       call refgen_th ( pre_ref(:), & ! [OUT]
                        tem_ref(:), & ! [OUT]
                        qv_ref (:)  ) ! [OUT]

    elseif( ref_type == 'TH-SP' ) then

       call refgen_th2( pre_ref(:), & ! [OUT]
                        tem_ref(:), & ! [OUT]
                        qv_ref (:)  ) ! [OUT]

    elseif( ref_type == 'OOYAMA' ) then

       call refgen_ooyama( sounding_fname, & ! [IN]
                           pre_ref(:),     & ! [OUT]
                           tem_ref(:),     & ! [OUT]
                           qv_ref (:)      ) ! [OUT]

    elseif(      ref_type == 'GCSS_CASE1' &
            .OR. ref_type == 'TWP-ICE'    ) then

       call refgen_gcss( sounding_fname, & ! [IN]
                         pre_ref(:),     & ! [OUT]
                         tem_ref(:),     & ! [OUT]
                         qv_ref (:)      ) ! [OUT]

    endif

    return
  end subroutine bsstate_generate

  !-----------------------------------------------------------------------------
  !> generation of basic state from reference state
  subroutine set_basicstate( &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_const, only: &
       Rdry => CONST_Rdry, &
       Rvap => CONST_Rvap
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_vmtr, only: &
       VMTR_PHI,   &
       VMTR_PHI_pl
    use mod_vintrpl, only: &
       VINTRPL_Z2Xi
    use mod_bndcnd, only: &
       BNDCND_thermo
    implicit none

    real(RP), intent(in) :: pre_ref(ADM_kall) ! reference pressure
    real(RP), intent(in) :: tem_ref(ADM_kall) ! reference temperature
    real(RP), intent(in) :: qv_ref (ADM_kall) ! water vapor

    real(RP) :: qv_bs   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qv_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer  :: k, l
    !---------------------------------------------------------------------------

    if ( ref_type == 'NOBASE' ) return

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       pre_bs(:,k,l) = pre_ref(k)
       tem_bs(:,k,l) = tem_ref(k)
       qv_bs (:,k,l) = qv_ref (k)
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          pre_bs_pl(:,k,l) = pre_ref(k)
          tem_bs_pl(:,k,l) = tem_ref(k)
          qv_bs_pl (:,k,l) = qv_ref (k)
       enddo
       enddo
    endif

    !-- from z-level to zstar-level
    call VINTRPL_Z2Xi( pre_bs(:,:,:), pre_bs_pl(:,:,:) )
    call VINTRPL_Z2Xi( tem_bs(:,:,:), tem_bs_pl(:,:,:) )
    call VINTRPL_Z2Xi( qv_bs (:,:,:), qv_bs_pl (:,:,:) )

    !--- set boundary conditions of basic state
    do l = 1, ADM_lall
       rho_bs(:,:,l) = pre_bs(:,:,l) / tem_bs(:,:,l) / ( ( 1.0_RP-qv_bs(:,:,l) ) * Rdry &
                                                       + (        qv_bs(:,:,l) ) * Rvap )

       call BNDCND_thermo( ADM_gall,        & ! [IN]
                           tem_bs  (:,:,l), & ! [INOUT]
                           rho_bs  (:,:,l), & ! [INOUT]
                           pre_bs  (:,:,l), & ! [INOUT]
                           VMTR_PHI(:,:,l)  ) ! [IN]
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          rho_bs_pl(:,:,l) = pre_bs_pl(:,:,l) / tem_bs_pl(:,:,l) / ( ( 1.0_RP-qv_bs_pl(:,:,l) ) * Rdry &
                                                                   + (        qv_bs_pl(:,:,l) ) * Rvap )

          call BNDCND_thermo( ADM_gall_pl,        & ! [IN]
                              tem_bs_pl  (:,:,l), & ! [INOUT]
                              rho_bs_pl  (:,:,l), & ! [INOUT]
                              pre_bs_pl  (:,:,l), & ! [INOUT]
                              VMTR_PHI_pl(:,:,l)  ) ! [IN]
       enddo
    endif

    return
  end subroutine set_basicstate

  !-----------------------------------------------------------------------------
  subroutine refgen_tem( &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall, &
       ADM_kmin, &
       ADM_kmax
    use mod_const, only: &
       GRAV => CONST_GRAV, &
       Rdry => CONST_Rdry
    use mod_grd, only: &
       GRD_gz,    &
       GRD_dgz,   &
       GRD_afact, &
       GRD_bfact
    implicit none

    real(RP), intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP), intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP), intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(RP) :: rho_ref(ADM_kall)
    real(RP) :: rho_sfc

    integer  :: k
    !---------------------------------------------------------------------------

    do k = 1, ADM_kall
       if ( GRD_gz(k) <= Z_TEM ) then
          tem_ref(k) = tem_sfc - lapse_rate * GRD_gz(k)
       else
          tem_ref(k) = tem_sfc - lapse_rate * Z_TEM
       endif
    enddo

    rho_sfc = pre_sfc / ( Rdry * tem_sfc )
    pre_ref(ADM_kmin-1) = pre_sfc    + 0.5_RP * rho_sfc * GRAV * GRD_dgz(ADM_kmin-1)
    pre_ref(ADM_kmin)   = pre_ref(ADM_kmin-1) - rho_sfc * GRAV * GRD_dgz(ADM_kmin-1)
    rho_ref(ADM_kmin)   = pre_ref(ADM_kmin) / ( Rdry * tem_ref(ADM_kmin) )

    do k = ADM_kmin+1, ADM_kmax+1
       rho_ref(k) = ( pre_ref(k-1) - rho_ref(k-1) * GRD_bfact(k) * GRAV * GRD_dgz(k-1) ) &
                  / (           Rdry * tem_ref(k) + GRD_afact(k) * GRAV * GRD_dgz(k-1) )
       pre_ref(k) = rho_ref(k) * Rdry * tem_ref(k)
    enddo

    qv_ref(:) = 0.0_RP

    return
  end subroutine refgen_tem

  !-----------------------------------------------------------------------------
  subroutine refgen_rho( &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall, &
       ADM_kmin, &
       ADM_kmax
    use mod_const, only: &
       EPS  => CONST_EPS,  &
       GRAV => CONST_GRAV, &
       Rdry => CONST_Rdry
    use mod_grd, only: &
       GRD_gz,    &
       GRD_dgz,   &
       GRD_afact, &
       GRD_bfact
    implicit none

    real(RP), intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP), intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP), intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(RP) :: rho_ref(ADM_kall)

    real(RP) :: total_mass0, total_mass, mass_diff_ratio
    real(RP) :: pre_s, rho_s

    integer  :: k
    !---------------------------------------------------------------------------

    total_mass0 = pre_sfc / GRAV
    pre_s       = pre_sfc
    do
       do k = ADM_kmin-1, ADM_kmax+1
          rho_ref(k) = pre_s / ( Rdry * tem_sfc ) / exp( GRAV * GRD_gz(k) / ( Rdry * tem_sfc ) )
       enddo

       total_mass = 0.0_RP
       do k = ADM_kmin, ADM_kmax
          total_mass = total_mass + rho_ref(k) * GRD_dgz(k)
       enddo

       mass_diff_ratio = total_mass0 / total_mass

       if ( abs( mass_diff_ratio-1.0_RP ) < EPS ) then
          exit
       else
          pre_s = pre_s * mass_diff_ratio
       endif
    enddo

    rho_s = pre_s / ( Rdry * tem_sfc )

    pre_ref(ADM_kmin-1) = pre_s      + 0.5_RP * rho_s * GRAV * GRD_dgz(ADM_kmin-1)
    pre_ref(ADM_kmin)   = pre_ref(ADM_kmin-1) - rho_s * GRAV * GRD_dgz(ADM_kmin-1)

    do k = ADM_kmin+1, ADM_kmax+1
       pre_ref(k) = pre_ref(k-1) - ( GRD_afact(k) * rho_ref(k)   &
                                   + GRD_bfact(k) * rho_ref(k-1) ) * GRAV * GRD_dgz(k-1)
    enddo

    tem_ref(:) = pre_ref(:) / ( Rdry * rho_ref(:) )
    qv_ref (:) = 0.0_RP

    return
  end subroutine refgen_rho

  !-----------------------------------------------------------------------------
  subroutine refgen_th( &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall, &
       ADM_kmin, &
       ADM_kmax
    use mod_const, only: &
       GRAV  => CONST_GRAV,  &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       PRE00 => CONST_PRE00
    use mod_grd, only: &
       GRD_gz,    &
       GRD_gzh,   &
       GRD_dgz,   &
       GRD_afact, &
       GRD_bfact
    implicit none

    real(RP), intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP), intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP), intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(RP) :: rho_ref(ADM_kall)
    real(RP) :: th_ref (ADM_kall)
    real(RP) :: rho_sfc
    real(RP) :: dpre_ref_k

    integer  :: k
    !---------------------------------------------------------------------------

    ! calculation of reference potential temperature
    th_ref(ADM_kmin-1) = tem_sfc / exp( BV_freq**2 / GRAV * ( GRD_gzh(ADM_kmin) - GRD_gz (ADM_kmin-1) ) )
    th_ref(ADM_kmin)   = tem_sfc * exp( BV_freq**2 / GRAV * ( GRD_gz (ADM_kmin) - GRD_gzh(ADM_kmin)   ) )

    do k = ADM_kmin+1, ADM_kmax+1
       if ( GRD_gz(k) <= Z_tem ) then
          th_ref(k) = th_ref(k-1) * exp( BV_freq**2 / GRAV * GRD_dgz(k-1) )
       else
          th_ref(k) = th_ref(k-1)
       endif
    enddo

    !--- calculation of density at the surface
    rho_sfc = pre_sfc / ( Rdry * tem_sfc )

    pre_ref(ADM_kmin-1) = pre_sfc    + 0.5_RP * rho_sfc * GRAV * GRD_dgz(ADM_kmin-1)
    tem_ref(ADM_kmin-1) = th_ref (ADM_kmin-1) / ( PRE00 / pre_ref(ADM_kmin-1) )**(Rdry/CPdry)
    rho_ref(ADM_kmin-1) = pre_ref(ADM_kmin-1) / ( Rdry  * tem_ref(ADM_kmin-1) )

    pre_ref(ADM_kmin)   = pre_ref(ADM_kmin-1) - rho_sfc * GRAV * GRD_dgz(ADM_kmin-1)
    tem_ref(ADM_kmin)   = th_ref (ADM_kmin)   / ( PRE00 / pre_ref(ADM_kmin)   )**(Rdry/CPdry)
    rho_ref(ADM_kmin)   = pre_ref(ADM_kmin)   / ( Rdry  * tem_ref(ADM_kmin)   )

    !--- Reference pressure and density at k level
    !---    In this caluculation, the hydrostatic balance equation
    !---    ( dP/dz=-RHO dPHI/dz ) is applied at the half integer
    !---    level ( k-1/2 ). RHO is obtained by extrapolation from
    !---    the values at lower level.

    ! fist guess
    do k = ADM_kmin+1, ADM_kmax+1
       pre_ref(k) = pre_ref(k-1) &
                  - GRAV * GRD_dgz(k-1) &
                  * ( rho_ref(k-1) + 0.5_RP * ( rho_ref(k-1) - rho_ref(k-2) ) * GRD_dgz(k-1) / GRD_dgz(k-2) )
       tem_ref(k) = th_ref (k) / ( PRE00 / pre_ref(k) )**(Rdry/CPdry)
       rho_ref(k) = pre_ref(k) / ( Rdry  * tem_ref(k) )
    enddo

    ! hydrostatic balance adjustment
    do k = ADM_kmin+1, ADM_kmax+1
       do
          tem_ref(k) = th_ref(k) / ( PRE00 / pre_ref(k) )**(Rdry/CPdry)
          rho_ref(k) = ( pre_ref(k-1) - rho_ref(k-1) * GRD_bfact(k) * GRAV * GRD_dgz(k-1) ) &
                     / (           Rdry * tem_ref(k) + GRD_afact(k) * GRAV * GRD_dgz(k-1) )

          dpre_ref_k = rho_ref(k) * Rdry * tem_ref(k) - pre_ref(k)
          pre_ref(k) = pre_ref(k) + dpre_ref_k

          if( abs(dpre_ref_k) < 1.E-10_RP ) exit
       enddo
    enddo

    qv_ref(:) = 0.0_RP

    return
  end subroutine refgen_th

  !-----------------------------------------------------------------------------
  subroutine refgen_th2( &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall, &
       ADM_kmin, &
       ADM_kmax
    use mod_const, only: &
       GRAV  => CONST_GRAV,  &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       PRE00 => CONST_PRE00
    use mod_grd, only: &
       GRD_gz,    &
       GRD_gzh,   &
       GRD_dgz,   &
       GRD_afact, &
       GRD_bfact
    implicit none

    real(RP), intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP), intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP), intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(RP) :: rho_ref(ADM_kall)
    real(RP) :: th_ref (ADM_kall)
    real(RP) :: rho_sfc
    real(RP) :: dpre_ref_k

    integer  :: k
    !---------------------------------------------------------------------------

    ! calculation of reference potential temperature
    th_ref(ADM_kmin-1) = tem_sfc / exp( BV_freq**2 / GRAV * ( GRD_gzh(ADM_kmin) - GRD_gz (ADM_kmin-1) ) )
    th_ref(ADM_kmin)   = tem_sfc * exp( BV_freq**2 / GRAV * ( GRD_gz (ADM_kmin) - GRD_gzh(ADM_kmin)   ) )

    do k = ADM_kmin+1, ADM_kmax+1
       if ( GRD_gz(k) <= 10000.0_RP ) then
          th_ref(k) = tem_sfc
       else
          th_ref(k) = tem_sfc + 10.0_RP / 1000.0_RP * ( GRD_gz(k) - 10000.0_RP )
       endif
    enddo

    !--- calculation of density at the surface
    rho_sfc = pre_sfc / ( Rdry * tem_sfc )

    pre_ref(ADM_kmin-1) = pre_sfc    + 0.5_RP * rho_sfc * GRAV * GRD_dgz(ADM_kmin-1)
    tem_ref(ADM_kmin-1) = th_ref (ADM_kmin-1) / ( PRE00 / pre_ref(ADM_kmin-1) )**(Rdry/CPdry)
    rho_ref(ADM_kmin-1) = pre_ref(ADM_kmin-1) / ( Rdry  * tem_ref(ADM_kmin-1) )

    pre_ref(ADM_kmin)   = pre_ref(ADM_kmin-1) - rho_sfc * GRAV * GRD_dgz(ADM_kmin-1)
    tem_ref(ADM_kmin)   = th_ref (ADM_kmin)   / ( PRE00 / pre_ref(ADM_kmin)   )**(Rdry/CPdry)
    rho_ref(ADM_kmin)   = pre_ref(ADM_kmin)   / ( Rdry  * tem_ref(ADM_kmin)   )

    !--- Reference pressure and density at k level
    !---    In this caluculation, the hydrostatic balance equation
    !---    ( dP/dz=-RHO dPHI/dz ) is applied at the half integer
    !---    level ( k-1/2 ). RHO is obtained by extrapolation from
    !---    the values at lower level.

    ! fist guess
    do k = ADM_kmin+1, ADM_kmax+1
       pre_ref(k) = pre_ref(k-1) &
                  - GRAV * GRD_dgz(k-1) &
                  * ( rho_ref(k-1) + 0.5_RP * ( rho_ref(k-1) - rho_ref(k-2) ) * GRD_dgz(k-1) / GRD_dgz(k-2) )
       tem_ref(k) = th_ref (k) / ( PRE00 / pre_ref(k) )**(Rdry/CPdry)
       rho_ref(k) = pre_ref(k) / ( Rdry  * tem_ref(k) )
    enddo

    ! hydrostatic balance adjustment
    do k = ADM_kmin+1, ADM_kmax+1
       do
          tem_ref(k) = th_ref(k) / ( PRE00 / pre_ref(k) )**(Rdry/CPdry)
          rho_ref(k) = ( pre_ref(k-1) - rho_ref(k-1) * GRD_bfact(k) * GRAV * GRD_dgz(k-1) ) &
                     / (           Rdry * tem_ref(k) + GRD_afact(k) * GRAV * GRD_dgz(k-1) )

          dpre_ref_k = rho_ref(k) * Rdry * tem_ref(k) - pre_ref(k)
          pre_ref(k) = pre_ref(k) + dpre_ref_k

          if( abs(dpre_ref_k) < 1.E-10_RP ) exit
       enddo
    enddo

    qv_ref(:) = 0.0_RP

    return
  end subroutine refgen_th2

  !-----------------------------------------------------------------------------
  subroutine refgen_ooyama( &
       fname,   &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall, &
       ADM_kmin, &
       ADM_kmax
    use mod_const, only: &
       Rdry => CONST_Rdry, &
       Rvap => CONST_Rvap
    use mod_grd, only: &
       GRD_gz
    implicit none

    character(len=*), intent(in)  :: fname
    real(RP),         intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP),         intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP),         intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(DP), allocatable :: z_s  (:)
    real(DP), allocatable :: rho_s(:)
    real(DP), allocatable :: tem_s(:)
    real(DP), allocatable :: qv_s (:)
    integer               :: kmax_s

    real(RP) :: rho_ref(ADM_kall)
    real(DP) :: gz     (ADM_kall)

    integer  :: fid
    integer  :: k, kk, kp

    real(DP) :: lag_intpl
    real(DP) :: z, z1, p1, z2, p2, z3, p3
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ((z-z2)*(z-z3))/((z1-z2)*(z1-z3))*p1 &
                                   + ((z-z1)*(z-z3))/((z2-z1)*(z2-z3))*p2 &
                                   + ((z-z1)*(z-z2))/((z3-z1)*(z3-z2))*p3
    !---------------------------------------------------------------------------

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'old'          )

       read(fid) kmax_s

       allocate( z_s  (kmax_s) )
       allocate( rho_s(kmax_s) )
       allocate( tem_s(kmax_s) )
       allocate( qv_s (kmax_s) )

       read(fid) z_s  (:)
       read(fid) rho_s(:)
       read(fid) tem_s(:)
       read(fid) qv_s (:)
    close(fid)

    gz(:) = real(GRD_gz(:),kind=DP)

    ! interpolation
    do k  = ADM_kmin, ADM_kmax
    do kk = 1, kmax_s-1
       if (       z_s(kk)   <= GRD_gz(k) &
            .AND. z_s(kk+1) >= GRD_gz(k) ) then

          kp = kk
          if( kk == 1 ) kp = 2

          rho_ref(k) = real(lag_intpl(gz(k), z_s(kp+1), rho_s(kp+1), &
                                             z_s(kp  ), rho_s(kp  ), &
                                             z_s(kp-1), rho_s(kp-1)  ),kind=RP)
          tem_ref(k) = real(lag_intpl(gz(k), z_s(kp+1), tem_s(kp+1), &
                                             z_s(kp  ), tem_s(kp  ), &
                                             z_s(kp-1), tem_s(kp-1)  ),kind=RP)
          qv_ref (k) = real(lag_intpl(gz(k), z_s(kp+1), qv_s (kp+1), &
                                             z_s(kp  ), qv_s (kp  ), &
                                             z_s(kp-1), qv_s (kp-1)  ),kind=RP)

          exit
       endif
    enddo
    enddo

    k  = ADM_kmin-1
    kp = 2
    rho_ref(k) = real(lag_intpl(gz(k), z_s(kp+1), rho_s(kp+1), &
                                       z_s(kp  ), rho_s(kp  ), &
                                       z_s(kp-1), rho_s(kp-1)  ),kind=RP)
    tem_ref(k) = real(lag_intpl(gz(k), z_s(kp+1), tem_s(kp+1), &
                                       z_s(kp  ), tem_s(kp  ), &
                                       z_s(kp-1), tem_s(kp-1)  ),kind=RP)
    qv_ref (k) = real(lag_intpl(gz(k), z_s(kp+1), qv_s (kp+1), &
                                       z_s(kp  ), qv_s (kp  ), &
                                       z_s(kp-1), qv_s (kp-1)  ),kind=RP)

    k  = ADM_kmax+1
    kp = kmax_s-1
    rho_ref(k) = real(lag_intpl(gz(k), z_s(kp+1), rho_s(kp+1), &
                                       z_s(kp  ), rho_s(kp  ), &
                                       z_s(kp-1), rho_s(kp-1)  ),kind=RP)
    tem_ref(k) = real(lag_intpl(gz(k), z_s(kp+1), tem_s(kp+1), &
                                       z_s(kp  ), tem_s(kp  ), &
                                       z_s(kp-1), tem_s(kp-1)  ),kind=RP)
    qv_ref (k) = real(lag_intpl(gz(k), z_s(kp+1), qv_s (kp+1), &
                                       z_s(kp  ), qv_s (kp  ), &
                                       z_s(kp-1), qv_s (kp-1)  ),kind=RP)

    deallocate(rho_s)
    deallocate(tem_s)
    deallocate(qv_s )

    pre_ref(:) = rho_ref(:) * tem_ref(:) * ( ( 1.0_RP - qv_ref(:) ) * Rdry &
                                           + (          qv_ref(:) ) * Rvap )

    return
  end subroutine refgen_ooyama

  !-----------------------------------------------------------------------------
  subroutine refgen_gcss( &
       fname,   &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall
    use mod_const, only: &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       PRE00 => CONST_PRE00
    implicit none

    character(len=*), intent(in)  :: fname
    real(RP),         intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP),         intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP),         intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(DP) :: th_ref_DP (ADM_kall)
    real(DP) :: pre_ref_DP(ADM_kall)
    real(DP) :: qv_ref_DP (ADM_kall)

    integer  :: fid
    !---------------------------------------------------------------------------

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'old'          )

       read(fid) th_ref_DP (:)
       read(fid) pre_ref_DP(:)
       read(fid) qv_ref_DP (:)

    close(fid)

    pre_ref(:) = real(pre_ref_DP(:),kind=RP)
    qv_ref (:) = real(qv_ref_DP (:),kind=RP)

    tem_ref(:) = real( th_ref_DP(:) / ( PRE00 / pre_ref_DP(:) )**(Rdry/CPdry), kind=RP )

    return
  end subroutine refgen_gcss

end module mod_bsstate
