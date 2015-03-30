!-------------------------------------------------------------------------------
!>
!! Energy/mass budget monitoring module
!!
!! @par Description
!!         This module is for monitoring the energy/mass budget
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2009-07-10 (H.Tomita)   [NEW]
!! @li      2009-07-28 (H.Tomita)   Bug fix
!!
!<
!-------------------------------------------------------------------------------
module mod_embudget
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
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
  public :: embudget_setup
  public :: embudget_monitor

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
  logical, private :: MNT_ON   = .false.
  integer, private :: MNT_INTV = 1
  integer, private :: MNT_m_fid
  integer, private :: MNT_e_fid

  real(RP), private, allocatable :: evap0     (:,:,:)
  real(RP), private, allocatable :: evap0_pl  (:,:,:)
  real(RP), private, allocatable :: precip0   (:,:,:)
  real(RP), private, allocatable :: precip0_pl(:,:,:)

  real(RP), private, allocatable :: sfcrad0          (:,:,:)
  real(RP), private, allocatable :: sfcrad0_pl       (:,:,:)
  real(RP), private, allocatable :: toarad0          (:,:,:)
  real(RP), private, allocatable :: toarad0_pl       (:,:,:)
  real(RP), private, allocatable :: evap_energy0     (:,:,:)
  real(RP), private, allocatable :: evap_energy0_pl  (:,:,:)
  real(RP), private, allocatable :: precip_energy0   (:,:,:)
  real(RP), private, allocatable :: precip_energy0_pl(:,:,:)
  real(RP), private, allocatable :: sh_flux_sfc0     (:,:,:)
  real(RP), private, allocatable :: sh_flux_sfc0_pl  (:,:,:)
  real(RP), private, allocatable :: lh_flux_sfc0     (:,:,:)
  real(RP), private, allocatable :: lh_flux_sfc0_pl  (:,:,:)

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine embudget_setup
    use mod_misc, only: &
       MISC_get_available_fid
    use mod_adm, only: &
       ADM_LOG_FID,        &
       ADM_CTL_FID,        &
       ADM_proc_stop,      &
       ADM_prc_me,         &
       ADM_prc_run_master, &
       ADM_gall,           &
       ADM_gall_pl,        &
       ADM_KNONE,          &
       ADM_lall,           &
       ADM_lall_pl
    implicit none

    integer :: ierr

    namelist / EMBUDGETPARAM / &
         MNT_INTV, &
         MNT_ON

    integer :: k0
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[embudget]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=EMBUDGETPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** EMBUDGETPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist EMBUDGETPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist EMBUDGETPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=EMBUDGETPARAM)

    if(.not.MNT_ON) return

    ! open budget.info file
    if ( ADM_prc_me == ADM_prc_run_master ) then
       MNT_m_fid  = MISC_get_available_fid()
       open( unit   = MNT_m_fid,          &
             file   = 'MASS_BUDGET.info', &
             form   = 'formatted',        &
             status = 'unknown'           )

          write(MNT_m_fid,'(A6)', ADVANCE='NO') '# STEP'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'Dry air[Kg]'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'Vapor[Kg]'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'Liquid water[Kg]'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'Ice water[Kg]'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'Total water[Kg]'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'diff.of Total water[Kg]'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'precipitaion[Kg]'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'evaporation[Kg]'
          write(MNT_m_fid,'(A22)',ADVANCE='NO') 'evap - precip[Kg]'
          write(MNT_m_fid,*)

       MNT_e_fid = MISC_get_available_fid()
       open( unit   = MNT_e_fid,            &
             file   = 'ENERGY_BUDGET.info', &
             form   = 'formatted',          &
             status = 'unknown'             )

          write(MNT_e_fid,'(A6)', ADVANCE='NO') '# STEP'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Int. E(moist)'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Potential'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Kinematic'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Total energy'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'diff. of tot. energy'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Radiation (SFC)'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Radiation (TOA)'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Sensible heat(SFC)'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Latent heat(SFC)'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Evap energy(SFC)'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Precip energy(SFC)'
          write(MNT_e_fid,'(A22)',ADVANCE='NO') 'Net SFC/TOA energy'
          write(MNT_e_fid,*)
    endif

    k0 = ADM_KNONE

    allocate( evap0     (ADM_gall,   k0,ADM_lall   ) )
    allocate( evap0_pl  (ADM_gall_pl,k0,ADM_lall_pl) )
    allocate( precip0   (ADM_gall,   k0,ADM_lall   ) )
    allocate( precip0_pl(ADM_gall_pl,k0,ADM_lall_pl) )

    allocate( sfcrad0          (ADM_gall,   k0,ADM_lall   ) )
    allocate( sfcrad0_pl       (ADM_gall_pl,k0,ADM_lall_pl) )
    allocate( toarad0          (ADM_gall,   k0,ADM_lall   ) )
    allocate( toarad0_pl       (ADM_gall_pl,k0,ADM_lall_pl) )
    allocate( evap_energy0     (ADM_gall,   k0,ADM_lall   ) )
    allocate( evap_energy0_pl  (ADM_gall_pl,k0,ADM_lall_pl) )
    allocate( precip_energy0   (ADM_gall,   k0,ADM_lall   ) )
    allocate( precip_energy0_pl(ADM_gall_pl,k0,ADM_lall_pl) )
    allocate( sh_flux_sfc0     (ADM_gall,   k0,ADM_lall   ) )
    allocate( sh_flux_sfc0_pl  (ADM_gall_pl,k0,ADM_lall_pl) )
    allocate( lh_flux_sfc0     (ADM_gall,   k0,ADM_lall   ) )
    allocate( lh_flux_sfc0_pl  (ADM_gall_pl,k0,ADM_lall_pl) )

    evap0     (:,:,:) = 0.0_RP
    evap0_pl  (:,:,:) = 0.0_RP
    precip0   (:,:,:) = 0.0_RP
    precip0_pl(:,:,:) = 0.0_RP

    sfcrad0          (:,:,:) = 0.0_RP
    sfcrad0_pl       (:,:,:) = 0.0_RP
    toarad0          (:,:,:) = 0.0_RP
    toarad0_pl       (:,:,:) = 0.0_RP
    evap_energy0     (:,:,:) = 0.0_RP
    evap_energy0_pl  (:,:,:) = 0.0_RP
    precip_energy0   (:,:,:) = 0.0_RP
    precip_energy0_pl(:,:,:) = 0.0_RP
    sh_flux_sfc0     (:,:,:) = 0.0_RP
    sh_flux_sfc0_pl  (:,:,:) = 0.0_RP
    lh_flux_sfc0     (:,:,:) = 0.0_RP
    lh_flux_sfc0_pl  (:,:,:) = 0.0_RP

    call diagnose_energy_mass

    return
  end subroutine embudget_setup

  !-----------------------------------------------------------------------------
  subroutine embudget_monitor
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_lall_pl
    use mod_time, only: &
       TIME_CSTEP, &
       TIME_DTL
    use mod_sfcvar, only: &
       sfcvar_get,      &
       I_PRECIP_TOT,    &
       I_EVAP_SFC,      &
       I_SFCRAD_ENERGY, &
       I_TOARAD_ENERGY, &
       I_EVAP_ENERGY,   &
       I_PRECIP_ENERGY, &
       I_SH_FLUX_SFC,   &
       I_LH_FLUX_SFC
    implicit none

    real(RP) :: evap     (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: evap_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: precip   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: precip_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)

    real(RP) :: sfcrad          (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: sfcrad_pl       (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: toarad          (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: toarad_pl       (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: evap_energy     (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: evap_energy_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: precip_energy   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: precip_energy_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: sh_flux_sfc     (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: sh_flux_sfc_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: lh_flux_sfc     (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: lh_flux_sfc_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    !---------------------------------------------------------------------------

    if( .NOT. MNT_ON ) return

    !--- mass
    call sfcvar_get( evap, evap_pl, vid = I_EVAP_SFC )
    evap0    = evap0    + evap    * TIME_DTL
    evap0_pl = evap0_pl + evap_pl * TIME_DTL

    call sfcvar_get( precip, precip_pl, vid = I_PRECIP_TOT    )
    precip0    = precip0    + precip    * TIME_DTL
    precip0_pl = precip0_pl + precip_pl * TIME_DTL

    !--- energy
    call sfcvar_get( sfcrad, sfcrad_pl, vid = I_SFCRAD_ENERGY )
    sfcrad0    = sfcrad0    + sfcrad    * TIME_DTL
    sfcrad0_pl = sfcrad0_pl + sfcrad_pl * TIME_DTL

    call sfcvar_get( toarad, toarad_pl, vid = I_TOARAD_ENERGY )
    toarad0    = toarad0    + toarad    * TIME_DTL
    toarad0_pl = toarad0_pl + toarad_pl * TIME_DTL

    call sfcvar_get( evap_energy, evap_energy_pl, vid = I_EVAP_ENERGY )
    evap_energy0    = evap_energy0    + evap_energy    * TIME_DTL
    evap_energy0_pl = evap_energy0_pl + evap_energy_pl * TIME_DTL

    call sfcvar_get( precip_energy, precip_energy_pl, vid = I_PRECIP_ENERGY )
    precip_energy0    = precip_energy0    + precip_energy    * TIME_DTL
    precip_energy0_pl = precip_energy0_pl + precip_energy_pl * TIME_DTL

    call sfcvar_get( sh_flux_sfc, sh_flux_sfc_pl, vid = I_SH_FLUX_SFC )
    sh_flux_sfc0    = sh_flux_sfc0    + sh_flux_sfc    * TIME_DTL
    sh_flux_sfc0_pl = sh_flux_sfc0_pl + sh_flux_sfc_pl * TIME_DTL

    call sfcvar_get( lh_flux_sfc, lh_flux_sfc_pl, vid = I_LH_FLUX_SFC )
    lh_flux_sfc0    = lh_flux_sfc0    + lh_flux_sfc    * TIME_DTL
    lh_flux_sfc0_pl = lh_flux_sfc0_pl + lh_flux_sfc_pl * TIME_DTL

    if ( mod(TIME_CSTEP,MNT_INTV) == 0 ) then
       call diagnose_energy_mass
    endif

    return
  end subroutine embudget_monitor

  !-----------------------------------------------------------------------------
  subroutine diagnose_energy_mass
    use mod_adm, only: &
       ADM_prc_me,         &
       ADM_prc_run_master, &
       ADM_prc_pl,         &
       ADM_gall,           &
       ADM_gall_pl,        &
       ADM_kall,           &
       ADM_lall,           &
       ADM_lall_pl
    use mod_cnst, only: &
       CNST_CV
    use mod_time, only: &
       TIME_CSTEP
    use mod_vmtr, only: &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_PHI,        &
       VMTR_PHI_pl
    use mod_gtl, only: &
       GTL_global_sum, &
       GTL_global_sum_srf
    use mod_runconf, only: &
       TRC_vmax, &
       NQW_STR,  &
       NQW_END,  &
       I_QV,     &
       I_QC,     &
       I_QR,     &
       I_QI,     &
       I_QS,     &
       I_QG,     &
       LHV,      &
       LHF,      &
       CVW
    use mod_prgvar, only: &
       prgvar_get_withdiag
    use mod_cnvvar, only: &
       cnvvar_rhogkin
    use mod_thrmdyn, only: &
       THRMDYN_qd
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
    real(RP) :: q        (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    real(RP) :: qd    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: qd_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: tmp   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: tmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rho_sum
    real(RP) :: rhoqd_sum
    real(RP) :: rhoq_sum(1:TRC_vmax)
    real(RP) :: rhoqw_sum
    real(RP) :: rhoqv_sum
    real(RP) :: rhoql_sum
    real(RP) :: rhoqi_sum

    real(RP) :: precip_sum
    real(RP) :: evap_sum

    real(RP) :: rhoetot_sum
    real(RP) :: rhophi_sum
    real(RP) :: rhokin_sum
    real(RP) :: rhoein_sum
    real(RP) :: rhoein_qd_sum
    real(RP) :: rhoein_qw_sum(1:TRC_vmax)

    real(RP) :: sfcrad_sum
    real(RP) :: toarad_sum
    real(RP) :: evap_energy_sum
    real(RP) :: precip_energy_sum
    real(RP) :: sh_flux_sfc_sum
    real(RP) :: lh_flux_sfc_sum

    real(RP), save :: rhoqw_sum_old
    real(RP), save :: rhoetot_sum_old
    logical, save :: iflag = .true.

    integer :: nq
    !---------------------------------------------------------------------------

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

    call THRMDYN_qd( ADM_gall,    & ! [IN]
                     ADM_kall,    & ! [IN]
                     ADM_lall,    & ! [IN]
                     q (:,:,:,:), & ! [IN]
                     qd(:,:,:)    ) ! [OUT]

    if ( ADM_prc_me == ADM_prc_pl ) then
       call THRMDYN_qd( ADM_gall_pl,    & ! [IN]
                        ADM_kall,       & ! [IN]
                        ADM_lall_pl,    & ! [IN]
                        q_pl (:,:,:,:), & ! [IN]
                        qd_pl(:,:,:)    ) ! [OUT]
    endif

    !----- Mass budget

    !--- total mass ( dry + water )
    tmp(:,:,:) = rho(:,:,:)
    if ( ADM_prc_me == ADM_prc_pl ) then
       tmp_pl(:,:,:) = rho_pl(:,:,:)
    end if
    rho_sum = GTL_global_sum( tmp, tmp_pl ) ! [kg/m3] -> [kg]

    !--- total mass (dry air)
    tmp(:,:,:) = rho(:,:,:) * qd(:,:,:)
    if ( ADM_prc_me == ADM_prc_pl ) then
       tmp_pl(:,:,:) = rho_pl(:,:,:) * qd_pl(:,:,:)
    endif
    rhoqd_sum = GTL_global_sum( tmp, tmp_pl )

    !--- total mass (each water category)
    do nq = NQW_STR, NQW_END
       tmp(:,:,:) = rho(:,:,:) * q(:,:,:,nq)
       if ( ADM_prc_me == ADM_prc_pl ) then
          tmp_pl(:,:,:) = rho_pl(:,:,:) * q_pl(:,:,:,nq)
       endif
       rhoq_sum(nq) = GTL_global_sum( tmp, tmp_pl )
    enddo

    !--- total mass (total/vapor/liquid/soild water)
    rhoqw_sum = 0.0_RP
    rhoqv_sum = 0.0_RP
    rhoql_sum = 0.0_RP
    rhoqi_sum = 0.0_RP
    do nq = NQW_STR, NQW_END
       rhoqw_sum = rhoqw_sum + rhoq_sum(nq)

       if    ( nq == I_QV ) then
          rhoqv_sum = rhoqv_sum + rhoq_sum(nq)
       elseif( nq == I_QC .OR. nq == I_QR ) then
          rhoql_sum = rhoql_sum + rhoq_sum(nq)
       elseif( nq == I_QI .OR. nq == I_QS  .OR. nq == I_QG ) then
          rhoqi_sum = rhoqi_sum + rhoq_sum(nq)
       endif
    enddo

    !--- total mass (precipitation/evapolation)
    precip_sum = GTL_global_sum_srf( precip0(:,:,:), precip0_pl(:,:,:) )
    evap_sum   = GTL_global_sum_srf( evap0  (:,:,:), evap0_pl  (:,:,:) )

    if ( iflag ) then
       rhoqw_sum_old = rhoqw_sum
    endif

    if ( ADM_prc_me == ADM_prc_run_master ) then
       write(MNT_m_fid,'(I6)',    ADVANCE='NO') TIME_CSTEP
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') rhoqd_sum
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') rhoqv_sum
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') rhoql_sum
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') rhoqi_sum
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') rhoqw_sum
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') rhoqw_sum - rhoqw_sum_old
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') precip_sum
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') evap_sum
       write(MNT_m_fid,'(E22.14)',ADVANCE='NO') evap_sum - precip_sum
       write(MNT_m_fid,*)
    endif


    !----- Energy budget

    rhoein_sum = 0.0_RP

    !--- internal energy (dry air)
    tmp = rho * qd * CNST_CV * tem
    if ( ADM_prc_me == ADM_prc_pl ) then
       tmp_pl = rho_pl * qd_pl * CNST_CV * tem_pl
    end if
    rhoein_qd_sum = GTL_global_sum( tmp, tmp_pl )
    rhoein_sum    = rhoein_sum + rhoein_qd_sum

    !--- internal energy (each water category)
    do nq = NQW_STR,NQW_END
       tmp(:,:,:) = rho(:,:,:) * q(:,:,:,nq) * CVW(nq) * tem(:,:,:)

       !--- correct latent heat
       if    ( nq == I_QV ) then
          tmp(:,:,:) = tmp(:,:,:) + rho(:,:,:) * q(:,:,:,nq) * LHV
       elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          tmp(:,:,:) = tmp(:,:,:) - rho(:,:,:) * q(:,:,:,nq) * LHF
       endif

       if ( ADM_prc_me == ADM_prc_pl ) then
          tmp_pl(:,:,:) = rho_pl(:,:,:) * q_pl(:,:,:,nq) * CVW(nq) * tem_pl(:,:,:)

          if    ( nq == I_QV ) then
             tmp_pl(:,:,:) = tmp_pl(:,:,:) + rho_pl(:,:,:) * q_pl(:,:,:,nq) * LHV
          elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
             tmp_pl(:,:,:) = tmp_pl(:,:,:) - rho_pl(:,:,:) * q_pl(:,:,:,nq) * LHF
          endif
       endif

       rhoein_qw_sum(nq) = GTL_global_sum( tmp, tmp_pl )
       rhoein_sum        = rhoein_sum + rhoein_qw_sum(nq)
    enddo

    !--- potential energy
    tmp(:,:,:) = rho(:,:,:) * VMTR_PHI(:,:,:)
    if ( ADM_prc_me == ADM_prc_pl ) then
       tmp_pl(:,:,:) = rho_pl(:,:,:) * VMTR_PHI_pl(:,:,:)
    endif
    rhophi_sum = GTL_global_sum( tmp, tmp_pl )

    !--- kinetic energy
    call cnvvar_rhogkin( rhog,   rhog_pl,   & !--- [IN]
                         rhogvx, rhogvx_pl, & !--- [IN]
                         rhogvy, rhogvy_pl, & !--- [IN]
                         rhogvz, rhogvz_pl, & !--- [IN]
                         rhogw,  rhogw_pl,  & !--- [IN]
                         tmp,    tmp_pl     ) !--- [OUT]

    tmp(:,:,:) = tmp(:,:,:) *  VMTR_RGSGAM2(:,:,:)
    if ( ADM_prc_me == ADM_prc_pl ) then
       tmp_pl(:,:,:) = tmp_pl(:,:,:) * VMTR_RGSGAM2_pl(:,:,:)
    endif
    rhokin_sum = GTL_global_sum( tmp, tmp_pl )

    !--- total energy
    rhoetot_sum = rhoein_sum + rhophi_sum + rhokin_sum

    if ( iflag ) then
       rhoetot_sum_old = rhoetot_sum
    endif

    sfcrad_sum        = GTL_global_sum_srf( sfcrad0(:,:,:), sfcrad0_pl(:,:,:) )
    toarad_sum        = GTL_global_sum_srf( toarad0(:,:,:), toarad0_pl(:,:,:) )

    evap_energy_sum   = GTL_global_sum_srf( evap_energy0  (:,:,:), evap_energy0_pl  (:,:,:) )
    precip_energy_sum = GTL_global_sum_srf( precip_energy0(:,:,:), precip_energy0_pl(:,:,:) )

    sh_flux_sfc_sum   = GTL_global_sum_srf( sh_flux_sfc0(:,:,:), sh_flux_sfc0_pl(:,:,:) )
    lh_flux_sfc_sum   = GTL_global_sum_srf( lh_flux_sfc0(:,:,:), lh_flux_sfc0_pl(:,:,:) )

    if ( ADM_prc_me == ADM_prc_run_master ) then

       write(MNT_e_fid,'(I6)',    ADVANCE='NO') TIME_CSTEP
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') rhoein_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') rhophi_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') rhokin_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') rhoetot_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') rhoetot_sum-rhoetot_sum_old
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') sfcrad_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') toarad_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') sh_flux_sfc_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') lh_flux_sfc_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') evap_energy_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') precip_energy_sum
       write(MNT_e_fid,'(E22.14)',ADVANCE='NO') sfcrad_sum-toarad_sum &
                                              + evap_energy_sum-precip_energy_sum &
                                              + sh_flux_sfc_sum+lh_flux_sfc_sum
       write(MNT_e_fid,*)
    endif

    if( iflag ) iflag = .false.
    rhoqw_sum_old   = rhoqw_sum
    rhoetot_sum_old = rhoetot_sum

    ! reset array
    evap0     (:,:,:) = 0.0_RP
    evap0_pl  (:,:,:) = 0.0_RP
    precip0   (:,:,:) = 0.0_RP
    precip0_pl(:,:,:) = 0.0_RP

    sfcrad0          (:,:,:) = 0.0_RP
    sfcrad0_pl       (:,:,:) = 0.0_RP
    toarad0          (:,:,:) = 0.0_RP
    toarad0_pl       (:,:,:) = 0.0_RP
    evap_energy0     (:,:,:) = 0.0_RP
    evap_energy0_pl  (:,:,:) = 0.0_RP
    precip_energy0   (:,:,:) = 0.0_RP
    precip_energy0_pl(:,:,:) = 0.0_RP
    sh_flux_sfc0     (:,:,:) = 0.0_RP
    sh_flux_sfc0_pl  (:,:,:) = 0.0_RP
    lh_flux_sfc0     (:,:,:) = 0.0_RP
    lh_flux_sfc0_pl  (:,:,:) = 0.0_RP

    return
  end subroutine diagnose_energy_mass

end module mod_embudget
