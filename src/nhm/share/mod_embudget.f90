!-------------------------------------------------------------------------------
!> Module budget monitoring
!!
!! @par Description
!!          This module is for monitoring the energy/mass budget
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_embudget
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
  private :: store_var

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private :: MNT_ON   = .false.
  integer,  private :: MNT_INTV = 1

  character(len=H_LONG), private :: MNT_m_fname = 'BUDGET_mass.log'
  character(len=H_LONG), private :: MNT_e_fname = 'BUDGET_energy.log'

  integer,               private :: MNT_m_fid
  integer,               private :: MNT_e_fid

  real(RP),              private :: Mass_budget_factor
  real(RP),              private :: Energy_budget_factor



  real(RP), private :: atm_mass_qdry_sum_old = 0.0_RP
  real(RP), private :: atm_mass_qvap_sum_old = 0.0_RP
  real(RP), private :: atm_mass_qliq_sum_old = 0.0_RP
  real(RP), private :: atm_mass_qice_sum_old = 0.0_RP
  real(RP), private :: atm_mass_qtot_sum_old = 0.0_RP
  real(RP), private :: atm_engy_pot_sum_old  = 0.0_RP
  real(RP), private :: atm_engy_int_sum_old  = 0.0_RP
  real(RP), private :: atm_engy_kin_sum_old  = 0.0_RP
  real(RP), private :: atm_engy_sum_old      = 0.0_RP

  logical,  private :: first = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine embudget_setup
    use mod_process, only: &
       PRC_IsMaster, &
       PRC_MPIstop
    use mod_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI
    use mod_time, only: &
       TIME_DTL
    implicit none

    namelist / EMBUDGETPARAM / &
       MNT_INTV, &
       MNT_ON

    integer  :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[embudget]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=EMBUDGETPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** EMBUDGETPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist EMBUDGETPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=EMBUDGETPARAM)

    if(.not.MNT_ON) return

    Mass_budget_factor   = 1.0_RP / ( 4.0_RP * PI * RADIUS * RADIUS )                                     ! [J]       -> [J/m2]
    Energy_budget_factor = 1.0_RP / ( TIME_DTL * real(MNT_INTV,kind=RP) * 4.0_RP * PI * RADIUS * RADIUS ) ! [J /step] -> [W/m2]

    if( IO_L ) write(IO_FID_LOG,*) "Mass_budget_factor   = ", Mass_budget_factor
    if( IO_L ) write(IO_FID_LOG,*) "Energy_budget_factor = ", Energy_budget_factor

    ! open budget.info file
    if ( PRC_IsMaster ) then
       MNT_m_fid = IO_get_available_fid()

       open( unit   = MNT_m_fid,         &
             file   = trim(MNT_m_fname), &
             status = 'replace',         &
             form   = 'formatted'        )

       MNT_e_fid = IO_get_available_fid()

       open( unit   = MNT_e_fid,         &
             file   = trim(MNT_e_fname), &
             status = 'replace',         &
             form   = 'formatted'        )
    endif

    call store_var

    return
  end subroutine embudget_setup

  !-----------------------------------------------------------------------------
  subroutine embudget_monitor
    use mod_time, only: &
       TIME_CSTEP
    implicit none
    !---------------------------------------------------------------------------

    if( .NOT. MNT_ON ) return

    if ( mod(TIME_CSTEP,MNT_INTV) == 0 ) then
       call store_var
    endif

    return
  end subroutine embudget_monitor

  !-----------------------------------------------------------------------------
  subroutine store_var
    use mod_process, only: &
       PRC_IsMaster
    use mod_const, only: &
       CVdry => CONST_CVdry,  &
       LHV   => CONST_LHV,    &
       LHF   => CONST_LHF
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_vmtr, only: &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_PHI,        &
       VMTR_PHI_pl
    use mod_time, only: &
       TIME_CSTEP, &
       TIME_DTL
    use mod_statistics, only: &
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
       CVW
    use mod_prgvar, only: &
       prgvar_get_withdiag
    use mod_cnvvar, only: &
       cnvvar_rhogkin
    use mod_thrmdyn, only: &
       THRMDYN_qd
    implicit none

    real(RP) :: rhog        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhoge       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhoge_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogq       (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: rhogq_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    real(RP) :: rho         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: pre         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: pre_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tem         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tem_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vx          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vx_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vy          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vy_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vz          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vz_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: w           (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: w_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q           (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: q_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    real(RP) :: qd          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qd_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tmp         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tmp_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhoq_sum    (TRC_vmax)

    real(RP) :: atm_engy_pot_sum
    real(RP) :: atm_engy_int_sum
    real(RP) :: atm_engy_kin_sum
    real(RP) :: atm_engy_sum

    real(RP) :: atm_engy_pot_sum_diff
    real(RP) :: atm_engy_int_sum_diff
    real(RP) :: atm_engy_kin_sum_diff
    real(RP) :: atm_engy_sum_diff

    real(RP) :: atm_mass_qdry_sum
    real(RP) :: atm_mass_qvap_sum
    real(RP) :: atm_mass_qliq_sum
    real(RP) :: atm_mass_qice_sum
    real(RP) :: atm_mass_qtot_sum

    real(RP) :: atm_mass_qdry_sum_diff
    real(RP) :: atm_mass_qvap_sum_diff
    real(RP) :: atm_mass_qliq_sum_diff
    real(RP) :: atm_mass_qice_sum_diff
    real(RP) :: atm_mass_qtot_sum_diff

    integer  :: nq
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

    if ( ADM_have_pl ) then
       call THRMDYN_qd( ADM_gall_pl,    & ! [IN]
                        ADM_kall,       & ! [IN]
                        ADM_lall_pl,    & ! [IN]
                        q_pl (:,:,:,:), & ! [IN]
                        qd_pl(:,:,:)    ) ! [OUT]
    endif

    !##### energy budget in the atmosphere #####

    !--- potential energy
    tmp(:,:,:) = rho(:,:,:) * VMTR_PHI(:,:,:)
    if ( ADM_have_pl ) then
       tmp_pl(:,:,:) = rho_pl(:,:,:) * VMTR_PHI_pl(:,:,:)
    endif

    atm_engy_pot_sum = GTL_global_sum( tmp(:,:,:), tmp_pl(:,:,:) )

    !--- internal energy
    tmp(:,:,:) = qd(:,:,:) * CVdry * tem(:,:,:)
    do nq = NQW_STR, NQW_END
       tmp(:,:,:) = tmp(:,:,:) + q(:,:,:,nq) * CVW(nq) * tem(:,:,:)

       if ( nq == I_QV ) then
          tmp(:,:,:) = tmp(:,:,:) + q(:,:,:,nq) * LHV
       endif

       if ( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          tmp(:,:,:) = tmp(:,:,:) - q(:,:,:,nq) * LHF
       endif
    enddo
    tmp(:,:,:) = rho(:,:,:) * tmp(:,:,:)

    if ( ADM_have_pl ) then
       tmp_pl(:,:,:) = qd_pl(:,:,:) * CVdry * tem_pl(:,:,:)
       do nq = NQW_STR, NQW_END
          tmp_pl(:,:,:) = tmp_pl(:,:,:) + q_pl(:,:,:,nq) * CVW(nq) * tem_pl(:,:,:)

          if ( nq == I_QV ) then
             tmp_pl(:,:,:) = tmp_pl(:,:,:) + q_pl(:,:,:,nq) * LHV
          endif

          if ( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
             tmp_pl(:,:,:) = tmp_pl(:,:,:) - q_pl(:,:,:,nq) * LHF
          endif
       enddo
       tmp_pl(:,:,:) = rho_pl(:,:,:) * tmp_pl(:,:,:)
    endif

    atm_engy_int_sum = GTL_global_sum( tmp(:,:,:), tmp_pl(:,:,:) )

    !--- kinetic energy
    call cnvvar_rhogkin( rhog,   rhog_pl,   & ! [IN]
                         rhogvx, rhogvx_pl, & ! [IN]
                         rhogvy, rhogvy_pl, & ! [IN]
                         rhogvz, rhogvz_pl, & ! [IN]
                         rhogw,  rhogw_pl,  & ! [IN]
                         tmp,    tmp_pl     ) ! [OUT]

    tmp(:,:,:) = tmp(:,:,:) * VMTR_RGSGAM2(:,:,:)

    if ( ADM_have_pl ) then
       tmp_pl(:,:,:) = tmp_pl(:,:,:) * VMTR_RGSGAM2_pl(:,:,:)
    endif

    atm_engy_kin_sum = GTL_global_sum( tmp(:,:,:), tmp_pl(:,:,:) )

    !--- total energy
    atm_engy_sum = atm_engy_pot_sum+ atm_engy_int_sum + atm_engy_kin_sum

    !##### mass budget in the atmosphere #####

    !--- dry air density
    tmp(:,:,:) = rho(:,:,:) * qd(:,:,:)
    if ( ADM_have_pl ) then
       tmp_pl(:,:,:) = rho_pl(:,:,:) * qd_pl(:,:,:)
    endif
    atm_mass_qdry_sum = GTL_global_sum( tmp(:,:,:), tmp_pl(:,:,:) )

    !--- water tracer density
    do nq = NQW_STR, NQW_END
       tmp(:,:,:) = rho(:,:,:) * q(:,:,:,nq)
       if ( ADM_have_pl ) then
          tmp_pl(:,:,:) = rho_pl(:,:,:) * q_pl(:,:,:,nq)
       endif
       rhoq_sum(nq) = GTL_global_sum( tmp, tmp_pl )
    enddo

    !--- total & subtotal water
    atm_mass_qtot_sum = 0.0_RP
    atm_mass_qvap_sum = 0.0_RP
    atm_mass_qliq_sum = 0.0_RP
    atm_mass_qice_sum = 0.0_RP
    do nq = NQW_STR, NQW_END
       atm_mass_qtot_sum = atm_mass_qtot_sum + rhoq_sum(nq)

       if ( nq == I_QV ) then
          atm_mass_qvap_sum = atm_mass_qvap_sum + rhoq_sum(nq)
       endif

       if ( nq == I_QC .OR. nq == I_QR ) then
          atm_mass_qliq_sum = atm_mass_qliq_sum + rhoq_sum(nq)
       endif

       if ( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          atm_mass_qice_sum = atm_mass_qice_sum + rhoq_sum(nq)
       endif
    enddo



    !##### File OUTPUT #####

    if ( first ) then
       ! [J/m2], absolute value
       atm_engy_pot_sum_diff  = atm_engy_pot_sum  * Mass_budget_factor
       atm_engy_int_sum_diff  = atm_engy_int_sum  * Mass_budget_factor
       atm_engy_kin_sum_diff  = atm_engy_kin_sum  * Mass_budget_factor
       atm_engy_sum_diff      = atm_engy_sum      * Mass_budget_factor
       ! [kg/m2], absolute value
       atm_mass_qdry_sum_diff = atm_mass_qdry_sum * Mass_budget_factor
       atm_mass_qvap_sum_diff = atm_mass_qvap_sum * Mass_budget_factor
       atm_mass_qliq_sum_diff = atm_mass_qliq_sum * Mass_budget_factor
       atm_mass_qice_sum_diff = atm_mass_qice_sum * Mass_budget_factor
       atm_mass_qtot_sum_diff = atm_mass_qtot_sum * Mass_budget_factor
    else
       ! [W/m2], difference from previous step
       atm_engy_pot_sum_diff  = ( atm_engy_pot_sum  - atm_engy_pot_sum_old  ) * Energy_budget_factor
       atm_engy_int_sum_diff  = ( atm_engy_int_sum  - atm_engy_int_sum_old  ) * Energy_budget_factor
       atm_engy_kin_sum_diff  = ( atm_engy_kin_sum  - atm_engy_kin_sum_old  ) * Energy_budget_factor
       atm_engy_sum_diff      = ( atm_engy_sum      - atm_engy_sum_old      ) * Energy_budget_factor
       ! [kg/m2/s], difference from previous step
       atm_mass_qdry_sum_diff = ( atm_mass_qdry_sum - atm_mass_qdry_sum_old ) * Energy_budget_factor
       atm_mass_qvap_sum_diff = ( atm_mass_qvap_sum - atm_mass_qvap_sum_old ) * Energy_budget_factor
       atm_mass_qliq_sum_diff = ( atm_mass_qliq_sum - atm_mass_qliq_sum_old ) * Energy_budget_factor
       atm_mass_qice_sum_diff = ( atm_mass_qice_sum - atm_mass_qice_sum_old ) * Energy_budget_factor
       atm_mass_qtot_sum_diff = ( atm_mass_qtot_sum - atm_mass_qtot_sum_old ) * Energy_budget_factor
    endif

    if ( PRC_IsMaster ) then

       !--- BUDGET_energy : energy budget in the atmosphere
       if ( first ) then
          write(MNT_e_fid,'(A6)' ,advance='no') '#STEP'
          write(MNT_e_fid,'(A16)',advance='no') 'Day'
          write(MNT_e_fid,'(A16)',advance='no') 'Eng(pot)'
          write(MNT_e_fid,'(A16)',advance='no') 'Eng(int)'
          write(MNT_e_fid,'(A16)',advance='no') 'Eng(kin)'
          write(MNT_e_fid,'(A16)',advance='no') 'Eng(tot)'
          write(MNT_e_fid,*)
       endif

       write(MNT_e_fid,'(I6)'    ,advance='no') TIME_CSTEP
       write(MNT_e_fid,'(ES16.8)',advance='no') TIME_CSTEP * TIME_DTL / 86400.0_RP
       write(MNT_e_fid,'(ES16.8)',advance='no') atm_engy_pot_sum_diff
       write(MNT_e_fid,'(ES16.8)',advance='no') atm_engy_int_sum_diff
       write(MNT_e_fid,'(ES16.8)',advance='no') atm_engy_kin_sum_diff
       write(MNT_e_fid,'(ES16.8)',advance='no') atm_engy_sum_diff
       write(MNT_e_fid,*)

       !--- BUDGET_mass : mass budget in the atmosphere
       if ( first ) then
          write(MNT_m_fid,'(A6)' ,advance='no') '#STEP'
          write(MNT_m_fid,'(A16)',advance='no') 'Day'
          write(MNT_m_fid,'(A16)',advance='no') 'dry air mass '
          write(MNT_m_fid,'(A16)',advance='no') 'water mass(g)'
          write(MNT_m_fid,'(A16)',advance='no') 'water mass(l)'
          write(MNT_m_fid,'(A16)',advance='no') 'water mass(s)'
          write(MNT_m_fid,'(A16)',advance='no') 'water mass(t)'
          write(MNT_m_fid,*)
       endif

       write(MNT_m_fid,'(I6)'    ,advance='no') TIME_CSTEP
       write(MNT_m_fid,'(ES16.8)',advance='no') TIME_CSTEP * TIME_DTL / 86400.0_RP
       write(MNT_m_fid,'(ES16.8)',advance='no') atm_mass_qdry_sum_diff
       write(MNT_m_fid,'(ES16.8)',advance='no') atm_mass_qvap_sum_diff
       write(MNT_m_fid,'(ES16.8)',advance='no') atm_mass_qliq_sum_diff
       write(MNT_m_fid,'(ES16.8)',advance='no') atm_mass_qice_sum_diff
       write(MNT_m_fid,'(ES16.8)',advance='no') atm_mass_qtot_sum_diff
       write(MNT_m_fid,*)
    endif

    atm_engy_pot_sum_old  = atm_engy_pot_sum
    atm_engy_int_sum_old  = atm_engy_int_sum
    atm_engy_kin_sum_old  = atm_engy_kin_sum
    atm_engy_sum_old      = atm_engy_sum

    atm_mass_qdry_sum_old = atm_mass_qdry_sum
    atm_mass_qvap_sum_old = atm_mass_qvap_sum
    atm_mass_qliq_sum_old = atm_mass_qliq_sum
    atm_mass_qice_sum_old = atm_mass_qice_sum
    atm_mass_qtot_sum_old = atm_mass_qtot_sum

    if ( first ) then
       first = .false.
    endif

    return
  end subroutine store_var

end module mod_embudget
