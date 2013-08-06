!-------------------------------------------------------------------------------
!
!+  Program driver
!
!-------------------------------------------------------------------------------
program prg_driver
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This program is a driver of non-hydrostatic model based on an 
  !       icosahedral grid system.
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !      0.03      04-05-31   Change by addtion of mod[onestep].
  !                05-12-01   M.Satoh add history_setup
  !                05-12-19   S.Iga moved output_timeinfo after output_all
  !                06-04-18   T.Mitsui add sfc_restart
  !                06-04-21   H.Tomita:  remove output_timeinfo due to 
  !                                      computational efficeincy.
  !                                      Instead, this process is move to 
  !                                      mod[mod_output].
  !                06-08-07   W.Yanase add history_vars
  !                06-09-27   S.Iga add history_vars_cfmip
  !                07-03-23   Y.Niwa add ndg_setup, ndg_do, FLAG_NUDGING
  !                07-06-27   Y.Niwa add history_vars_setup
  !                07-07-24   K.Suzuki: implementing SPRINTARS aerosol model
  !                07-08-06   Y.Niwa: add history_obs, history_vars_obs
  !                07-11-07   T.Mitsui: add option to omit output_all
  !                08-03-10   T.Mitsui: add intermediate output of restart file
  !                08-05-24   T.Mitsui: trivial fix
  !                08-09-09   Y.Niwa : modfied for nudging
  !                09-01-23   H.Tomita: a) abolish mod_output, mod_extdata
  !                                     mod_history_vars_cfmip, mod_o3var.
  !                                     b) introduce mod_extdata2.
  !                09-04-14   T.Mitsui: arrange initialization of aerosols
  !                09-07-10   H.Tomita: Add the module [mod_embudget].
  !                09-08-05   S.Iga: remove latlon_setup (suggested by T.Mitsui)
  !                09-08-05   T.Mitsui: add conditioning by ADM_myprc_is_run
  !                                     to keep out extra-processes from main routines.
  !                09-08-18   T.Mitsui: change of 09-08-05 is not enough.
  !                10-03-08   C.Kodama: Modify for overwrite_restart option
  !                10-04-30   M.Satoh: move diagvar_setup
  !                11-09-03   H.Yashiro : New I/O
  !                11-11-28   Y.Yamada : merge Terai-san timer
  !                12-06-07   T.Seiki  : Application to Multi-job System
  !                12-10-12   R.Yoshida  : Modify for Dynamical Core test
  !                12-10-22   R.Yoshida  : add papi instructions
  !      ----------------------------------------------------------------------- 
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_debug
  use mod_adm, only: &
     ADM_MULTI_PRC,      &
     ADM_LOG_FID,        &
     ADM_prc_me,         &
     ADM_prc_run_master, &
     ADM_proc_init,      &
     ADM_proc_stop,      &
     ADM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_cnst, only: &
     CNST_setup
  use mod_calendar, only: &
     calendar_setup
  use mod_time, only: &
     TIME_setup,          &
     TIME_report,         &
     TIME_advance,        &
     TIME_LSTEP_MAX,      &
     cstep => TIME_CSTEP, &
     ctime => TIME_CTIME, &
     dtime => TIME_DTL
  use mod_grd, only: &
     GRD_setup
  use mod_gmtr, only: &
     GMTR_setup
  use mod_oprt, only: &
     OPRT_setup
  use mod_vmtr, only: &
     VMTR_setup

  use mod_runconf, only: &
     runconf_setup, &
     FLAG_NUDGING
  use mod_prgvar, only: &
     prgvar_setup,            &
     restart_input_basename,  &
     restart_output_basename, &
     restart_input,           &
     restart_output
  use mod_diagvar, only: &
       diagvar_setup, &
       diagvar_restart_output
  use mod_sfcvar, only: &
       sfcvar_setup
  use mod_bsstate, only: &
       bsstate_setup
  use mod_bndcnd, only   :   &
       bndcnd_setup
  use mod_numfilter, only  : &
       numfilter_setup
  use mod_forcing_driver, only : &
       forcing_init, &
       forcing
  use mod_ndg, only: &
       ndg_setup
  use mod_dynstep, only : &
       dynstep
  use mod_history, only: &
       history_setup, &
       history_out,   &
       HIST_output_step0
  use mod_history_vars, only: &
       history_vars_setup, &
       history_vars
  use mod_embudget, only: &
       embudget_setup, &
       embudget_monitor
  implicit none

  character(len=14) :: cdate

  integer :: n
  !-----------------------------------------------------------------------------

  call ADM_proc_init(ADM_MULTI_PRC)

  !---< admin module setup >---
  call ADM_setup('nhm_driver.cnf')

  !#############################################################################

  write(ADM_LOG_FID,*) '##### start  setup     #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### start  setup     #####'
  endif

  call DEBUG_rapstart('Total')
  call DEBUG_rapstart('Setup ALL')

  !---< cnst module setup >---
  call CNST_setup

  !---< calendar module setup >---
  call calendar_setup

  !---< I/O module setup >---
  call FIO_setup

  !---< comm module setup >---
  call COMM_setup

  !---< grid module setup >---
  call GRD_setup

  !---< geometrics module setup >---
  call GMTR_setup

  !---< operator module setup >---
  call OPRT_setup

  !---< vertical metrics module setup >---
  call VMTR_setup

  !---< time module setup >---
  call TIME_setup


  !---< nhm_runconf module setup >---
  call runconf_setup

  !---< prognostic variable module setup >---
  call prgvar_setup
  call restart_input( restart_input_basename )

  !---< diagnostic variable module setup >---
  call diagvar_setup

  !---< surface variable module setup >---
  call sfcvar_setup


  !---< boundary condition module setup >---
  call bndcnd_setup

  !---< basic state module setup >---
  call bsstate_setup

  !---< numerical filter module setup >---
  call numfilter_setup

  !---< forcing module setup >---
  call forcing_init

  !---< energy&mass budget module setup >---
  call embudget_setup

  !---< history module setup >---
  call history_setup

  !---< history variable module setup >---
  call history_vars_setup

  !---< nudging module setup >---
  if( FLAG_NUDGING ) call ndg_setup( ctime, dtime )

  write(ADM_LOG_FID,*) '##### finish setup     #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### finish setup     #####'
  endif

  call DEBUG_rapend('Setup ALL')

  !#############################################################################
#ifdef _FIPP_
  call fipp_start()
#endif
  call DEBUG_rapstart('Main ALL')

  write(ADM_LOG_FID,*) '##### start  main loop #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### start  main loop #####'
  endif

  call TIME_report

  !--- history output at initial time
  if ( HIST_output_step0 ) then
     cstep = -1
     ctime = -dtime
     call history_vars
     call TIME_advance
     call history_out
  endif

  do n = 1, TIME_LSTEP_MAX

     call DEBUG_rapstart('+Atmos')
     call dynstep
     call forcing
     call DEBUG_rapend  ('+Atmos')

     call DEBUG_rapstart('+History')
     call history_vars
     call TIME_advance

     !--- budget monitor
     call embudget_monitor
     call history_out

     if (n == TIME_LSTEP_MAX) then
        cdate = ""
        call restart_output( restart_output_basename )
        call diagvar_restart_output ( cdate )
     endif
     call DEBUG_rapend  ('+History')

  enddo

  write(ADM_LOG_FID,*) '##### finish main loop #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### finish main loop #####'
  endif

  call DEBUG_rapend('Main ALL')
#ifdef _FIPP_
  call fipp_stop()
#endif
  !#############################################################################

  call DEBUG_rapend('Total')
  call DEBUG_rapreport

  !--- all processes stop
  call ADM_proc_stop

  stop
end program prg_driver
!-------------------------------------------------------------------------------
