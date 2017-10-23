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
  !                07-03-23   Y.Niwa add NDG_setup, ndg_do, FLAG_NUDGING
  !                07-06-27   Y.Niwa add history_vars_setup
  !                07-07-24   K.Suzuki: implementing SPRINTARS aerosol model
  !                07-08-06   Y.Niwa: add history_obs, history_vars_obs
  !                07-11-07   T.Mitsui: add option to omit output_all
  !                08-03-10   T.Mitsui: add intermediate output of restart file
  !                08-05-24   T.Mitsui: trivial fix
  !                08-09-09   Y.Niwa : modfied for nudging
  !                09-01-23   H.Tomita: a) abolish mod_output, mod_extdata
  !                                     mod_history_vars_cfmip, mod_o3var.
  !                                     b) introduce mod_extdata.
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
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_process, only: &
     PRC_IsMaster,    &
     PRC_MPIstart,    &
     PRC_LOCAL_setup, &
     PRC_MPIfinish
  use mod_const, only: &
     CONST_setup
  use mod_calendar, only: &
     CALENDAR_setup
  use mod_random, only: &
     RANDOM_setup
  use mod_adm, only: &
     ADM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_hio, only: &
     HIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_grd, only: &
     GRD_setup
  use mod_gmtr, only: &
     GMTR_setup
  use mod_oprt, only: &
     OPRT_setup
  use mod_vmtr, only: &
     VMTR_setup
  use mod_time, only: &
     TIME_setup,     &
     TIME_report,    &
     TIME_advance,   &
     TIME_LSTEP_MAX, &
     TIME_CSTEP,     &
     TIME_CTIME,     &
     TIME_DTL
  use mod_extdata, only: &
     extdata_setup
  use mod_runconf, only: &
     runconf_setup
  use mod_saturation, only: &
     saturation_setup
  use mod_prgvar, only: &
     prgvar_setup,            &
     restart_input_basename,  &
     restart_output_basename, &
     restart_input,           &
     restart_output
  use mod_dynamics, only: &
     dynamics_setup, &
     dynamics_step
  use mod_forcing_driver, only: &
     forcing_setup, &
     forcing_step
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

  !##### OpenACC (for data copy) #####
  use mod_adm, only: &
     RGNMNG_l2r,  &
     RGNMNG_vert_num
  use mod_comm, only: &
     Send_info_r2r, Send_info_p2r, Send_info_r2p, &
     Send_list_r2r, Send_list_p2r, Send_list_r2p, &
     Recv_info_r2r, Recv_info_p2r, Recv_info_r2p, &
     Recv_list_r2r, Recv_list_p2r, Recv_list_r2p, &
     Copy_info_r2r, Copy_info_p2r, Copy_info_r2p, &
     Copy_list_r2r, Copy_list_p2r, Copy_list_r2p, &
     sendbuf_r2r_SP, sendbuf_p2r_SP, sendbuf_r2p_SP, &
     recvbuf_r2r_SP, recvbuf_p2r_SP, recvbuf_r2p_SP, &
     sendbuf_r2r_DP, sendbuf_p2r_DP, sendbuf_r2p_DP, &
     recvbuf_r2r_DP, recvbuf_p2r_DP, recvbuf_r2p_DP, &
     Singular_info, Singular_list, REQ_list
  use mod_grd, only: &
     GRD_x,     &
     GRD_xt,    &
     GRD_zs,    &
     GRD_rdgz,  &
     GRD_rdgzh, &
     GRD_vz
  use mod_gmtr, only: &
     GMTR_p, &
     GMTR_t, &
     GMTR_a
  use mod_vmtr, only: &
     VMTR_GAM2H,     &
     VMTR_GSGAM2,    &
     VMTR_GSGAM2H,   &
     VMTR_RGSQRTH,   &
     VMTR_RGAM,      &
     VMTR_RGAMH,     &
     VMTR_RGSGAM2,   &
     VMTR_RGSGAM2H,  &
     VMTR_W2Cfact,   &
     VMTR_C2Wfact,   &
     VMTR_C2WfactGz, &
     VMTR_PHI
  use mod_runconf, only: &
     CVW
  use mod_prgvar, only: &
     PRG_var,  &
     DIAG_var
  use mod_bsstate, only: &
     rho_bs, &
     pre_bs, &
     tem_bs
  use mod_numfilter, only: &
     Kh_coef,      &
     Kh_coef_lap1, &
     divdamp_coef
  use mod_vi, only: &
     Mc, &
     Ml, &
     Mu
  use mod_history, only: &
     ksumstr,     &
     cnvpre_klev, &
     cnvpre_fac1, &
     cnvpre_fac2
  !##### OpenACC #####
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  integer  :: comm_world
  integer  :: myrank
  logical  :: ismaster

  integer  :: n
  !=============================================================================

  !---< MPI start >---
  call PRC_MPIstart( comm_world ) ! [OUT]

  !---< STDIO setup >---
  call IO_setup( 'NICAM-DC',      & ! [IN]
                 'nhm_driver.cnf' ) ! [IN]

  !---< Local process management setup >---
  call PRC_LOCAL_setup( comm_world, & ! [IN]
                        myrank,     & ! [OUT]
                        ismaster    ) ! [OUT]

  !---< Logfile setup >---
  call IO_LOG_setup( myrank,  & ! [IN]
                     ismaster ) ! [IN]

  !---< profiler module setup >---
  call PROF_setup

  !#############################################################################
  call PROF_setprefx('INIT')
  call PROF_rapstart('Initialize',0)

  if( IO_L ) write(IO_FID_LOG,*) '##### start  setup     #####'
  if( PRC_IsMaster ) write(*,*)  '##### start  setup     #####'

  !---< cnst module setup >---
  call CONST_setup

  !---< calendar module setup >---
  call CALENDAR_setup

  !---< radom module setup >---
  call RANDOM_setup

  !---< admin module setup >---
  call ADM_setup

  !---< I/O module setup >---
  call FIO_setup
  call HIO_setup

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

  !---< external data module setup >---
  call extdata_setup


  !---< nhm_runconf module setup >---
  call runconf_setup

  !---< saturation module setup >---
  call saturation_setup

  !---< prognostic variable module setup >---
  call prgvar_setup
  call restart_input( restart_input_basename )


  !---< dynamics module setup >---
  call dynamics_setup

  !---< forcing module setup >---
  call forcing_setup

  !---< energy&mass budget module setup >---
  call embudget_setup

  !---< history module setup >---
  call history_setup

  !---< history variable module setup >---
  call history_vars_setup

  if( IO_L ) write(IO_FID_LOG,*) '##### finish setup     #####'
  if( PRC_IsMaster ) write(*,*)  '##### finish setup     #####'

  call PROF_rapend('Initialize',0)
  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Main_Loop',0)

  if( IO_L ) write(IO_FID_LOG,*) '##### start  main loop #####'
  if( PRC_IsMaster ) write(*,*)  '##### start  main loop #####'

#ifdef _FIPP_
  call fipp_start()
#endif

  !$acc data &
  !$acc& pcopyin(RGNMNG_l2r,RGNMNG_vert_num) &
  !$acc& pcopyin(Send_info_r2r,Send_info_p2r,Send_info_r2p) &
  !$acc& pcopyin(Send_list_r2r,Send_list_p2r,Send_list_r2p) &
  !$acc& pcopyin(Recv_info_r2r,Recv_info_p2r,Recv_info_r2p) &
  !$acc& pcopyin(Recv_list_r2r,Recv_list_p2r,Recv_list_r2p) &
  !$acc& pcopyin(Copy_info_r2r,Copy_info_p2r,Copy_info_r2p) &
  !$acc& pcopyin(Copy_list_r2r,Copy_list_p2r,Copy_list_r2p) &
  !$acc& pcopyin(sendbuf_r2r,sendbuf_p2r,sendbuf_r2p) &
  !$acc& pcopyin(recvbuf_r2r,recvbuf_p2r,recvbuf_r2p) &
  !$acc& pcopyin(Singular_info,Singular_list,REQ_list) &
  !$acc& pcopyin(GRD_rdgz,GRD_rdgzh,GRD_x,GRD_xt,GRD_vz,GRD_zs) &
  !$acc& pcopyin(GMTR_p,GMTR_t,GMTR_a) &
  !$acc& pcopyin(cdiv,cgrad,clap,cinterp_TN,cinterp_HN,cinterp_TRA,cinterp_PRA) &
  !$acc& pcopyin(VMTR_GAM2,VMTR_GAM2H,VMTR_GSGAM2,VMTR_GSGAM2H) &
  !$acc& pcopyin(VMTR_RGSQRTH,VMTR_RGAM,VMTR_RGAMH,VMTR_RGSGAM2,VMTR_RGSGAM2H) &
  !$acc& pcopyin(VMTR_W2Cfact,VMTR_C2Wfact,VMTR_C2WfactGz,VMTR_PHI) &
  !$acc& pcopyin(CVW) &
  !$acc& pcopyin(rho_bs,pre_bs,tem_bs) &
  !$acc& pcopyin(divdamp_coef,Kh_coef,Kh_coef_lap1) &
  !$acc& pcopyin(Mc,Mu,Ml) &
  !$acc& pcopyin(ksumstr,cnvpre_klev,cnvpre_fac1,cnvpre_fac2) &
  !$acc& pcopy  (PRG_var,DIAG_var)

  !--- history output at initial time
  if ( HIST_output_step0 ) then
     TIME_CSTEP = TIME_CSTEP - 1
     TIME_CTIME = TIME_CTIME - TIME_DTL
     call history_vars
     call TIME_advance
     call history_out
  else
     call TIME_report
  endif

  do n = 1, TIME_LSTEP_MAX

     call PROF_rapstart('_Atmos',1)
     call dynamics_step
     call forcing_step
     call PROF_rapend  ('_Atmos',1)

     call PROF_rapstart('_History',1)
     call history_vars
     call TIME_advance

     !--- budget monitor
     call embudget_monitor
     call history_out

     if ( n == TIME_LSTEP_MAX ) then
        call restart_output( restart_output_basename )
     endif
     call PROF_rapend  ('_History',1)

  enddo

  !$acc end data

#ifdef _FIPP_
  call fipp_stop()
#endif

  if( IO_L ) write(IO_FID_LOG,*) '##### finish main loop #####'
  if( PRC_IsMaster ) write(*,*)  '##### finish main loop #####'

  call PROF_rapend('Main_Loop',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

  stop
end program prg_driver
!-------------------------------------------------------------------------------
