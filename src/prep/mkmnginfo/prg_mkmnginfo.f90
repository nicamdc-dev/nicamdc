!-------------------------------------------------------------------------------
!>
!! Program mkmnginfo
!!
!! @par Description
!!         This program makes managing infomation file for paralell computation.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2007-10-22 (T.Mitsui)  change value of rgn_nlim
!! @li      2010-06-07 (S.Iga)     new grid (Iga 2010) is implemented. (see string XTMS)
!! @li      2011-07-21 (T.Ohno)    2 new grid systems (1DMD-ON-SPHERE are added by Hara-san@JAMSTEC)
!! @li      2012-06-11 (H.Yashiro) Milestone-project, code cleanup
!!
!<
Program prg_mkmnginfo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_debug
  use mod_adm, only: &
     ADM_SINGLE_PRC, &
     ADM_proc_init,  &
     ADM_proc_stop,  &
     ADM_setup,      &
     ADM_LOG_FID,    &
     ADM_CTL_FID,    &
     ADM_MAXFNAME,   &
     ADM_rgn_etab,   &
     ADM_prc_rnum,   &
     ADM_prc_tab
  use mod_rgnmng, only: &
     RGNMNG_output
  !-----------------------------------------------------------------------------
  implicit none

  character(len=ADM_MAXFNAME) :: rgnmng_out_fname = '' !< output name of region-management file

  namelist / MKMNGINFOPARAM / &
       rgnmng_out_fname

  !=============================================================================

  call ADM_proc_init(ADM_MULTI_PRC)

  call DEBUG_rapstart('Total')

  !---< admin module setup >---
  call ADM_setup('mkmnginfo.cnf')

  !--- read parameters
  rewind(ADM_CTL_FID)
  read(ADM_CTL_FID,nml=MKMNGINFOPARAM,iostat=ierr)
  if ( ierr < 0 ) then
     write(ADM_LOG_FID,*) '*** MKMNGINFOPARAM is not specified. use default.'
  elseif( ierr > 0 ) then
     write(*,          *) 'xxx Not appropriate names in namelist MKMNGINFOPARAM. STOP.'
     write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist MKMNGINFOPARAM. STOP.'
     call ADM_proc_stop
  endif
  write(ADM_LOG_FID,MKMNGINFOPARAM)

  if ( rgnmng_out_fname == '' ) then
     write(*,          *) 'xxx rgnmng_out_fname is not specified. STOP.'
     write(ADM_LOG_FID,*) 'xxx rgnmng_out_fname is not specified. STOP.'
     call ADM_proc_stop
  endif

  call RGNMNG_output( rgnmng_out_fname, &
                      ADM_rgn_etab,     &
                      ADM_prc_rnum,     &
                      ADM_prc_tab       )

  call DEBUG_rapend('Total')
  call DEBUG_rapreport

  !--- all processes stop
  call ADM_proc_stop

  stop
end program prg_mkmnginfo
!-------------------------------------------------------------------------------
