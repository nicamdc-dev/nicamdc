!-------------------------------------------------------------------------------
!> Program mkllmap
!!
!! @par Description
!!          Making remapping coefficient between lat-lon and icosahedral grid
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
program prg_mkllmap
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_process, only: &
     PRC_MPIstart,    &
     PRC_LOCAL_setup, &
     PRC_MPIstop,     &
     PRC_MPIfinish
  use mod_const, only: &
     CONST_setup
  use mod_adm, only: &
     ADM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_grd, only: &
     GRD_setup
  use mod_latlon, only: &
     LATLON_setup, &
     LATLON_ico_setup
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  integer  :: comm_world
  integer  :: myrank
  logical  :: ismaster

  character(len=H_LONG) :: output_dir   = './'

  namelist /mkllmap_param/ &
     output_dir

  integer  :: ierr
  !=============================================================================

  !---< MPI start >---
  call PRC_MPIstart( comm_world ) ! [OUT]

  !---< STDIO setup >---
  call IO_setup( 'NICAM-DC',   & ! [IN]
                 'mkllmap.cnf' ) ! [IN]

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

  !--- < cnst module setup > ---
  call CONST_setup

  !--- < admin module setup > ---
  call ADM_setup

  !---< I/O module setup >---
  call FIO_setup

  !--- < comm module setup > ---
  call COMM_setup

  !--- < grid module setup > ---
  call GRD_setup

  !--- read parameters
  if( IO_L ) write(IO_FID_LOG,*)
  if( IO_L ) write(IO_FID_LOG,*) '+++ Program[mkllmap]/Category[tool]'
  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=MKLLMAP_PARAM,iostat=ierr)
  if ( ierr < 0 ) then
     if( IO_L ) write(IO_FID_LOG,*) '*** MKLLMAP_PARAM is not specified. use default.'
  elseif( ierr > 0 ) then
     if( IO_L ) write(*,*) 'xxx Not appropriate names in namelist MKLLMAP_PARAM. STOP.'
     call PRC_MPIstop
  endif
  if( IO_NML ) write(IO_FID_LOG,nml=MKLLMAP_PARAM)

  call LATLON_ico_setup

  call PROF_rapend('Initialize',0)
  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Main_MKLLMAP',0)

  call LATLON_setup( output_dir )

  call PROF_rapend('Main_MKLLMAP',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

  stop
end program prg_mkllmap
