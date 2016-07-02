!-------------------------------------------------------------------------------
!
!+  Program mkllmap
!
!-------------------------------------------------------------------------------
program prg_mkllmap
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !
  !++ Current Corresponding Author :
  !
  !++ History:
  !       This program originate from 'cnvlatlon.f90'(by H.Tomita) ver 4.38 .
  !
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !                 11-11-09  H.Yashiro [mod] Avoid arc-cos, precise calculation
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules (shared)
  !
  use mod_precision
  use mod_stdio
  use mod_debug
  use mod_process, only: &
     PRC_MPIstart,    &
     PRC_LOCAL_setup, &
     PRC_MPIstop,     &
     PRC_MPIfinish
  use mod_adm, only: &
     ADM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_cnst, only: &
     CNST_setup
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
  integer :: comm_world
  integer :: myrank
  logical :: ismaster

  character(len=H_LONG) :: output_dir   = './'

  namelist /mkllmap_param/ &
     output_dir

  integer :: ierr
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

  !--- < admin module setup > ---
  call ADM_setup
  !
  call FIO_setup
  !
  !--- < comm module setup > ---
  call COMM_setup
  !
  !--- < cnst module setup > ---
  call CNST_setup
  !
  !--- < grid module setup > ---
  call GRD_setup

  !--- read parameters
  write(IO_FID_LOG,*)
  write(IO_FID_LOG,*) '+++ Program[mkllmap]/Category[tool]'
  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=MKLLMAP_PARAM,iostat=ierr)
  if ( ierr < 0 ) then
     write(IO_FID_LOG,*) '*** MKLLMAP_PARAM is not specified. use default.'
  elseif( ierr > 0 ) then
     write(*         ,*) 'xxx Not appropriate names in namelist MKLLMAP_PARAM. STOP.'
     write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist MKLLMAP_PARAM. STOP.'
     call PRC_MPIstop
  endif
  write(IO_FID_LOG,nml=MKLLMAP_PARAM)

  call LATLON_ico_setup

  call LATLON_setup( output_dir )

  call PRC_MPIfinish

end program prg_mkllmap
!-------------------------------------------------------------------------------
