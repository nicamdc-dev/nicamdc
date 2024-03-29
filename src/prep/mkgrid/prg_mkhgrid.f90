!-------------------------------------------------------------------------------
!> Program mkhgrid
!!
!! @par Description
!!          Making horizontal grid systems based on the icosahedral grid configuration
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
program mkhgrid
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
     PRC_MPIfinish
  use mod_const, only: &
     RADIUS => CONST_RADIUS, &
     CONST_setup
  use mod_adm, only: &
     ADM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_grd, only: &
     GRD_input_hgrid,  &
     GRD_output_hgrid, &
     GRD_makelatlon,   &
     GRD_scaling
  use mod_gmtr, only: &
     GMTR_setup
  use mod_mkgrd, only: &
     MKGRD_setup,        &
     MKGRD_prerotate,    &
     MKGRD_stretch,      &
     MKGRD_shrink,       &
     MKGRD_rotate,       &
     MKGRD_gravcenter,   &
     MKGRD_IN_BASENAME,  &
     MKGRD_IN_io_mode,   &
     MKGRD_OUT_BASENAME, &
     MKGRD_OUT_io_mode
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  integer  :: comm_world
  integer  :: myrank
  logical  :: ismaster

  !=============================================================================

  !---< MPI start >---
  call PRC_MPIstart( comm_world ) ! [OUT]

  !---< STDIO setup >---
  call IO_setup( 'NICAM-DC',   & ! [IN]
                 'mkhgrid.cnf' ) ! [IN]

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

  !---< cnst module setup >---
  call CONST_setup

  !---< admin module setup >---
  call ADM_setup

  !---< I/O module setup >---
  call FIO_setup

  !---< comm module setup >---
  call COMM_setup

  !---< mkgrid module setup >---
  call MKGRD_setup

  call PROF_rapend('Initialize',0)
  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Main_MKGRD',0)

  call GRD_input_hgrid( basename     = MKGRD_IN_BASENAME, & ! [IN]
                        input_vertex = .false.,           & ! [IN]
                        io_mode      = MKGRD_IN_io_mode   ) ! [IN]

  call PROF_rapstart('MKGRD_prerotate',0)
  call MKGRD_prerotate
  call PROF_rapend  ('MKGRD_prerotate',0)

  call PROF_rapstart('MKGRD_stretch',0)
  call MKGRD_stretch
  call PROF_rapend  ('MKGRD_stretch',0)

  call PROF_rapstart('MKGRD_shrink',0)
  call MKGRD_shrink
  call PROF_rapend  ('MKGRD_shrink',0)

  call PROF_rapstart('MKGRD_rotate',0)
  call MKGRD_rotate
  call PROF_rapend  ('MKGRD_rotate',0)

  call PROF_rapstart('MKGRD_gravcenter',0)
  call MKGRD_gravcenter
  call PROF_rapend  ('MKGRD_gravcenter',0)

  call GRD_output_hgrid( basename      = MKGRD_OUT_BASENAME, & ! [IN]
                         output_vertex = .true.,             & ! [IN]
                         io_mode       = MKGRD_OUT_io_mode   ) ! [IN]


  !---< gmtr module setup >---
  call PROF_rapstart('GMTR_setup',0)
  call GRD_makelatlon
  call GRD_scaling( RADIUS )
  call GMTR_setup
  call PROF_rapend  ('GMTR_setup',0)

  call PROF_rapend('Main_MKGRD',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

  stop
end program mkhgrid
