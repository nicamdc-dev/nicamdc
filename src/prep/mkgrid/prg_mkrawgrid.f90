!-------------------------------------------------------------------------------
!>
!! Program mkrawgrid
!!
!! @par Description
!!          Making grid systems based on the icosahedral grid configuration
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2013-05-1  (H.Yashiro) NICAM-DC
!!
!<
program mkrawgrid
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
  use mod_adm, only: &
     ADM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_const, only: &
     CONST_setup
  use mod_grd, only: &
     GRD_output_hgrid
  use mod_mkgrd, only: &
     MKGRD_setup,        &
     MKGRD_standard,     &
     MKGRD_spring,       &
     MKGRD_OUT_BASENAME, &
     MKGRD_OUT_io_mode
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  integer :: comm_world
  integer :: myrank
  logical :: ismaster

  !=============================================================================

  !---< MPI start >---
  call PRC_MPIstart( comm_world ) ! [OUT]

  !---< STDIO setup >---
  call IO_setup( 'NICAM-DC',     & ! [IN]
                 'mkrawgrid.cnf' ) ! [IN]

  !---< Local process management setup >---
  call PRC_LOCAL_setup( comm_world, & ! [IN]
                        myrank,     & ! [OUT]
                        ismaster    ) ! [OUT]

  !---< Logfile setup >---
  call IO_LOG_setup( myrank,  & ! [IN]
                     ismaster ) ! [IN]

  !---< admin module setup >---
  call ADM_setup

  !---< I/O module setup >---
  call FIO_setup

  !---< comm module setup >---
  call COMM_setup

  !---< cnst module setup >---
  call CONST_setup

  !---< mkgrid module setup >---
  call MKGRD_setup

  !########## main ##########

  call MKGRD_standard

  call MKGRD_spring

  call GRD_output_hgrid( basename      = MKGRD_OUT_BASENAME, & ! [IN]
                         output_vertex = .false.,            & ! [IN]
                         io_mode       = MKGRD_OUT_io_mode   ) ! [IN]

  !########## Finalize ##########

  !--- all processes stop
  call PRC_MPIfinish

  stop
end program mkrawgrid
!-------------------------------------------------------------------------------
