!-------------------------------------------------------------------------------
!> module PROCESS
!!
!! @par Description
!!          MPI/non-MPI management module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-11 (R.Yoshida)  [new]
!! @li      2011-11-11 (H.Yashiro)  [mod] Integrate to SCALE-LES ver.3
!!
!<
module mod_process
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use mod_precision
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PRC_MPIstart
  public :: PRC_LOCAL_setup
  public :: PRC_MPIstop
  public :: PRC_MPIfinish

  public :: PRC_MPIbarrier
  public :: PRC_MPItime
  public :: PRC_MPItimestat

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !                          [ communicator system ]
  !    MPI_COMM_WORLD
  !          |
  ! PRC_LOCAL_COMM_WORLD --split--> BULK_COMM_WORLD
  !                                     |
  !                            PRC_GLOBAL_COMM_WORLD --split--> PRC_LOCAL_COMM_WORLD
  !-----------------------------------------------------------------------------
  integer, public, parameter :: PRC_masterrank      = 0    !< master process in each communicator

  ! local world
  integer, public :: PRC_LOCAL_COMM_WORLD     = -1      !< local communicator
  integer, public :: PRC_nprocs               = 1       !< myrank         in local communicator
  integer, public :: PRC_myrank               = 0       !< process num    in local communicator
  logical, public :: PRC_IsMaster             = .false. !< master process in local communicator?
  logical, public :: PRC_mpi_alive            = .false. !< MPI is alive?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Start MPI
  subroutine PRC_MPIstart( &
       comm )
    implicit none

    integer, intent(out) :: comm ! communicator

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_Init(ierr)

    PRC_mpi_alive = .true.

    comm = MPI_COMM_WORLD

    return
  end subroutine PRC_MPIstart

  !-----------------------------------------------------------------------------
  !> Setup MPI
  subroutine PRC_LOCAL_setup( &
       comm,    &
       myrank,  &
       ismaster )
    implicit none

    integer, intent(in)  :: comm     ! communicator
    integer, intent(out) :: myrank   ! myrank         in this communicator
    logical, intent(out) :: ismaster ! master process in this communicator?

    integer :: ierr
    !---------------------------------------------------------------------------

    PRC_LOCAL_COMM_WORLD = comm

    call MPI_COMM_RANK(PRC_LOCAL_COMM_WORLD,PRC_myrank,ierr)
    call MPI_COMM_SIZE(PRC_LOCAL_COMM_WORLD,PRC_nprocs,ierr)

    if ( PRC_myrank == PRC_masterrank ) then
       PRC_IsMaster = .true.
    else
       PRC_IsMaster = .false.
    endif

    myrank   = PRC_myrank
    ismaster = PRC_IsMaster

    return
  end subroutine PRC_LOCAL_setup

  !-----------------------------------------------------------------------------
  !> Abort MPI
  subroutine PRC_MPIstop
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    ! flush 1kbyte
    write(IO_FID_LOG,'(32A32)') '                                '

    write(IO_FID_LOG,*)           '+++ Abort MPI'
    if( PRC_IsMaster ) write(*,*) '+++ Abort MPI'

    if ( IO_L ) then
       if( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    if ( PRC_mpi_alive ) then
       ! Abort MPI
       call MPI_Abort(PRC_LOCAL_COMM_WORLD, 1, ierr)
    endif

    stop
  end subroutine PRC_MPIstop

  !-----------------------------------------------------------------------------
  !> Stop MPI peacefully
  subroutine PRC_MPIfinish
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    ! Stop MPI
    if ( PRC_mpi_alive ) then
       if ( IO_L ) then
          write(IO_FID_LOG,*)
          write(IO_FID_LOG,*) '++++++ Stop MPI'
          write(IO_FID_LOG,*)
       endif

       call MPI_Barrier(PRC_LOCAL_COMM_WORLD,ierr)

       call MPI_Finalize(ierr)
       if( IO_L ) write(IO_FID_LOG,*) '*** MPI is peacefully finalized'
    endif

    ! Close logfile, configfile
    if ( IO_L ) then
       if( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    ! Stop program
    stop
  end subroutine PRC_MPIfinish

  !-----------------------------------------------------------------------------
  !> Barrier MPI
  subroutine PRC_MPIbarrier
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( PRC_mpi_alive ) then
       call MPI_barrier(PRC_LOCAL_COMM_WORLD,ierr)
    endif

  end subroutine PRC_MPIbarrier

  !-----------------------------------------------------------------------------
  !> Get MPI time
  !> @return time
  function PRC_MPItime() result(time)
    implicit none

    real(DP) :: time
    !---------------------------------------------------------------------------

    if ( PRC_mpi_alive ) then
       time = real(MPI_WTIME(), kind=DP)
    else
       call CPU_TIME(time)
    endif

  end function PRC_MPItime

  !-----------------------------------------------------------------------------
  !> Calc global statistics for timer
  subroutine PRC_MPItimestat( &
      avgvar, &
      maxvar, &
      minvar, &
      maxidx, &
      minidx, &
      var     )
    implicit none

    real(DP), intent(out) :: avgvar(:) !< average
    real(DP), intent(out) :: maxvar(:) !< maximum
    real(DP), intent(out) :: minvar(:) !< minimum
    integer,  intent(out) :: maxidx(:) !< index of maximum
    integer,  intent(out) :: minidx(:) !< index of minimum
    real(DP), intent(in)  :: var(:)    !< values for statistics

    real(DP), allocatable :: statval(:,:)
    integer               :: vsize

    real(DP) :: totalvar
    integer  :: ierr
    integer  :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:))

    allocate( statval(vsize,0:PRC_nprocs-1) )
    statval(:,:) = 0.0_DP

    do v = 1, vsize
       statval(v,PRC_myrank) = var(v)
    enddo

    ! MPI broadcast
    do p = 0, PRC_nprocs-1
       call MPI_Bcast( statval(1,p),         &
                       vsize,                &
                       MPI_DOUBLE_PRECISION, &
                       p,                    &
                       PRC_LOCAL_COMM_WORLD, &
                       ierr                  )
    enddo

    do v = 1, vsize
       totalvar = 0.0_DP
       do p = 0, PRC_nprocs-1
          totalvar = totalvar + statval(v,p)
       enddo
       avgvar(v) = totalvar / PRC_nprocs

       maxvar(v)   = maxval(statval(v,:))
       minvar(v)   = minval(statval(v,:))
       maxidx(v:v) = maxloc(statval(v,:))
       minidx(v:v) = minloc(statval(v,:))
    enddo

    deallocate( statval )

    return
  end subroutine PRC_MPItimestat

end module mod_process
!-------------------------------------------------------------------------------
