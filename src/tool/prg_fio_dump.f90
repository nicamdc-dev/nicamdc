!-------------------------------------------------------------------------------
!> Program FIO dump
!!
!! @par Description
!!          header/data veiwer for new format data
!!          ( packaged NICAM data format : PaNDa )
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
program prg_fio_dump
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use iso_c_binding
  use mod_fio_common
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ C interface
  !
  include 'mod_fio_panda.inc'

  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  character(len=H_LONG) :: fname     = ''
  integer               :: mode      = FIO_DUMP_HEADER
  integer               :: endian    = FIO_BIG_ENDIAN
  logical               :: filelok   = .false.
  logical               :: modelok   = .false.
  logical               :: endianlok = .false.

  character(len=H_LONG) :: argstr
  integer  :: n, narg

#ifdef _NOF2003
  integer  :: IARGC
#else
  integer  :: command_argument_count
#endif

  integer  :: fid, ierr ! return from C program
  !=============================================================================

  call MPI_Init(ierr)

#ifdef _NOF2003
  narg = IARGC()
#else
  narg = command_argument_count()
#endif

  if ( narg == 0 ) then
     write(*,*) 'Usage : fio_dump [option] [file]'
     write(*,*) '  -h show header only'
     write(*,*) '  -d dump all data   '
     write(*,*) '  -e dump all data (60 digit mode)'
     write(*,*) '  -b force dump with big-endian'
     write(*,*) '  -l force dump with little-endian'
     stop
  endif

  do n = 1, narg

#ifdef _NOF2003
     call GETARG(n,argstr)
#else
     call get_command_argument(n,argstr)
#endif

     if ( argstr(1:1) == '-' ) then
        select case(argstr(2:2))
        case('h')
           if(.not. modelok) mode = FIO_DUMP_HEADER
           modelok = .true.
        case('d')
           if(.not. modelok) mode = FIO_DUMP_ALL
           modelok = .true.
        case('e') ! [add] 20120621 H.Yashiro
           if(.not. modelok) mode = FIO_DUMP_ALL_MORE
           modelok = .true.
        case('b')
           if(.not. endianlok) endian = FIO_BIG_ENDIAN
           endianlok = .true.
        case('l')
           if(.not. endianlok) endian = FIO_LITTLE_ENDIAN
           endianlok = .true.
        endselect
     else
        if(.not. filelok) fname = trim(argstr)
        filelok = .true.
     endif
  enddo

  ierr = fio_syscheck()

  fid = fio_register_file(cstr(fname))

  ierr = fio_dump_finfo(fid,endian,mode)

  call MPI_Finalize(ierr)

end program prg_fio_dump
