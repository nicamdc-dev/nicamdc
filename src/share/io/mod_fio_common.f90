!-------------------------------------------------------------------------------
!> Module file I/O common
!!
!! @par Description
!!         File I/O module (common parameter)
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_fio_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use iso_c_binding

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: cstr
  public :: cstr2
  public :: fstr

  interface fstr
     module procedure fstr1
     module procedure fstr2
  end interface fstr

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !--- character length
  integer, public, parameter :: FIO_HSHORT =   16 !< character length for short var.
  integer, public, parameter :: FIO_HMID   =   64 !< character length for middle var.
  integer, public, parameter :: FIO_HLONG  =  256 !< character length for long var.
  integer, public, parameter :: FIO_HFILE  = 1024 !< character length for long var.

  !--- data type
  integer, public, parameter :: FIO_REAL4    = 0 !< ID for 4byte real
  integer, public, parameter :: FIO_REAL8    = 1 !< ID for 8byte real
  integer, public, parameter :: FIO_INTEGER4 = 2 !< ID for 4byte int
  integer, public, parameter :: FIO_INTEGER2 = 3 !< ID for 8byte int

  !--- data endian
  integer, public, parameter :: FIO_UNKNOWN_ENDIAN = 0 !< ID for unknown endian
  integer, public, parameter :: FIO_LITTLE_ENDIAN  = 1 !< ID for little endian
  integer, public, parameter :: FIO_BIG_ENDIAN     = 2 !< ID for big endian

  !--- topology
  integer, public, parameter :: FIO_ICOSAHEDRON = 0 !< ID for ico grid
  integer, public, parameter :: FIO_IGA_LCP     = 1 !< ID for LCP grid
  integer, public, parameter :: FIO_IGA_MLCP    = 2 !< ID for MLCP grid

  !--- file mode (partial or complete)
  integer, public, parameter :: FIO_SPLIT_FILE = 0 !< ID for split(partical) file
  integer, public, parameter :: FIO_INTEG_FILE = 1 !< ID for integrated(complete) file

  !--- proccessor type
  integer, public, parameter :: FIO_SINGLE_PROC = 0 !< ID for single processor
  integer, public, parameter :: FIO_MULTI_PROC  = 1 !< ID for multi processor

  !--- action type
  integer, public, parameter :: FIO_FREAD    =  0 !< ID for read file
  integer, public, parameter :: FIO_FWRITE   =  1 !< ID for write file
  integer, public, parameter :: FIO_FAPPEND  =  2 !< ID for append file

  !--- data dump type
  integer, public, parameter :: FIO_DUMP_OFF      = 0 !< Dumping off
  integer, public, parameter :: FIO_DUMP_HEADER   = 1 !< Dump header only
  integer, public, parameter :: FIO_DUMP_ALL      = 2 !< Dump all
  integer, public, parameter :: FIO_DUMP_ALL_MORE = 3 !< Dump all and more

  integer, public, parameter :: FIO_file_nlim     = 1024
  integer, public, parameter :: FIO_data_nlim     = 2500 !< max time step num
  integer, public, parameter :: FIO_preclist(0:3) = (/ 4, 8, 4, 2 /)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  function cstr(str)
    implicit none

    character(*), intent(in) :: str

    character(:,c_char), allocatable, target :: cstr
    !---------------------------------------------------------------------------

    cstr = trim(str)//c_null_char

  end function cstr

  !-----------------------------------------------------------------------------
  subroutine cstr2(cstr, fstr)
    implicit none

    character(c_char), intent(out) :: cstr(:)
    character(len=*),  intent(in)  :: fstr

    integer :: i, j, l
    !---------------------------------------------------------------------------

    l = size(cstr)

    cstr(l) = c_null_char ! at least last digit

    do i = l-1, 1, -1
       if ( i > len(fstr) ) then
          cstr(i) = c_null_char
       else
          if ( fstr(i:i) == " " ) then
             cstr(i) = c_null_char
          else
             exit
          endif
       endif
    enddo

    do j = i, 1, -1
       cstr(j) = fstr(j:j)
    enddo

    return
  end subroutine cstr2

  !-----------------------------------------------------------------------------
  subroutine fstr1(str)
    implicit none

    character(len=*), intent(inout) :: str

    integer :: i, j
    !---------------------------------------------------------------------------

    do i = 1, len(str)
       if( str(i:i) == c_null_char ) exit
    enddo

    do j = i, len(str)
       str(j:j) = " "
    enddo

    return
  end subroutine fstr1

  !-----------------------------------------------------------------------------
  subroutine fstr2(fstr, cstr)
    implicit none

    character(len=*),  intent(out) :: fstr
    character(c_char), intent(in)  :: cstr(:)

    integer :: i, j, l
    !---------------------------------------------------------------------------

    l = min( len(fstr), size(cstr) )

    do i = 1, l
       if( cstr(i) == c_null_char ) exit
       fstr(i:i) = cstr(i)
    enddo

    do j = i, len(fstr)
       fstr(j:j) = " "
    enddo

    return
  end subroutine fstr2

end module mod_fio_common
