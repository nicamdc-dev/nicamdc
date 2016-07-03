!-------------------------------------------------------------------------------
!>
!! File I/O module
!!
!! @par Description
!!         This module is continer for file I/O (PaNDa on HDF5, POH5 format)
!!
!! @author T.Inoue
!!
!! @par History
!! @li      2016-02-19 (T.Inoue)   [NEW]
!!
!<
module mod_hio
  !
  !++ Used modules
  !
  use mod_precision
  use mod_io_param
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: HIO_setup
  public :: HIO_input
  public :: HIO_seek
  public :: HIO_output
  public :: HIO_close
  public :: HIO_finalize

  interface HIO_input
     module procedure HIO_input_SP
     module procedure HIO_input_DP
  end interface HIO_input

  interface HIO_output
     module procedure HIO_output_SP
     module procedure HIO_output_DP
  end interface HIO_output

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !> struct for package infomation
  type, public :: headerinfo
     character(len=H_LONG)  :: fname         !< file name
     character(len=H_MID)   :: description   !< file description
     character(len=H_LONG)  :: note          !< longer note of file
     integer                :: num_of_var    !< number of data
     integer                :: fmode         !< file mode(0,1,2)
     integer                :: endiantype    !< endian type(0,1,2)
     integer                :: grid_topology !< grid topology(0,1,2)
     integer                :: glevel        !< glevel
     integer                :: rlevel        !< rlevel
     integer                :: num_of_rgn    !< number of region
     integer, pointer       :: rgnid(:)      !< array of region id
  endtype headerinfo

  !> struct for data infomation
  type, public :: datainfo
     character(len=H_SHORT) :: varname      !< variable name
     character(len=H_MID)   :: description  !< variable description
     character(len=H_SHORT) :: unit         !< unit of variable
     character(len=H_SHORT) :: layername    !< layer name
     character(len=H_LONG)  :: note         !< longer note of variable
     integer(DP)            :: datasize     !< data size
     integer                :: datatype     !< data type(0,1,2,3)
     integer                :: num_of_layer !< number of layer
     integer                :: num_of_step  !< number of step
  endtype datainfo

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: max_num_of_data = 2500 !--- max time step num

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup this module.
  !!
  !! Must be called first.
  !!
  subroutine HIO_setup
    implicit none
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[hio]/Category[common share]'
    write(IO_FID_LOG,*) '+++ HIO is disabled.'

    return
  end subroutine HIO_setup

  !-----------------------------------------------------------------------------
  subroutine HIO_input_SP( &
       var,           &
       basename,      &
       varname,       &
       layername,     &
       k_start,       &
       k_end,         &
       step,          &
       allow_missingq ) !--- optional
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(SP),         intent(out) :: var(:,:,:) !< variable(ij,k,l)
    character(len=*), intent(in)  :: basename   !< basename of file
    character(len=*), intent(in)  :: varname    !< variable name
    character(len=*), intent(in)  :: layername  !< layer name
    integer,          intent(in)  :: k_start    !< start index of vertical level
    integer,          intent(in)  :: k_end      !< end   index of vertical level
    integer,          intent(in)  :: step       !< step to be read

    logical, intent(in), optional :: allow_missingq !< if data is missing, set value to zero, else execution stops.
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) 'xxx HIO is disabled. STOP!', trim(varname), trim(basename)
    call PRC_MPIstop

    return
  end subroutine HIO_input_SP

  !-----------------------------------------------------------------------------
  subroutine HIO_input_DP( &
       var,           &
       basename,      &
       varname,       &
       layername,     &
       k_start,       &
       k_end,         &
       step,          &
       allow_missingq ) !--- optional
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(DP),         intent(out) :: var(:,:,:) !< variable(ij,k,l)
    character(len=*), intent(in)  :: basename   !< basename of file
    character(len=*), intent(in)  :: varname    !< variable name
    character(len=*), intent(in)  :: layername  !< layer name
    integer,          intent(in)  :: k_start    !< start index of vertical level
    integer,          intent(in)  :: k_end      !< end   index of vertical level
    integer,          intent(in)  :: step       !< step to be read

    logical, intent(in), optional :: allow_missingq !< if data is missing, set value to zero, else execution stops.
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) 'xxx HIO is disabled. STOP!', trim(varname), trim(basename)
    call PRC_MPIstop

    return
  end subroutine HIO_input_DP

  !-----------------------------------------------------------------------------
  subroutine HIO_seek( &
       start_step,       &
       num_of_step,      &
       data_date,        &
       prec,             &
       basename,         &
       varname,          &
       layername,        &
       k_start,          &
       k_end,            &
       ctime,            &
       cdate,            &
       opt_periodic_year )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    integer,          intent(inout) :: start_step
    integer,          intent(inout) :: num_of_step
    integer,          intent(inout) :: data_date(6,max_num_of_data)
    integer,          intent(inout) :: prec
    character(len=*), intent(in)    :: basename
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: layername      ! for verification only
    integer,          intent(in)    :: k_start, k_end ! for verification only
    real(DP),         intent(in)    :: ctime
    integer,          intent(in)    :: cdate(6)       ! cdate(1) is only used only when opt_periodic_year is T.
    logical,          intent(in)    :: opt_periodic_year
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) 'xxx HIO is disabled. STOP!', trim(varname), trim(basename)
    call PRC_MPIstop

    return
  end subroutine HIO_seek

  !-----------------------------------------------------------------------------
  !> Append data with data header
  subroutine HIO_output_SP( &
       var,       &
       basename,  &
       pkg_desc,  &
       pkg_note,  &
       varname,   &
       data_desc, &
       data_note, &
       unit,      &
       dtype,     &
       layername, &
       k_start,   &
       k_end,     &
       step,      &
       t_start,   &
       t_end,     &
       append     )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(SP),         intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: basename
    character(len=*), intent(in) :: pkg_desc
    character(len=*), intent(in) :: pkg_note
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: data_desc
    character(len=*), intent(in) :: data_note
    character(len=*), intent(in) :: unit
    integer,          intent(in) :: dtype
    character(len=*), intent(in) :: layername
    integer,          intent(in) :: k_start, k_end
    integer,          intent(in) :: step
    real(DP),         intent(in) :: t_start, t_end

    logical,intent(in), optional :: append
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) 'xxx HIO is disabled. STOP!', trim(varname), trim(basename)
    call PRC_MPIstop

    return
  end subroutine HIO_output_SP

  !-----------------------------------------------------------------------------
  !> Append data with data header
  subroutine HIO_output_DP( &
       var,       &
       basename,  &
       pkg_desc,  &
       pkg_note,  &
       varname,   &
       data_desc, &
       data_note, &
       unit,      &
       dtype,     &
       layername, &
       k_start,   &
       k_end,     &
       step,      &
       t_start,   &
       t_end,     &
       append     )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(DP),         intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: basename
    character(len=*), intent(in) :: pkg_desc
    character(len=*), intent(in) :: pkg_note
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: data_desc
    character(len=*), intent(in) :: data_note
    character(len=*), intent(in) :: unit
    integer,          intent(in) :: dtype
    character(len=*), intent(in) :: layername
    integer,          intent(in) :: k_start, k_end
    integer,          intent(in) :: step
    real(DP),         intent(in) :: t_start, t_end

    logical,intent(in), optional :: append
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) 'xxx HIO is disabled. STOP!', trim(varname), trim(basename)
    call PRC_MPIstop

    return
  end subroutine HIO_output_DP

  !-----------------------------------------------------------------------------
  subroutine HIO_close( &
       basename )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: basename
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) 'xxx HIO is disabled. STOP!', trim(basename)
    call PRC_MPIstop

    return
  end subroutine HIO_close

  !-----------------------------------------------------------------------------
  subroutine HIO_finalize
    use mod_process, only: &
       PRC_MPIstop
    implicit none
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) 'xxx HIO is disabled. STOP!'
    call PRC_MPIstop

    return
  end subroutine HIO_finalize

end module mod_hio
