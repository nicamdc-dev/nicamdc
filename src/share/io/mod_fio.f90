!-------------------------------------------------------------------------------
!> Module file I/O
!!
!! @par Description
!!         File I/O module (interface)
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_fio
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof

  use mod_fio_common
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FIO_setup
  public :: FIO_input
  public :: FIO_seek
  public :: FIO_output
  public :: FIO_close
  public :: FIO_finalize

  interface FIO_input
     module procedure FIO_input_SP
     module procedure FIO_input_DP
  end interface FIO_input

  interface FIO_output
     module procedure FIO_output_SP
     module procedure FIO_output_DP
  end interface FIO_output

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT),   private :: FIO_FORMAT = 'PANDA'

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine FIO_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_fio_panda, only: &
       FIO_PANDA_setup
    implicit none

    namelist / FIOPARAM / &
       FIO_FORMAT

    integer :: ierr
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[fileio]/Category[common share]'

    !--- read parameters
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=FIOPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** FIOPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,         *) 'xxx Not appropriate names in namelist FIOPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist FIOPARAM. STOP.'
       call PRC_MPIstop
    endif
    write(IO_FID_LOG,nml=FIOPARAM)

    select case(FIO_FORMAT)
    case('PANDA')
       write(IO_FID_LOG,*) '*** File I/O format : Packaged NICAM DAta format (PaNDa)'
       call FIO_PANDA_setup
    case('HDF5')
       write(IO_FID_LOG,*) '*** File I/O format : HDF5 format'
       write(*,*)           'xxx HDF5 format is not implemented yet! STOP.'
       call PRC_MPIstop
    case('NETCDF')
       write(IO_FID_LOG,*) '*** File I/O format : NetCDF(4) format'
       write(*,*)           'xxx NetCDF format is not implemented yet! STOP.'
       call PRC_MPIstop
    case default
       write(*,*) 'xxx Not appropriate FIO_FORMAT. STOP : ', trim(FIO_FORMAT)
       call PRC_MPIstop
    end select

    return
  end subroutine FIO_setup

  !-----------------------------------------------------------------------------
  subroutine FIO_input_SP( &
       var,            &
       basename,       &
       varname,        &
       layername,      &
       k_start,        &
       k_end,          &
       step,           &
       allow_missingq, & !--- optional
       did_out         ) !--- optional
    use mod_fio_panda, only: &
       FIO_PANDA_input_SP
    implicit none

    real(SP),         intent(out) :: var(:,:,:)
    character(len=*), intent(in)  :: basename
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: layername
    integer,          intent(in)  :: k_start, k_end
    integer,          intent(in)  :: step

    logical, intent(in),  optional :: allow_missingq !--- if data is missing, set value to zero
    integer, intent(out), optional :: did_out

    logical :: allow_missingq_
    integer :: did
    !---------------------------------------------------------------------------

    allow_missingq_ = .false.
    if ( present(allow_missingq) ) then
       allow_missingq_ = allow_missingq
    endif

    select case(FIO_FORMAT)
    case('PANDA')

       call FIO_PANDA_input_SP ( var(:,:,:),      & ! [OUT]
                                 basename,        & ! [IN]
                                 varname,         & ! [IN]
                                 layername,       & ! [IN]
                                 k_start,         & ! [IN]
                                 k_end,           & ! [IN]
                                 step,            & ! [IN]
                                 allow_missingq_, & ! [IN]
                                 did              ) ! [IN]

!   case('HDF5')
!
!      call FIO_HDF5_input_SP  ( var(:,:,:),      & ! [OUT]
!                                basename,        & ! [IN]
!                                varname,         & ! [IN]
!                                layername,       & ! [IN]
!                                k_start,         & ! [IN]
!                                k_end,           & ! [IN]
!                                step,            & ! [IN]
!                                allow_missingq_, & ! [IN]
!                                did              ) ! [IN]
!
!   case('NETCDF')
!
!      call FIO_NETCDF_input_SP( var(:,:,:),     & ! [OUT]
!                                basename,       & ! [IN]
!                                varname,        & ! [IN]
!                                layername,      & ! [IN]
!                                k_start,        & ! [IN]
!                                k_end,          & ! [IN]
!                                step,           & ! [IN]
!                                allow_missingq_ ) ! [IN]

    end select

    if ( present(did_out) )then
      did_out = did
    endif

    return
  end subroutine FIO_input_SP

  !-----------------------------------------------------------------------------
  subroutine FIO_input_DP( &
       var,            &
       basename,       &
       varname,        &
       layername,      &
       k_start,        &
       k_end,          &
       step,           &
       allow_missingq, & !--- optional
       did_out         )
    use mod_fio_panda, only: &
       FIO_PANDA_input_DP
    implicit none

    real(DP),         intent(out) :: var(:,:,:)
    character(len=*), intent(in)  :: basename
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: layername
    integer,          intent(in)  :: k_start, k_end
    integer,          intent(in)  :: step

    logical, intent(in),  optional :: allow_missingq !--- if data is missing, set value to zero
    integer, intent(out), optional :: did_out

    logical :: allow_missingq_
    integer :: did
    !---------------------------------------------------------------------------

    allow_missingq_ = .false.
    if ( present(allow_missingq) ) then
       allow_missingq_ = allow_missingq
    endif

    select case(FIO_FORMAT)
    case('PANDA')

       call FIO_PANDA_input_DP ( var(:,:,:),      & ! [OUT]
                                 basename,        & ! [IN]
                                 varname,         & ! [IN]
                                 layername,       & ! [IN]
                                 k_start,         & ! [IN]
                                 k_end,           & ! [IN]
                                 step,            & ! [IN]
                                 allow_missingq_, &
                                 did              ) ! [IN]

!   case('HDF5')
!
!      call FIO_HDF5_input_DP  ( var(:,:,:),      & ! [OUT]
!                                basename,        & ! [IN]
!                                varname,         & ! [IN]
!                                layername,       & ! [IN]
!                                k_start,         & ! [IN]
!                                k_end,           & ! [IN]
!                                step,            & ! [IN]
!                                allow_missingq_, & ! [IN]
!                                did              ) ! [IN]
!
!   case('NETCDF')
!
!      call FIO_NETCDF_input_DP( var(:,:,:),     & ! [OUT]
!                                basename,       & ! [IN]
!                                varname,        & ! [IN]
!                                layername,      & ! [IN]
!                                k_start,        & ! [IN]
!                                k_end,          & ! [IN]
!                                step,           & ! [IN]
!                                allow_missingq_ ) ! [IN]

    end select

    if ( present(did_out) )then
      did_out = did
    endif

    return
  end subroutine FIO_input_DP

  !-----------------------------------------------------------------------------
  subroutine FIO_seek( &
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
    use mod_fio_panda, only: &
       FIO_PANDA_seek
    implicit none

    integer,          intent(inout) :: start_step
    integer,          intent(inout) :: num_of_step
    integer,          intent(inout) :: data_date(6,FIO_data_nlim)
    integer,          intent(inout) :: prec
    character(len=*), intent(in)    :: basename
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: layername
    integer,          intent(in)    :: k_start, k_end
    real(DP),         intent(in)    :: ctime
    integer,          intent(in)    :: cdate(6)
    logical,          intent(in)    :: opt_periodic_year
    !---------------------------------------------------------------------------

    select case(FIO_FORMAT)
    case('PANDA')

       call FIO_PANDA_seek ( start_step,       & ! [INOUT]
                             num_of_step,      & ! [INOUT]
                             data_date(:,:),   & ! [INOUT]
                             prec,             & ! [INOUT]
                             basename,         & ! [IN]
                             varname,          & ! [IN]
                             layername,        & ! [IN]
                             k_start,          & ! [IN]
                             k_end,            & ! [IN]
                             ctime,            & ! [IN]
                             cdate(:),         & ! [IN]
                             opt_periodic_year ) ! [IN]

!   case('HDF5')
!
!      call FIO_HDF5_seek  ( start_step,       & ! [INOUT]
!                            num_of_step,      & ! [INOUT]
!                            data_date(:,:),   & ! [INOUT]
!                            prec,             & ! [INOUT]
!                            basename,         & ! [IN]
!                            varname,          & ! [IN]
!                            layername,        & ! [IN]
!                            k_start,          & ! [IN]
!                            k_end,            & ! [IN]
!                            ctime,            & ! [IN]
!                            cdate(:),         & ! [IN]
!                            opt_periodic_year ) ! [IN]
!
!   case('NETCDF')
!
!      call FIO_NETCDF_seek( start_step,       & ! [INOUT]
!                            num_of_step,      & ! [INOUT]
!                            data_date(:,:),   & ! [INOUT]
!                            prec,             & ! [INOUT]
!                            basename,         & ! [IN]
!                            varname,          & ! [IN]
!                            layername,        & ! [IN]
!                            k_start,          & ! [IN]
!                            k_end,            & ! [IN]
!                            ctime,            & ! [IN]
!                            cdate(:),         & ! [IN]
!                            opt_periodic_year ) ! [IN]

    end select

    return
  end subroutine FIO_seek

  !-----------------------------------------------------------------------------
  subroutine FIO_output_SP( &
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
       digits,    &
       abstol,    &
       reltol,    &
       append     )
    use mod_fio_panda, only: &
       FIO_PANDA_output_SP
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
    integer,          intent(in), optional :: digits ! [SCIL] significant digits       for lossy compression
    real(DP),         intent(in), optional :: abstol ! [SCIL] absolute error tolerance for lossy compression
    real(DP),         intent(in), optional :: reltol ! [SCIL] relative error tolerance for lossy compression
    logical,          intent(in), optional :: append

    integer  :: digits_
    real(DP) :: abstol_, reltol_
    logical  :: append_
    !---------------------------------------------------------------------------

    append_ = .false.
    if ( present(append) ) then
       append_ = append
    endif

    select case(FIO_FORMAT)
    case('PANDA')

       call FIO_PANDA_output_SP ( var(:,:,:), & ! [IN]
                                  basename,   & ! [IN]
                                  pkg_desc,   & ! [IN]
                                  pkg_note,   & ! [IN]
                                  varname,    & ! [IN]
                                  data_desc,  & ! [IN]
                                  data_note,  & ! [IN]
                                  unit,       & ! [IN]
                                  dtype,      & ! [IN]
                                  layername,  & ! [IN]
                                  k_start,    & ! [IN]
                                  k_end,      & ! [IN]
                                  step,       & ! [IN]
                                  t_start,    & ! [IN]
                                  t_end,      & ! [IN]
                                  append      ) ! [IN]

!   case('HDF5')
!
!      !=> [SCIL]
!      if ( present(digits) ) then
!         digits_ = digits
!      else
!         digits_ = -1
!      endif
!      if ( present(abstol) ) then
!         abstol_ = abstol
!      else
!         abstol_ = -999.0_DP
!      endif
!      if ( present(reltol) ) then
!         reltol_ = reltol
!      else
!         reltol_ = -999.0_DP
!      endif
!      !<= [SCIL]
!
!      call FIO_HDF5_output_SP  ( var(:,:,:), & ! [IN]
!                                 basename,   & ! [IN]
!                                 pkg_desc,   & ! [IN]
!                                 pkg_note,   & ! [IN]
!                                 varname,    & ! [IN]
!                                 data_desc,  & ! [IN]
!                                 data_note,  & ! [IN]
!                                 unit,       & ! [IN]
!                                 dtype,      & ! [IN]
!                                 layername,  & ! [IN]
!                                 k_start,    & ! [IN]
!                                 k_end,      & ! [IN]
!                                 step,       & ! [IN]
!                                 t_start,    & ! [IN]
!                                 t_end,      & ! [IN]
!                                 digits_,    & ! [IN] [SCIL]
!                                 abstol_,    & ! [IN] [SCIL]
!                                 reltol_,    & ! [IN] [SCIL]
!                                 append      ) ! [IN]
!
!   case('NETCDF')
!
!      call FIO_NETCDF_output_SP( var(:,:,:), & ! [IN]
!                                 basename,   & ! [IN]
!                                 pkg_desc,   & ! [IN]
!                                 pkg_note,   & ! [IN]
!                                 varname,    & ! [IN]
!                                 data_desc,  & ! [IN]
!                                 data_note,  & ! [IN]
!                                 unit,       & ! [IN]
!                                 dtype,      & ! [IN]
!                                 layername,  & ! [IN]
!                                 k_start,    & ! [IN]
!                                 k_end,      & ! [IN]
!                                 step,       & ! [IN]
!                                 t_start,    & ! [IN]
!                                 t_end,      & ! [IN]
!                                 append      ) ! [IN]

    end select

    return
  end subroutine FIO_output_SP

  !-----------------------------------------------------------------------------
  subroutine FIO_output_DP( &
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
       digits,    &
       abstol,    &
       reltol,    &
       append     )
    use mod_fio_panda, only: &
       FIO_PANDA_output_DP
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
    integer,          intent(in), optional :: digits ! [SCIL] significant digits       for lossy compression
    real(DP),         intent(in), optional :: abstol ! [SCIL] absolute error tolerance for lossy compression
    real(DP),         intent(in), optional :: reltol ! [SCIL] relative error tolerance for lossy compression
    logical,          intent(in), optional :: append

    integer  :: digits_
    real(DP) :: abstol_, reltol_
    logical  :: append_
    !---------------------------------------------------------------------------

    append_ = .false.
    if ( present(append) ) then
       append_ = append
    endif

    select case(FIO_FORMAT)
    case('PANDA')

       call FIO_PANDA_output_DP ( var(:,:,:), & ! [IN]
                                  basename,   & ! [IN]
                                  pkg_desc,   & ! [IN]
                                  pkg_note,   & ! [IN]
                                  varname,    & ! [IN]
                                  data_desc,  & ! [IN]
                                  data_note,  & ! [IN]
                                  unit,       & ! [IN]
                                  dtype,      & ! [IN]
                                  layername,  & ! [IN]
                                  k_start,    & ! [IN]
                                  k_end,      & ! [IN]
                                  step,       & ! [IN]
                                  t_start,    & ! [IN]
                                  t_end,      & ! [IN]
                                  append      ) ! [IN]

!   case('HDF5')
!
!      !=> [SCIL]
!      if ( present(digits) ) then
!         digits_ = digits
!      else
!         digits_ = -1
!      endif
!      if ( present(abstol) ) then
!         abstol_ = abstol
!      else
!         abstol_ = -999.0_DP
!      endif
!      if ( present(reltol) ) then
!         reltol_ = reltol
!      else
!         reltol_ = -999.0_DP
!      endif
!      !<= [SCIL]
!
!      call FIO_HDF5_output_DP  ( var(:,:,:), & ! [IN]
!                                 basename,   & ! [IN]
!                                 pkg_desc,   & ! [IN]
!                                 pkg_note,   & ! [IN]
!                                 varname,    & ! [IN]
!                                 data_desc,  & ! [IN]
!                                 data_note,  & ! [IN]
!                                 unit,       & ! [IN]
!                                 dtype,      & ! [IN]
!                                 layername,  & ! [IN]
!                                 k_start,    & ! [IN]
!                                 k_end,      & ! [IN]
!                                 step,       & ! [IN]
!                                 t_start,    & ! [IN]
!                                 t_end,      & ! [IN]
!                                 digits_,    & ! [IN] [SCIL]
!                                 abstol_,    & ! [IN] [SCIL]
!                                 reltol_,    & ! [IN] [SCIL]
!                                 append      ) ! [IN]
!
!   case('NETCDF')
!
!      call FIO_NETCDF_output_DP( var(:,:,:), & ! [IN]
!                                 basename,   & ! [IN]
!                                 pkg_desc,   & ! [IN]
!                                 pkg_note,   & ! [IN]
!                                 varname,    & ! [IN]
!                                 data_desc,  & ! [IN]
!                                 data_note,  & ! [IN]
!                                 unit,       & ! [IN]
!                                 dtype,      & ! [IN]
!                                 layername,  & ! [IN]
!                                 k_start,    & ! [IN]
!                                 k_end,      & ! [IN]
!                                 step,       & ! [IN]
!                                 t_start,    & ! [IN]
!                                 t_end,      & ! [IN]
!                                 append      ) ! [IN]

    end select

    return
  end subroutine FIO_output_DP

  !-----------------------------------------------------------------------------
  subroutine FIO_close( &
       basename )
    use mod_fio_panda, only: &
       FIO_PANDA_close
    implicit none

    character(len=*), intent(in) :: basename
    !---------------------------------------------------------------------------

    select case(FIO_FORMAT)
    case('PANDA')

       call FIO_PANDA_close ( basename ) ! [IN]

!   case('HDF5')
!
!      call FIO_HDF5_close  ( basename ) ! [IN]
!
!   case('NETCDF')
!
!      call FIO_NETCDF_close( basename ) ! [IN]

    end select

    return
  end subroutine FIO_close

  !-----------------------------------------------------------------------------
  subroutine FIO_finalize
    use mod_fio_panda, only: &
       FIO_PANDA_finalize
    implicit none
    !---------------------------------------------------------------------------

    select case(FIO_FORMAT)
    case('PANDA')

       call FIO_PANDA_finalize

!   case('HDF5')
!
!      call FIO_HDF5_finalize
!
!   case('NETCDF')
!
!      call FIO_NETCDF_finalize

    end select

    return
  end subroutine FIO_finalize

end module mod_fio
