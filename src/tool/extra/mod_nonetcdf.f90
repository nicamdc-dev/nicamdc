! dummy module

module mod_netcdf
  use mod_precision
  implicit none
  private
  !
  !++ Private parameters
  !
  integer,parameter :: NOT_OPENED     = 0
  integer,parameter :: OPEN_FOR_READ  = 1
  integer,parameter :: OPEN_FOR_WRITE = -1
  integer,parameter :: CLEN           = 1024
  !
  !++ Public procedures
  !
  public :: netcdf_set_logfid
  public :: netcdf_open_for_write
  public :: netcdf_open_for_read
  public :: netcdf_write
  public :: netcdf_read
  public :: netcdf_read_dim
  public :: netcdf_close
  public :: netcdf_create_grads_ctl

  !--- interfaces for multiple array dimensions
  interface netcdf_write
     module procedure netcdf_write_0d
     module procedure netcdf_write_1d
     module procedure netcdf_write_2d
     module procedure netcdf_write_3d
     module procedure netcdf_write_4d
  end interface

  interface netcdf_read
     module procedure netcdf_read_0d
     module procedure netcdf_read_1d
     module procedure netcdf_read_2d
     module procedure netcdf_read_3d
     module procedure netcdf_read_4d
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !--- Handler for each netCDF file
  !   Do not access each member outside mod_netcdf.f90.
  type,public :: netcdf_handler

     !--- general
     integer             :: id_nc
     integer             :: count(4)
     integer             :: status = NOT_OPENED

     !--- global attribute
     character(CLEN)     :: title
     character(CLEN)     :: history
     character(CLEN)     :: comment

     !--- longitude
     integer             :: id_dim_lon
     integer             :: id_var_lon
     integer             :: imax = -1
     real(8),allocatable :: lon(:)
     character(CLEN)     :: lon_units

     !--- latitude
     integer             :: id_dim_lat
     integer             :: id_var_lat
     integer             :: jmax = -1
     real(8),allocatable :: lat(:)
     character(CLEN)     :: lat_units

     !--- level
     integer             :: id_dim_lev
     integer             :: id_var_lev
     integer             :: kmax = -1
     real(8),allocatable :: lev(:)
     character(CLEN)     :: lev_units

     !--- time
     integer             :: id_dim_time
     integer             :: id_var_time
     integer             :: tmax = -1
     real(8),allocatable :: time(:)
     character(CLEN)     :: time_units

     !--- variable (basic)
     integer             :: id_var
     character(CLEN)     :: var_name
     character(CLEN)     :: var_desc
     character(CLEN)     :: var_units
     real(4)             :: var_missing

     !--- variable (related to compression)
     integer             :: nf90_type
     integer             :: nf90_cmode
     real(4)             :: var_valid_min
     real(4)             :: var_valid_max
     real(4)             :: var_scale
     real(4)             :: var_offset
     integer(2)          :: var_missing_int2
     real(4)             :: var_force_to_set_valid_min    ! maximal value to force to set var_valid_min
     real(4)             :: var_force_to_set_valid_max    ! minimal value to force to set var_valid_max
     integer             :: chunksizes(4)
     logical             :: shuffle       = .true.   ! shuffle filter (affect compression rate)
     logical             :: fletcher32    = .false.  ! checksum
     integer             :: deflate_level = 1        ! level of compression

     !--- statistics monitor (for check)
     real(4)             :: var_org_min  ! actual minimal value of the original data
     real(4)             :: var_org_max  ! actual maximal value of the original data
     real(4)             :: var_out_min  ! actual minimal value of the output data
     real(4)             :: var_out_max  ! actual maximal value of the output data

  end type netcdf_handler

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tell unit number for log file (default is 6=STDOUT)
  subroutine netcdf_set_logfid( &
      fid )
    implicit none

    integer, intent(in) :: fid
    !---------------------------------------------------------------------------

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1

  end subroutine netcdf_set_logfid

  !-----------------------------------------------------------------------------
  subroutine netcdf_open_for_write( &
       nc,                         &
       ncfile,                     &
       count,                      &
       title,                      &
       history,                    &
       comment,                    &
       imax,                       &
       jmax,                       &
       kmax,                       &
       tmax,                       &
       lon,                        &
       lat,                        &
       lev,                        &
       time,                       &
       lon_units,                  &
       lat_units,                  &
       lev_units,                  &
       time_units,                 &
       var_name,                   &
       var_desc,                   &
       var_units,                  &
       var_missing,                &
       var_try_comp_2byte,         &
       var_comp_2byte,             &
       var_comp_hdf5,              &
       var_valid_min,              &
       var_valid_max,              &
       var_force_to_set_valid_min, &
       var_force_to_set_valid_max, &
       chunksizes,                 &
       var_comp_table_file         &  ! [add] C.Kodama 2014.05.09
       )
    type(netcdf_handler),intent(out)          :: nc        ! handler
    character(*),        intent(in)           :: ncfile    ! output netcdf filename
    integer,             intent(in) ,optional :: count(4)  ! unit array size to write
    character(*),        intent(in) ,optional :: title
    character(*),        intent(in) ,optional :: history
    character(*),        intent(in) ,optional :: comment

    integer,             intent(in) ,optional :: imax
    integer,             intent(in) ,optional :: jmax
    integer,             intent(in) ,optional :: kmax
    integer,             intent(in) ,optional :: tmax

    real(8),             intent(in) ,optional :: lon(:)
    real(8),             intent(in) ,optional :: lat(:)
    real(8),             intent(in) ,optional :: lev(:)
    real(8),             intent(in) ,optional :: time(:)

    character(*),        intent(in) ,optional :: lon_units
    character(*),        intent(in) ,optional :: lat_units
    character(*),        intent(in) ,optional :: lev_units
    character(*),        intent(in) ,optional :: time_units

    character(*),        intent(in) ,optional :: var_name  ! e.g. 'ms_tem', 'sa_u10m'
    character(*),        intent(in) ,optional :: var_desc  ! description of the variable
    character(*),        intent(in) ,optional :: var_units
    real(4),             intent(in) ,optional :: var_missing

    logical,             intent(in) ,optional :: var_try_comp_2byte ! just try to
    logical,             intent(in) ,optional :: var_comp_2byte     ! force to
    logical,             intent(in) ,optional :: var_comp_hdf5      ! force to
    real(4),             intent(in) ,optional :: var_valid_min
    real(4),             intent(in) ,optional :: var_valid_max
    real(4),             intent(in) ,optional :: var_force_to_set_valid_min
    real(4),             intent(in) ,optional :: var_force_to_set_valid_max
    integer,             intent(in) ,optional :: chunksizes(4)
    character(*),        intent(in) ,optional :: var_comp_table_file  ! [add] C.Kodama 2014.05.09

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_open_for_write


  subroutine netcdf_open_for_read( &
       nc,         &
       ncfile,     &
       count,      &
       title,      &
       history,    &
       comment,    &
       imax,       &
       jmax,       &
       kmax,       &
       tmax,       &
!       lon,        &
!       lat,        &
!       lev,        &
!       time,       &
       lon_units,  &
       lat_units,  &
       lev_units,  &
       time_units, &
       var_name,   &
       var_desc,   &
       var_units,  &
       var_missing &
       )
    type(netcdf_handler),intent(out)          :: nc
    character(*),        intent(in)           :: ncfile
    integer,             intent(out),optional :: count(4)
    character(*),        intent(out),optional :: title
    character(*),        intent(out),optional :: history
    character(*),        intent(out),optional :: comment

    integer,             intent(out),optional :: imax
    integer,             intent(out),optional :: jmax
    integer,             intent(out),optional :: kmax
    integer,             intent(out),optional :: tmax

!    real(8),allocatable, intent(out),optional :: lon(:)
!    real(8),allocatable, intent(out),optional :: lat(:)
!    real(8),allocatable, intent(out),optional :: lev(:)
!    real(8),allocatable, intent(out),optional :: time(:)

    character(*),        intent(out),optional :: lon_units
    character(*),        intent(out),optional :: lat_units
    character(*),        intent(out),optional :: lev_units
    character(*),        intent(out),optional :: time_units

    character(*),        intent(out),optional :: var_name
    character(*),        intent(out),optional :: var_desc
    character(*),        intent(out),optional :: var_units
    real(4),             intent(out),optional :: var_missing

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_open_for_read


  subroutine netcdf_set_count( nc, count )
    type(netcdf_handler),intent(inout) :: nc
    integer,             intent(in)    :: count(4)

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_set_count


  subroutine netcdf_write_0d( nc, var, i, j, k, t )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var
    integer,             intent(   in) :: i, j, k, t
    real(4),allocatable                :: wrk_var_real4(:,:,:,:)

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_write_0d

  subroutine netcdf_write_1d( nc, var, j, k, t )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var(:)
    integer,             intent(   in) :: j, k, t
    real(4),allocatable                :: wrk_var_real4(:,:,:,:)

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_write_1d

  subroutine netcdf_write_2d( nc, var, k, t )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var(:,:)
    integer,             intent(   in) :: k, t

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_write_2d

  subroutine netcdf_write_3d( nc, var, t )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var(:,:,:)
    integer,             intent(   in) :: t

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_write_3d

  subroutine netcdf_write_4d( nc, var )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var(:,:,:,:)

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_write_4d

  subroutine netcdf_write_main( nc, wrk_var_real4, i, j, k, t )
    type(netcdf_handler),intent(inout)          :: nc
    real(4),             intent(inout)          :: wrk_var_real4(:,:,:,:)
    integer,             intent(   in),optional :: i, j, k, t  ! record position to write

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_write_main


  subroutine netcdf_read_0d( nc, var, i, j, k, t )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var
    integer,             intent(in)  :: i, j, k, t

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_read_0d

  subroutine netcdf_read_1d( nc, var, j, k, t )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var(:)
    integer,             intent(in)  :: j, k, t

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_read_1d

  subroutine netcdf_read_2d( nc, var, k, t )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var(:,:)
    integer,             intent(in)  :: k, t

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_read_2d

  subroutine netcdf_read_3d( nc, var, t )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var(:,:,:)
    integer,             intent(in)  :: t

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_read_3d

  subroutine netcdf_read_4d( nc, var )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var(:,:,:,:)

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_read_4d

  subroutine netcdf_read_main( nc, wrk_var_real4, i, j, k, t )
    type(netcdf_handler),intent(in)           :: nc
    real(4),             intent(out)          :: wrk_var_real4(:,:,:,:)
    integer,             intent(in) ,optional :: i, j, k, t  ! record position to read

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_read_main


  subroutine netcdf_read_dim( &
       nc,         &
       lon,        &
       lat,        &
       lev,        &
       time        &
       )
    type(netcdf_handler),intent(in) :: nc
    real(8),intent(out),optional    :: lon(1:nc%imax)
    real(8),intent(out),optional    :: lat(1:nc%jmax)
    real(8),intent(out),optional    :: lev(1:nc%kmax)
    real(8),intent(out),optional    :: time(1:nc%tmax)

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_read_dim


  subroutine netcdf_close( nc )
    type(netcdf_handler),intent(inout) :: nc

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_close


  ! based on mod_grads.f90 and prg_ico2ll.f90 in NICAM
  subroutine netcdf_create_grads_ctl( nc, fid_ctl, endian )
    type(netcdf_handler),intent(in)           :: nc
    integer,             intent(in)           :: fid_ctl
    character(*),        intent(in) ,optional :: endian

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine netcdf_create_grads_ctl


  subroutine check( status )
    implicit none
    integer,intent(in)  :: status

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine check

  subroutine handle_error( status )
    implicit none
    integer, intent (in) :: status

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine handle_error


  subroutine var_comp_def_table( nc, var_name )
    type(netcdf_handler),intent(inout) :: nc
    character(*),        intent(   in) :: var_name

    write(*,*) 'netcdf is NOT supported in this build.'
    stop 1
  end subroutine var_comp_def_table

end module mod_netcdf
