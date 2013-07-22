!-------------------------------------------------------------------------------
!>
!! File I/O module
!!
!! @par Description
!!         This module is continer for file I/O
!!
!! @author H.Tomita, H.Yashiro
!!
!! @par History
!! @li      2011-07-27 (H.Tomita)  [NEW]
!! @li      2011-08-19 (H.Yashiro) Incorporate into NICAM
!! @li      2011-09-03 (H.Yashiro) Complete format specification
!! @li      2011-12-14 (T.Seiki)   allocatable => pointer in type structure
!! @li      2012-02-01 (T.Seiki)   fix array size over in array assignment
!!
!<
module mod_fio
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_debug
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
  public :: FIO_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !--- character length 
  integer, parameter, public :: FIO_HSHORT =  16
  integer, parameter, public :: FIO_HMID   =  64
  integer, parameter, public :: FIO_HLONG  = 256

  !--- data type 
  integer, parameter, public :: FIO_REAL4    = 0
  integer, parameter, public :: FIO_REAL8    = 1
  integer, parameter, public :: FIO_INTEGER4 = 2
  integer, parameter, public :: FIO_INTEGER8 = 3

  !--- data endian 
  integer, parameter, public :: FIO_UNKNOWN_ENDIAN = 0
  integer, parameter, public :: FIO_LITTLE_ENDIAN  = 1
  integer, parameter, public :: FIO_BIG_ENDIAN     = 2

  !--- topology 
  integer, parameter, public :: FIO_ICOSAHEDRON = 0
  integer, parameter, public :: FIO_IGA_LCP     = 1
  integer, parameter, public :: FIO_IGA_MLCP    = 2

  !--- file mode (partial or complete) 
  integer, parameter, public :: FIO_SPLIT_FILE = 0
  integer, parameter, public :: FIO_INTEG_FILE = 1

  !--- proccessor type 
  integer, parameter, public :: FIO_SINGLE_PROC = 0
  integer, parameter, public :: FIO_MULTI_PROC  = 1

  !--- action type 
  integer, parameter, public :: FIO_FREAD   = 0
  integer, parameter, public :: FIO_FWRITE  = 1
  integer, parameter, public :: FIO_FAPPEND = 2 ! [add] H.Yashiro 20110907 overwrite mode

  !--- data dump type 
  integer, parameter, public :: FIO_DUMP_OFF      = 0
  integer, parameter, public :: FIO_DUMP_HEADER   = 1
  integer, parameter, public :: FIO_DUMP_ALL      = 2
  integer, parameter, public :: FIO_DUMP_ALL_MORE = 3

  !--- struct for package infomation
  type, public :: headerinfo
     character(LEN=FIO_HLONG) :: fname
     character(LEN=FIO_HMID)  :: description
     character(LEN=FIO_HLONG) :: note
     integer                  :: num_of_data
     integer                  :: fmode
     integer                  :: endiantype
     integer                  :: grid_topology
     integer                  :: glevel
     integer                  :: rlevel
     integer                  :: num_of_rgn
     ! [Mod] 2011/12/14, T.Seiki
!!$  integer,allocatable      :: rgnid(:)
     integer, pointer         :: rgnid(:)
  endtype headerinfo

  !--- struct for data infomation
  type, public :: datainfo
     character(LEN=FIO_HSHORT) :: varname
     character(LEN=FIO_HMID)   :: description
     character(LEN=FIO_HSHORT) :: unit
     character(LEN=FIO_HSHORT) :: layername
     character(LEN=FIO_HLONG)  :: note
     integer(8)                :: datasize
     integer                   :: datatype
     integer                   :: num_of_layer
     integer                   :: step
     integer(8)                :: time_start
     integer(8)                :: time_end
  endtype datainfo

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,             parameter, private :: FIO_nmaxfile = 64
  character(LEN=FIO_HLONG), save, private :: FIO_fname_list(FIO_nmaxfile)
  integer,                  save, private :: FIO_fid_list  (FIO_nmaxfile)
  integer,                  save, private :: FIO_fid_count = 1

  type(headerinfo), private :: hinfo 
  type(datainfo),   private :: dinfo 

  integer, private, parameter :: max_num_of_data = 2500 !--- max time step num
  integer, parameter, private :: preclist(0:3) = (/ 4, 8, 4, 8 /)

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  subroutine FIO_setup
    use mod_adm, only : &
      ADM_LOG_FID,   &
      ADM_prc_me,    &
      ADM_prc_tab,   &
      ADM_glevel,    &
      ADM_rlevel,    &
      ADM_lall
    implicit none

    integer, allocatable :: prc_tab(:) 
    !---------------------------------------------------------------------------

    ! dummy call
    call DEBUG_rapstart('FILEIO in')
    call DEBUG_rapend  ('FILEIO in')
    call DEBUG_rapstart('FILEIO out')
    call DEBUG_rapend  ('FILEIO out')

    allocate( prc_tab(ADM_lall) )
    ! [fix] 20120201 T.Seiki
    !!$ prc_tab(:) = ADM_prc_tab(:,ADM_prc_me)-1
    prc_tab(1:ADM_lall) = ADM_prc_tab(1:ADM_lall,ADM_prc_me)-1
    
    call fio_syscheck()
    call fio_put_commoninfo( FIO_SPLIT_FILE,  &
                             FIO_BIG_ENDIAN,  &
                             FIO_ICOSAHEDRON, &
                             ADM_glevel,      &
                             ADM_rlevel,      &
                             ADM_lall,        &
                             prc_tab          )

    deallocate(prc_tab)

    allocate( hinfo%rgnid(ADM_lall) )

    return
  end subroutine FIO_setup

  !-----------------------------------------------------------------------------
  subroutine FIO_getfid( &
      fid,      &
      basename, &
      rwtype,   &
      pkg_desc, &
      pkg_note  )
    use mod_adm, only : &
      ADM_LOG_FID, &
      ADM_prc_me
    implicit none

    integer,          intent(out) :: fid
    character(LEN=*), intent( in) :: basename
    integer,          intent( in) :: rwtype
    character(LEN=*), intent( in) :: pkg_desc
    character(LEN=*), intent( in) :: pkg_note

    character(LEN=FIO_HSHORT) :: rwname(0:2)
    data rwname / 'READ','WRITE','APPEND' / ! [fix] H.Yashiro 20110912

    character(LEN=FIO_HLONG) :: fname
    integer                  :: n
    !---------------------------------------------------------------------------

    !--- search existing file
    fid = -1
    do n = 1, FIO_fid_count
       if ( trim(basename)==trim(FIO_fname_list(n)) ) fid = FIO_fid_list(n)
    enddo

    if ( fid < 0 ) then ! file registration
       !--- register new file and open
       call fio_mk_fname(fname,trim(basename),'pe',ADM_prc_me-1,6)
       call fio_register_file(n,fname)

       if ( rwtype == FIO_FREAD ) then

!          call fio_dump_finfo(n,FIO_BIG_ENDIAN,FIO_DUMP_HEADER) ! dump to stdout(check)
          call fio_fopen(n,FIO_FREAD)
          call fio_read_allinfo(n)

       elseif( rwtype == FIO_FWRITE ) then

          call fio_fopen(n,FIO_FWRITE)
          call fio_put_write_pkginfo(n,pkg_desc,pkg_note)

       endif

       write(ADM_LOG_FID,*) '*** [FIO] File registration : ',trim(rwname(rwtype)),'-', n
       write(ADM_LOG_FID,*) '*** filename: ', trim(fname)

       FIO_fname_list(FIO_fid_count) = trim(basename)
       FIO_fid_list  (FIO_fid_count) = n
       FIO_fid_count = FIO_fid_count + 1
       fid = n
    endif

    return
  end subroutine FIO_getfid

  !-----------------------------------------------------------------------------
  subroutine FIO_input( &
      var,           &
      basename,      &
      varname,       &
      layername,     &
      k_start,       &
      k_end,         &
      step,          &
      allow_missingq ) !--- optional
    use mod_adm, only : &
      ADM_proc_stop, &
      ADM_LOG_FID, &
      ADM_gall,    &
      ADM_lall
    implicit none

    real(8),          intent(out) :: var(:,:,:)
    character(LEN=*), intent( in) :: basename
    character(LEN=*), intent( in) :: varname
    character(LEN=*), intent( in) :: layername
    integer,          intent( in) :: k_start, k_end
    integer,          intent( in) :: step

    logical, intent(in), optional :: allow_missingq !--- if data is missing, set value to zero

    real(4) :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(8) :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('FILEIO in')

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FREAD, "", "" )

    !--- seek data ID and get information
    call fio_seek_datainfo(did,fid,varname,step)
    call fio_get_datainfo(fid,did,dinfo)

!    write(ADM_LOG_FID,*) dinfo%varname
!    write(ADM_LOG_FID,*) dinfo%description
!    write(ADM_LOG_FID,*) dinfo%unit
!    write(ADM_LOG_FID,*) dinfo%layername
!    write(ADM_LOG_FID,*) dinfo%note
!    write(ADM_LOG_FID,*) dinfo%datasize
!    write(ADM_LOG_FID,*) dinfo%datatype
!    write(ADM_LOG_FID,*) dinfo%num_of_layer
!    write(ADM_LOG_FID,*) dinfo%step
!    write(ADM_LOG_FID,*) dinfo%time_start
!    write(ADM_LOG_FID,*) dinfo%time_end

    !--- verify
    if ( did == -1 ) then
       if ( present(allow_missingq) ) then ! [bugfix] H.Yashiro 20110912
          if ( allow_missingq ) then
             write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] data not found! : ', &
                                  'varname= ',trim(varname),', step=',step
             write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] Q Value is set to 0.'

             var(:,k_start:k_end,:) = 0.D0

             return
          endif
       else
          write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] data not found! : ', &
                               'varname= ',trim(varname),', step=',step
          call ADM_proc_stop
       endif
    endif

    if ( trim(dinfo%layername) /= trim(layername) ) then
       write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                            "[",trim(dinfo%layername),":",trim(layername),"]"
       call ADM_proc_stop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch! ', &
                            dinfo%num_of_layer,k_end-k_start+1
       call ADM_proc_stop
    endif

    !--- read data
    if ( dinfo%datatype == FIO_REAL4 ) then

       call fio_read_data(fid,did,var4(:,:,:))
       var(:,k_start:k_end,:) = real(var4(:,1:dinfo%num_of_layer,:),kind=8)

    elseif( dinfo%datatype == FIO_REAL8 ) then

       call fio_read_data(fid,did,var8(:,:,:))
       var(:,k_start:k_end,:) = var8(:,1:dinfo%num_of_layer,:)

    endif

    call DEBUG_rapend('FILEIO in')

    return
  end subroutine FIO_input

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
    use mod_adm, only : &
      ADM_proc_stop, &
      ADM_LOG_FID
    use mod_calendar, only :&
      calendar_ss2yh, &
      calendar_yh2ss
    implicit none

    integer,          intent(inout) :: start_step
    integer,          intent(inout) :: num_of_step
    integer,          intent(inout) :: data_date(6,max_num_of_data)
    integer,          intent(inout) :: prec
    character(len=*), intent(in) :: basename
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: layername
    integer,          intent(in) :: k_start, k_end
    real(8),          intent(in) :: ctime
    integer,          intent(in) :: cdate(6)
    logical,          intent(in) :: opt_periodic_year

    real(8) :: midtime !--- [sec]
    logical :: startflag
    integer :: did, fid
    integer :: i
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('FILEIO in')

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FREAD, "", "" )

    startflag = .false.

    do i = 1, max_num_of_data
       !--- seek data ID and get information
       call fio_seek_datainfo(did,fid,varname,i)
       call fio_get_datainfo(fid,did,dinfo)

!       write(ADM_LOG_FID,*) dinfo%varname
!       write(ADM_LOG_FID,*) dinfo%description
!       write(ADM_LOG_FID,*) dinfo%unit
!       write(ADM_LOG_FID,*) dinfo%layername
!       write(ADM_LOG_FID,*) dinfo%note
!       write(ADM_LOG_FID,*) dinfo%datasize
!       write(ADM_LOG_FID,*) dinfo%datatype
!       write(ADM_LOG_FID,*) dinfo%num_of_layer
!       write(ADM_LOG_FID,*) dinfo%step
!       write(ADM_LOG_FID,*) dinfo%time_start
!       write(ADM_LOG_FID,*) dinfo%time_end

       if ( did == -1 ) then
          num_of_step = i - 1
          exit
       endif

       !--- verify
       if ( trim(dinfo%layername) /= trim(layername) ) then
          write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                               "[",trim(dinfo%layername),":",trim(layername),"]"
          call ADM_proc_stop
       elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
          write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch!', &
                               dinfo%num_of_layer,k_end-k_start+1
          call ADM_proc_stop
       endif

       ! [fix] H.Yashiro 20111011 : specify int kind=8
       midtime = dble( int( (dinfo%time_start+dinfo%time_end)*0.5D0+1.D0,kind=8 ) )
       call calendar_ss2yh( data_date(:,i), midtime )

       if ( opt_periodic_year ) then
          data_date(1,i) = cdate(1)
          call calendar_yh2ss( midtime, data_date(:,i) )
       endif

       if (       ( .not. startflag ) &
            .AND. ( ctime < midtime ) ) then
          startflag  = .true.
          start_step = i
          prec       = preclist(dinfo%datatype)
       endif
    enddo

    call DEBUG_rapend('FILEIO in')

    return
  end subroutine FIO_seek

  !-----------------------------------------------------------------------------
  subroutine FIO_output( &
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
      t_end      )
    use mod_adm, only : &
      ADM_proc_stop, &
      ADM_LOG_FID, &
      ADM_prc_me, &
      ADM_gall,   &
      ADM_lall
    use mod_cnst, only : &
      CNST_UNDEF4
    implicit none

    real(8),          intent(in) :: var(:,:,:)
    character(LEN=*), intent(in) :: basename
    character(LEN=*), intent(in) :: pkg_desc
    character(LEN=*), intent(in) :: pkg_note
    character(LEN=*), intent(in) :: varname
    character(LEN=*), intent(in) :: data_desc
    character(LEN=*), intent(in) :: data_note
    character(LEN=*), intent(in) :: unit
    integer,          intent(in) :: dtype
    character(LEN=*), intent(in) :: layername
    integer,          intent(in) :: k_start, k_end
    integer,          intent(in) :: step
    real(8),          intent(in) :: t_start, t_end

    real(4) :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(8) :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('FILEIO out')

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FWRITE, pkg_desc, pkg_note )

    !--- append data to the file
    dinfo%varname      = varname
    dinfo%description  = data_desc
    dinfo%unit         = unit
    dinfo%layername    = layername
    dinfo%note         = data_note
    dinfo%datasize     = int( ADM_gall * ADM_lall * (k_end-k_start+1) * preclist(dtype), kind=8 )
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1
    dinfo%step         = step
    dinfo%time_start   = int( t_start, kind=8 )
    dinfo%time_end     = int( t_end,   kind=8 )

!    write(ADM_LOG_FID,*) dinfo%varname
!    write(ADM_LOG_FID,*) dinfo%description
!    write(ADM_LOG_FID,*) dinfo%unit
!    write(ADM_LOG_FID,*) dinfo%layername
!    write(ADM_LOG_FID,*) dinfo%note
!    write(ADM_LOG_FID,*) dinfo%datasize
!    write(ADM_LOG_FID,*) dinfo%datatype
!    write(ADM_LOG_FID,*) dinfo%num_of_layer
!    write(ADM_LOG_FID,*) dinfo%step
!    write(ADM_LOG_FID,*) dinfo%time_start
!    write(ADM_LOG_FID,*) dinfo%time_end

    if ( dtype == FIO_REAL4 ) then

       var4(:,k_start:k_end,:)=real(var(:,k_start:k_end,:),kind=4)
       where( var4(:,:,:) < (CNST_UNDEF4+1.0) )
          var4(:,:,:) = CNST_UNDEF4
       endwhere

       call fio_put_write_datainfo_data(did,fid,dinfo,var4(:,:,:))

    elseif( dtype == FIO_REAL8 ) then

       var8(:,k_start:k_end,:)=var(:,k_start:k_end,:)

       call fio_put_write_datainfo_data(did,fid,dinfo,var8(:,:,:))
    else
       write(ADM_LOG_FID,*) 'xxx [OUTPUT]/[FIO] Unsupported datatype!', dtype
       call ADM_proc_stop
    endif

    call DEBUG_rapend('FILEIO out')

    return
  end subroutine FIO_output

  !-----------------------------------------------------------------------------
  subroutine FIO_finalize
    use mod_adm, only : &
      ADM_LOG_FID, &
      ADM_prc_me
    implicit none

    character(LEN=FIO_HLONG) :: fname
    integer                  :: n
    !---------------------------------------------------------------------------

    do n = 1, FIO_fid_count
       call fio_fclose(FIO_fid_list(n))

       write(ADM_LOG_FID,*) '*** [FIO] File Close : NO.', n
       call fio_mk_fname(fname,trim(FIO_fname_list(n)),'pe',ADM_prc_me-1,6)
       write(ADM_LOG_FID,*) '*** closed filename: ', trim(fname)
    enddo

    return
  end subroutine FIO_finalize

end module mod_fio
!-------------------------------------------------------------------------------
