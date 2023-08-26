!-------------------------------------------------------------------------------
!> Module file I/O PaNDA
!!
!! @par Description
!!         File I/O module (Packaged NICAM Data format)
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_fio_panda
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof

  use iso_c_binding
  use mod_fio_common
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FIO_PANDA_setup

  public :: FIO_PANDA_input_SP
  public :: FIO_PANDA_input_DP
  public :: FIO_PANDA_seek
  public :: FIO_PANDA_output_SP
  public :: FIO_PANDA_output_DP
  public :: FIO_PANDA_close
  public :: FIO_PANDA_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ C interface
  !
  include 'mod_fio_panda.inc'

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: FIO_PANDA_getfid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=FIO_HLONG), private :: FIO_PANDA_fname_list(FIO_file_nlim) = ''
  integer,                  private :: FIO_PANDA_fid_list  (FIO_file_nlim) = -1
  integer,                  private :: FIO_PANDA_fid_count = 1

  type(datainfo_panda), private :: dinfo

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup this module.
  !!
  !! Must be called first.
  !!
  subroutine FIO_PANDA_setup
    use mod_adm, only: &
       RGNMNG_lp2r, &
       ADM_prc_me,  &
       ADM_glevel,  &
       ADM_rlevel,  &
       ADM_lall
    implicit none

    integer, allocatable :: prc_tab(:)

    integer :: ierr
    !---------------------------------------------------------------------------

    allocate( prc_tab(ADM_lall) )

    prc_tab(1:ADM_lall) = RGNMNG_lp2r(1:ADM_lall,ADM_prc_me)-1

    ierr = fio_syscheck()
    ierr = fio_put_commoninfo( FIO_SPLIT_FILE,  &
                               FIO_BIG_ENDIAN,  &
                               FIO_ICOSAHEDRON, &
                               ADM_glevel,      &
                               ADM_rlevel,      &
                               ADM_lall,        &
                               prc_tab          )

    deallocate(prc_tab)

    return
  end subroutine FIO_PANDA_setup

  !-----------------------------------------------------------------------------
  !> Get file ID of given basename.
  !!
  !! Open it if not opened yet.
  !!
  subroutine FIO_PANDA_getfid( &
       fid,      &
       basename, &
       rwtype,   &
       pkg_desc, &
       pkg_note  )
    use mod_adm, only: &
       GLOBAL_prefix_dir,    &
       GLOBAL_extension_ens, &
       ADM_prc_me
    implicit none

    integer,          intent(out) :: fid      !< file ID
    character(len=*), intent(in)  :: basename !< basename of file
    integer,          intent(in)  :: rwtype   !< file access type
    character(len=*), intent(in)  :: pkg_desc !< package(file) description
    character(len=*), intent(in)  :: pkg_note !< package(file) note

    character(len=FIO_HSHORT) :: rwname(0:2)
    data rwname / 'READ','WRITE','APPEND' / ! [fix] H.Yashiro 20110912

    character(len=FIO_HLONG) :: basename_mod
    character(len=FIO_HLONG) :: fname
    integer                  :: n

    integer :: ierr
    !---------------------------------------------------------------------------

    basename_mod = trim(GLOBAL_prefix_dir)//trim(basename)//trim(GLOBAL_extension_ens)

    !--- search existing file
    fid = -1
    do n = 1, FIO_PANDA_fid_count
       if ( basename_mod == FIO_PANDA_fname_list(n) ) fid = FIO_PANDA_fid_list(n)
    enddo

    if ( fid < 0 ) then ! file registration
       !--- register new file and open
       call fio_mk_fname(fname,cstr(basename_mod),cstr('pe'),ADM_prc_me-1,6)
       fid = fio_register_file(fname)
       call fstr(fname)

       if ( rwtype == FIO_FREAD ) then

!          ierr = fio_dump_finfo(fid,FIO_BIG_ENDIAN,FIO_DUMP_HEADER) ! dump to stdout(check)
          ierr = fio_fopen(fid,FIO_FREAD)
          ierr = fio_read_allinfo(fid)

       elseif( rwtype == FIO_FWRITE ) then

          ierr = fio_fopen(fid,FIO_FWRITE)
          ierr = fio_put_write_pkginfo(fid,cstr(pkg_desc),cstr(pkg_note))

       elseif( rwtype == FIO_FAPPEND ) then

          ierr = fio_fopen(fid,FIO_FAPPEND)
          ierr = fio_read_pkginfo(fid)
          ierr = fio_write_pkginfo(fid)

       endif

       write(IO_FID_LOG,'(1x,A,A,A,I4)') '*** [FIO_PANDA] File registration (ADVANCED) : ', &
                                         trim(rwname(rwtype)),' - ', FIO_PANDA_fid_count
       write(IO_FID_LOG,'(1x,A,I4,A,A)') '***       fid= ', fid, ', name: ', trim(fname)

       FIO_PANDA_fname_list(FIO_PANDA_fid_count) = trim(basename_mod)
       FIO_PANDA_fid_list  (FIO_PANDA_fid_count) = fid
       FIO_PANDA_fid_count                       = FIO_PANDA_fid_count + 1
    endif

    return
  end subroutine FIO_PANDA_getfid

  !-----------------------------------------------------------------------------
  !> Input(read) one variable at one step.
  subroutine FIO_PANDA_input_SP( &
       var,            &
       basename,       &
       varname,        &
       layername,      &
       k_start,        &
       k_end,          &
       step,           &
       allow_missingq, &
       did_out         )
    use mod_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_gall, &
       ADM_lall
    implicit none

    real(SP),         intent(out) :: var(:,:,:)     !< variable(ij,k,l)
    character(len=*), intent(in)  :: basename       !< basename of file
    character(len=*), intent(in)  :: varname        !< variable name
    character(len=*), intent(in)  :: layername      !< layer name
    integer,          intent(in)  :: k_start        !< start index of vertical level
    integer,          intent(in)  :: k_end          !< end index of vertical level
    integer,          intent(in)  :: step           !< step to be read
    logical,          intent(in)  :: allow_missingq !< if data is missing, set value to zero, else execution stops.
    integer,          intent(out) :: did_out

    real(SP), allocatable, target :: var4(:,:,:)
    real(DP), allocatable, target :: var8(:,:,:)

    character(len=FIO_HSHORT) :: layername_file

    integer :: did, fid
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_in')

    !--- search/register file
    call FIO_PANDA_getfid( fid, basename, FIO_FREAD, '', '' )

    !--- seek data ID and get information
    did = fio_seek_datainfo(fid,cstr(varname),step)

    did_out = did

    !--- verify
    if ( did == -1 ) then
       if ( allow_missingq ) then
          write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] data not found! : ', &
                              'varname= ', trim(varname), ', step=', step
          write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] Q Value is set to 0.'

          var(:,k_start:k_end,:) = 0.0_SP

          call PROF_rapend  ('FILEIO_in')
          return
       else
          write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] data not found! : ', &
                              'varname= ', trim(varname), ', step=', step
          call PRC_MPIstop
       endif
    endif

    ierr = fio_get_datainfo(dinfo,fid,did)
    call fstr(layername_file,dinfo%layername)

!    write(IO_FID_LOG,*) "dinfo%varname     : ", dinfo%varname
!    write(IO_FID_LOG,*) "dinfo%description : ", dinfo%description
!    write(IO_FID_LOG,*) "dinfo%unit        : ", dinfo%unit
!    write(IO_FID_LOG,*) "dinfo%layername   : ", dinfo%layername
!    write(IO_FID_LOG,*) "dinfo%note        : ", dinfo%note
!    write(IO_FID_LOG,*) "dinfo%datasize    : ", dinfo%datasize
!    write(IO_FID_LOG,*) "dinfo%datatype    : ", dinfo%datatype
!    write(IO_FID_LOG,*) "dinfo%num_of_layer: ", dinfo%num_of_layer
!    write(IO_FID_LOG,*) "dinfo%step        : ", dinfo%step
!    write(IO_FID_LOG,*) "dinfo%time_start  : ", dinfo%time_start
!    write(IO_FID_LOG,*) "dinfo%time_end    : ", dinfo%time_end
!    call flush(IO_FID_LOG)

    if ( layername_file /= layername ) then
       write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] layername mismatch! ', &
                           '[', trim(layername_file), ':', trim(layername), '] ', trim(varname)
       call PRC_MPIstop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] num_of_layer mismatch! ', &
                           dinfo%num_of_layer, k_end-k_start+1, ' ', trim(varname)
       call PRC_MPIstop
    endif

    !--- read data
    if ( dinfo%datatype == FIO_REAL4 ) then

       allocate( var4(ADM_gall,dinfo%num_of_layer,ADM_lall) )

       ierr = fio_read_data(fid,did,c_loc(var4(:,:,:)))
       var(:,k_start:k_end,:) = real(var4(:,1:dinfo%num_of_layer,:),kind=SP)

       deallocate( var4 )

    elseif( dinfo%datatype == FIO_REAL8 ) then

       allocate( var8(ADM_gall,dinfo%num_of_layer,ADM_lall) )

       ierr = fio_read_data(fid,did,c_loc(var8(:,:,:)))
       var(:,k_start:k_end,:) = real(var8(:,1:dinfo%num_of_layer,:),kind=SP)

       deallocate( var8 )

    endif

    call PROF_rapend  ('FILEIO_in')

    return
  end subroutine FIO_PANDA_input_SP

  !-----------------------------------------------------------------------------
  !> Input(read) one variable at one step.
  subroutine FIO_PANDA_input_DP( &
       var,            &
       basename,       &
       varname,        &
       layername,      &
       k_start,        &
       k_end,          &
       step,           &
       allow_missingq, &
       did_out         )
    use mod_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_gall, &
       ADM_lall
    implicit none

    real(DP),         intent(out) :: var(:,:,:)     !< variable(ij,k,l)
    character(len=*), intent(in)  :: basename       !< basename of file
    character(len=*), intent(in)  :: varname        !< variable name
    character(len=*), intent(in)  :: layername      !< layer name
    integer,          intent(in)  :: k_start        !< start index of vertical level
    integer,          intent(in)  :: k_end          !< end index of vertical level
    integer,          intent(in)  :: step           !< step to be read
    logical,          intent(in)  :: allow_missingq !< if data is missing, set value to zero, else execution stops.
    integer,          intent(out) :: did_out

    real(SP), allocatable, target :: var4(:,:,:)
    real(DP), allocatable, target :: var8(:,:,:)

    character(len=FIO_HSHORT) :: layername_file

    integer :: did, fid
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_in')

    !--- search/register file
    call FIO_PANDA_getfid( fid, basename, FIO_FREAD, '', '' )

    !--- seek data ID and get information
    did = fio_seek_datainfo(fid,cstr(varname),step)

    !--- verify
    if ( did == -1 ) then
       if ( allow_missingq ) then
          write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] data not found! : ', &
                              'varname= ', trim(varname), ', step=', step
          write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] Q Value is set to 0.'

          var(:,k_start:k_end,:) = 0.0_DP

          call PROF_rapend  ('FILEIO_in')
          return
       else
          write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] data not found! : ', &
                              'varname= ', trim(varname), ', step=', step
          call PRC_MPIstop
       endif
    endif

    ierr = fio_get_datainfo(dinfo,fid,did)
    call fstr(layername_file,dinfo%layername)

!    write(IO_FID_LOG,*) "dinfo%varname     : ", dinfo%varname
!    write(IO_FID_LOG,*) "dinfo%description : ", dinfo%description
!    write(IO_FID_LOG,*) "dinfo%unit        : ", dinfo%unit
!    write(IO_FID_LOG,*) "dinfo%layername   : ", dinfo%layername
!    write(IO_FID_LOG,*) "dinfo%note        : ", dinfo%note
!    write(IO_FID_LOG,*) "dinfo%datasize    : ", dinfo%datasize
!    write(IO_FID_LOG,*) "dinfo%datatype    : ", dinfo%datatype
!    write(IO_FID_LOG,*) "dinfo%num_of_layer: ", dinfo%num_of_layer
!    write(IO_FID_LOG,*) "dinfo%step        : ", dinfo%step
!    write(IO_FID_LOG,*) "dinfo%time_start  : ", dinfo%time_start
!    write(IO_FID_LOG,*) "dinfo%time_end    : ", dinfo%time_end
!    call flush(IO_FID_LOG)

    did_out = did

    if ( layername_file /= layername ) then
       write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] layername mismatch! ', &
                           '[', trim(layername_file), ':', trim(layername), '] ', trim(varname)
       call PRC_MPIstop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] num_of_layer mismatch! ', &
                           dinfo%num_of_layer, k_end-k_start+1, ' ', trim(varname)
       call PRC_MPIstop
    endif

    !--- read data
    if ( dinfo%datatype == FIO_REAL4 ) then

       allocate( var4(ADM_gall,dinfo%num_of_layer,ADM_lall) )

       ierr = fio_read_data(fid,did,c_loc(var4(:,:,:)))
       var(:,k_start:k_end,:) = real(var4(:,1:dinfo%num_of_layer,:),kind=DP)

       deallocate( var4 )

    elseif( dinfo%datatype == FIO_REAL8 ) then

       allocate( var8(ADM_gall,dinfo%num_of_layer,ADM_lall) )

       ierr = fio_read_data(fid,did,c_loc(var8(:,:,:)))
       var(:,k_start:k_end,:) = real(var8(:,1:dinfo%num_of_layer,:),kind=DP)

       deallocate( var8 )

    endif

    call PROF_rapend  ('FILEIO_in')

    return
  end subroutine FIO_PANDA_input_DP

  !!-----------------------------------------------------------------------------
  !! Read in all steps of given `varname`, returns total data size(`num_of_step`)
  !! and mid of (time_start+time_end) of each step(`data_date`).
  !! `start_step` is maximum step where `ctime < 0.5*(ts(step)+te(step))` is true.
  !! `prec` is presicion, 4 or 8.
  !!
  !! If `opt_periodic_year` is T, data_date(:,1) is set as cdate(1) on return,
  !! else cdate(:) is neglected.
  !!
  subroutine FIO_PANDA_seek( &
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
    use mod_calendar, only: &
       CALENDAR_ss2yh, &
       CALENDAR_yh2ss
    implicit none

    integer,          intent(inout) :: start_step
    integer,          intent(inout) :: num_of_step
    integer,          intent(inout) :: data_date(6,FIO_data_nlim)
    integer,          intent(inout) :: prec
    character(len=*), intent(in)    :: basename
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: layername      ! for verification only
    integer,          intent(in)    :: k_start, k_end ! for verification only
    real(DP),         intent(in)    :: ctime
    integer,          intent(in)    :: cdate(6)       ! cdate(1) is only used only when opt_periodic_year is T.
    logical,          intent(in)    :: opt_periodic_year

    character(len=FIO_HSHORT) :: layername_file

    real(DP) :: midtime ! [sec]
    logical  :: startflag
    integer  :: did, fid
    integer  :: i
    integer  :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_in')

    !--- search/register file
    call FIO_PANDA_getfid( fid, basename, FIO_FREAD, '', '' )

    startflag      = .false.
    num_of_step    = -1
    data_date(:,:) = -1
    prec           = -1

    do i = 1, FIO_data_nlim
       !--- seek data ID and get information
       did = fio_seek_datainfo(fid,cstr(varname),i)

       if ( did == -1 ) then
          num_of_step = i - 1
          exit
       endif

       ierr = fio_get_datainfo(dinfo,fid,did)
       call fstr(layername_file,dinfo%layername)

!       write(IO_FID_LOG,*) "dinfo%varname     : ", dinfo%varname
!       write(IO_FID_LOG,*) "dinfo%description : ", dinfo%description
!       write(IO_FID_LOG,*) "dinfo%unit        : ", dinfo%unit
!       write(IO_FID_LOG,*) "dinfo%layername   : ", dinfo%layername
!       write(IO_FID_LOG,*) "dinfo%note        : ", dinfo%note
!       write(IO_FID_LOG,*) "dinfo%datasize    : ", dinfo%datasize
!       write(IO_FID_LOG,*) "dinfo%datatype    : ", dinfo%datatype
!       write(IO_FID_LOG,*) "dinfo%num_of_layer: ", dinfo%num_of_layer
!       write(IO_FID_LOG,*) "dinfo%step        : ", dinfo%step
!       write(IO_FID_LOG,*) "dinfo%time_start  : ", dinfo%time_start
!       write(IO_FID_LOG,*) "dinfo%time_end    : ", dinfo%time_end
!       call flush(IO_FID_LOG)

       !--- verify
       if ( layername_file /= layername ) then
          write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] layername mismatch! ', &
                              '[', trim(layername_file), ':', trim(layername), ']'
          call PRC_MPIstop
       elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
          write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] num_of_layer mismatch!', &
                              dinfo%num_of_layer, k_end-k_start+1
          call PRC_MPIstop
       endif

       midtime = real( int( (dinfo%time_start+dinfo%time_end)*0.5_DP+1.0_DP,kind=8 ),kind=DP ) ! specify int kind=8
       call CALENDAR_ss2yh( data_date(:,i), midtime )

       if ( opt_periodic_year ) then
          data_date(1,i) = cdate(1)
          call CALENDAR_yh2ss( midtime, data_date(:,i) )
       endif

       if (       ( .NOT. startflag ) &
            .AND. ( ctime < midtime ) ) then
          startflag  = .true.
          start_step = i
          prec       = FIO_preclist(dinfo%datatype)
       endif
    enddo

    if ( did /= -1 ) then
       write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO_PANDA] Date limit over! ', trim(varname), &
                           ' contains temporal data more than ', FIO_data_nlim
       call PRC_MPIstop
    endif

    call PROF_rapend  ('FILEIO_in')

    return
  end subroutine FIO_PANDA_seek

  !-----------------------------------------------------------------------------
  !! Append data with data header
  subroutine FIO_PANDA_output_SP( &
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
    use mod_adm, only: &
       ADM_gall, &
       ADM_lall
    use mod_const, only: &
       CONST_UNDEF4
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
    logical,          intent(in), optional  :: append

    real(SP), allocatable, target :: var4(:,:,:)
    real(DP), allocatable, target :: var8(:,:,:)

    integer :: file_mode
    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_out')

    !--- search/register file
    file_mode = FIO_FWRITE
    if ( present(append) ) then
       if ( append ) then
          file_mode = FIO_FAPPEND
       endif
    endif

    call FIO_PANDA_getfid( fid, basename, file_mode, pkg_desc, pkg_note )

    !--- append data to the file
    call cstr2(dinfo%varname    , varname  )
    call cstr2(dinfo%description, data_desc)
    call cstr2(dinfo%unit       , unit     )
    call cstr2(dinfo%layername  , layername)
    call cstr2(dinfo%note       , data_note)
    dinfo%datasize     = int( ADM_gall * ADM_lall * (k_end-k_start+1) * FIO_preclist(dtype), kind=8 )
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1
    dinfo%step         = step
    dinfo%time_start   = int( t_start, kind=8 )
    dinfo%time_end     = int( t_end,   kind=8 )

    !write(IO_FID_LOG,*) dinfo%varname
    !write(IO_FID_LOG,*) dinfo%description
    !write(IO_FID_LOG,*) dinfo%unit
    !write(IO_FID_LOG,*) dinfo%layername
    !write(IO_FID_LOG,*) dinfo%note
    !write(IO_FID_LOG,*) dinfo%datasize
    !write(IO_FID_LOG,*) dinfo%datatype
    !write(IO_FID_LOG,*) dinfo%num_of_layer
    !write(IO_FID_LOG,*) dinfo%step
    !write(IO_FID_LOG,*) dinfo%time_start
    !write(IO_FID_LOG,*) dinfo%time_end
    !call flush(IO_FID_LOG)

    if ( dtype == FIO_REAL4 ) then

       allocate( var4(ADM_gall,dinfo%num_of_layer,ADM_lall) )

       var4(:,1:dinfo%num_of_layer,:) = real(var(:,k_start:k_end,:),kind=SP)
       where( var4(:,:,:) < (CONST_UNDEF4+1.0_SP) )
          var4(:,:,:) = CONST_UNDEF4
       endwhere

       did = fio_put_write_datainfo_data(fid,dinfo,c_loc(var4(:,:,:)))

       deallocate( var4 )

    elseif( dtype == FIO_REAL8 ) then

       allocate( var8(ADM_gall,dinfo%num_of_layer,ADM_lall) )

       var8(:,1:dinfo%num_of_layer,:) = real(var(:,k_start:k_end,:),kind=DP)

       did = fio_put_write_datainfo_data(fid,dinfo,c_loc(var8(:,:,:)))

       deallocate( var8 )

    else
       write(IO_FID_LOG,*) 'xxx [OUTPUT]/[FIO_PANDA] Unsupported datatype!', dtype
       call PRC_MPIstop
    endif

    call PROF_rapend  ('FILEIO_out')

    return
  end subroutine FIO_PANDA_output_SP

  !-----------------------------------------------------------------------------
  !! Append data with data header
  subroutine FIO_PANDA_output_DP( &
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
    use mod_adm, only: &
       ADM_gall, &
       ADM_lall
    use mod_const, only: &
       CONST_UNDEF4
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
    logical,          intent(in), optional  :: append

    real(SP), allocatable, target :: var4(:,:,:)
    real(DP), allocatable, target :: var8(:,:,:)

    integer :: file_mode
    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_out')

    !--- search/register file
    file_mode = FIO_FWRITE
    if ( present(append) ) then
       if ( append ) then
          file_mode = FIO_FAPPEND
       endif
    endif

    call FIO_PANDA_getfid( fid, basename, file_mode, pkg_desc, pkg_note )

    !--- append data to the file
    call cstr2(dinfo%varname    , varname  )
    call cstr2(dinfo%description, data_desc)
    call cstr2(dinfo%unit       , unit     )
    call cstr2(dinfo%layername  , layername)
    call cstr2(dinfo%note       , data_note)
    dinfo%datasize     = int( ADM_gall * ADM_lall * (k_end-k_start+1) * FIO_preclist(dtype), kind=8 )
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1
    dinfo%step         = step
    dinfo%time_start   = int( t_start, kind=8 )
    dinfo%time_end     = int( t_end,   kind=8 )

    !write(IO_FID_LOG,*) dinfo%varname
    !write(IO_FID_LOG,*) dinfo%description
    !write(IO_FID_LOG,*) dinfo%unit
    !write(IO_FID_LOG,*) dinfo%layername
    !write(IO_FID_LOG,*) dinfo%note
    !write(IO_FID_LOG,*) dinfo%datasize
    !write(IO_FID_LOG,*) dinfo%datatype
    !write(IO_FID_LOG,*) dinfo%num_of_layer
    !write(IO_FID_LOG,*) dinfo%step
    !write(IO_FID_LOG,*) dinfo%time_start
    !write(IO_FID_LOG,*) dinfo%time_end
    !call flush(IO_FID_LOG)

    if ( dtype == FIO_REAL4 ) then

       allocate( var4(ADM_gall,dinfo%num_of_layer,ADM_lall) )

       var4(:,1:dinfo%num_of_layer,:) = real(var(:,k_start:k_end,:),kind=SP)
       where( var4(:,:,:) < (CONST_UNDEF4+1.0_SP) )
          var4(:,:,:) = CONST_UNDEF4
       endwhere

       did = fio_put_write_datainfo_data(fid,dinfo,c_loc(var4(:,:,:)))

       deallocate( var4 )

    elseif( dtype == FIO_REAL8 ) then

       allocate( var8(ADM_gall,dinfo%num_of_layer,ADM_lall) )

       var8(:,1:dinfo%num_of_layer,:) = real(var(:,k_start:k_end,:),kind=DP)

       did = fio_put_write_datainfo_data(fid,dinfo,c_loc(var8(:,:,:)))

       deallocate( var8 )

    else
       write(IO_FID_LOG,*) 'xxx [OUTPUT]/[FIO_PANDA] Unsupported datatype!', dtype
       call PRC_MPIstop
    endif

    call PROF_rapend  ('FILEIO_out')

    return
  end subroutine FIO_PANDA_output_DP

  !-----------------------------------------------------------------------------
  subroutine FIO_PANDA_close( &
       basename )
    use mod_adm, only: &
       GLOBAL_prefix_dir,    &
       GLOBAL_extension_ens, &
       ADM_prc_me
    implicit none

    character(len=*), intent(in) :: basename

    character(len=FIO_HLONG) :: basename_mod
    character(len=FIO_HLONG) :: fname

    integer :: fid, ierr
    integer :: n
    !---------------------------------------------------------------------------

    basename_mod = trim(GLOBAL_prefix_dir)//trim(basename)//trim(GLOBAL_extension_ens)

    !--- search/close file
    do n = 1, FIO_PANDA_fid_count
       if ( basename_mod == FIO_PANDA_fname_list(n) ) then
          fid = FIO_PANDA_fid_list(n)

          ierr = fio_fclose(fid)
          call fio_mk_fname(fname,cstr(FIO_PANDA_fname_list(n)),cstr('pe'),ADM_prc_me-1,6)
          call fstr(fname)

          write(IO_FID_LOG,'(1x,A,I4,A,A)') &
          '*** [FIO_PANDA] File close (ADVANCED) fid= ', fid, ', name: ', trim(fname)

          ! remove closed file info from the list
          FIO_PANDA_fname_list(n) = ''
          FIO_PANDA_fid_list  (n) = -1
       endif
    enddo

    return
  end subroutine FIO_PANDA_close

  !-----------------------------------------------------------------------------
  subroutine FIO_PANDA_finalize
    use mod_adm, only: &
       ADM_prc_me
    implicit none

    character(len=FIO_HLONG) :: fname
    integer                  :: n, fid, ierr
    !---------------------------------------------------------------------------

    do n = 1, FIO_PANDA_fid_count
       fid = FIO_PANDA_fid_list(n)

       ierr = fio_fclose(fid)
       call fio_mk_fname(fname,cstr(FIO_PANDA_fname_list(n)),cstr('pe'),ADM_prc_me-1,6)
       call fstr(fname)

       write(IO_FID_LOG,'(1x,A,I4,A,A)') &
       '*** [FIO_PANDA] File close (ADVANCED) fid= ', fid, ', name: ', trim(fname)
    enddo

    return
  end subroutine FIO_PANDA_finalize

end module mod_fio_panda
