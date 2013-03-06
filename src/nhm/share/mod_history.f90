!-------------------------------------------------------------------------------
!>
!! History module
!!
!! @par Description
!!         This module is for managing the output variables
!!
!! @author M.Satoh
!!
!! @par History
!! @li      2005-11-29 (M.Satoh)    [new]
!! @li      2005-12-14 (S.Iga)      ADM_gall
!! @li      2006-02-16 (M.Satoh)    T. Mitsui: correct timing
!! @li      2006-08-07 (W.Yanase)   v_save=0 in history_setup, NO_VINTRPL
!! @li      2007-01-19 (K.Suzuki)   higher vectorized rate and allowing undefined value in average
!! @li      2007-06-27 (Y.Niwa)     add MONTHLY_AVERAGE option add ktype 'GL' 'GO' options
!! @li      2007-07-02 (Y.Niwa)     bug fix
!! @li      2007-11-30 (Y.Niwa)     add option for output at pressure levels
!! @li      2007-12-05 (T.Mitsui)   bug fix
!! @li      2008-05-30 (T.Mitsui)   distinguish w-grid, and option of v_interpolation
!! @li      2009-07-13 (S.Iga)      check_count is added. (nmhist miswriting checker)
!! @li      2010-05-11 (M.Satoh)    add l_region in history_in 
!! @li      2011-04-26 (C.Kodama)   support >10000 time steps
!! @li      2011-09-03 (H.Yashiro)  New I/O
!! @li      2012-01-26 (Y.Yamada)   trivial bug fix
!! @li      2012-03-28 (T.Seiki)    fix undefined reference
!! @li      2012-06-07 (T.Seiki)    add output_path for multi-job run
!! @li      2012-11-05 (H.Yashiro)  NICAM milestone project (Phase I:cleanup of shared module)
!!
!<
module mod_history
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_MAXFNAME, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: history_setup
  public :: history_in
  public :: history_out

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,                 public,              save :: HIST_req_nmax
  character(len=ADM_NSYS), public, allocatable, save :: item_save(:)
  logical,                 public,              save :: HIST_output_step0 = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: history_outlist
  private :: history_timeinfo
  private :: get_log_pres

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: HIST_req_limit = 1000

  character(len=ADM_MAXFNAME), private, save :: HIST_io_fname  = ''
  character(len=ADM_NSYS),     private, save :: HIST_io_desc   = ''
  integer,                     private, save :: HIST_dtype     = -1
  character(len=ADM_MAXFNAME), private, save :: output_path    = ''
  character(len=ADM_NSYS),     private, save :: histall_fname  = ''
  character(len=ADM_MAXFNAME), private, save :: output_io_mode != 'LEGACY'
  logical,                     private, save :: direct_access  = .false.
  integer,                     private, save :: output_size    = 4
  integer,                     private, save :: npreslev       = 1
  real(8),                     private, save :: pres_levs(60)  != CNST_PRE00
  logical,                     private, save :: check_flag     = .true.

  integer,                     private, save :: ksum
  logical,                     private, save :: calc_pressure = .false.

  character(len=ADM_MAXFNAME), private, allocatable, save :: file_save (:)
  character(len=ADM_NSYS),     private, allocatable, save :: desc_save (:)
  character(len=ADM_NSYS),     private, allocatable, save :: unit_save (:)
  integer,                     private, allocatable, save :: step_save (:)
  character(len=ADM_NSYS),     private, allocatable, save :: ktype_save(:)
  integer,                     private, allocatable, save :: kstr_save (:)
  integer,                     private, allocatable, save :: kend_save (:)
  integer,                     private, allocatable, save :: kmax_save (:)
  character(len=ADM_NSYS),     private, allocatable, save :: output_type_save  (:)
  logical,                     private, allocatable, save :: out_prelev_save   (:)
  logical,                     private, allocatable, save :: out_vintrpl_save  (:)
  logical,                     private, allocatable, save :: opt_wgrid_save    (:)
  logical,                     private, allocatable, save :: opt_lagintrpl_save(:)

  character(len=ADM_NSYS),     private, allocatable, save :: lname_save   (:)
  integer,                     private, allocatable, save :: tmax_save    (:)
  real(8),                     private, allocatable, save :: tstr_save    (:)
  real(8),                     private, allocatable, save :: tend_save    (:)
  integer,                     private, allocatable, save :: month_old    (:)
  integer,                     private, allocatable, save :: l_region_save(:)

  integer,                     private, allocatable, save :: ksumstr  (:)
  integer,                     private, allocatable, save :: ksumend  (:)
  real(8),                     private, allocatable, save :: tsum_save(:,:)
  logical,                     private, allocatable, save :: flag_save(:)

  real(8),                     private, allocatable, save :: v_save   (:,:,:,:)
  real(8),                     private, allocatable, save :: v_save_pl(:,:,:,:)
  real(8),                     private, allocatable, save :: zlev_save(:)

  real(8),                     private, allocatable, save :: pres_levs_ln(:)
  integer,                     private, allocatable, save :: cnvpre_klev(:,:,:)
  real(8),                     private, allocatable, save :: cnvpre_fac1(:,:,:)
  real(8),                     private, allocatable, save :: cnvpre_fac2(:,:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine history_setup
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_kmin,      &
       ADM_kmax,      &
       ADM_vlayer
    use mod_cnst, only: &
       CNST_PRE00
    use mod_fio, only: &
       FIO_REAL8, &
       FIO_REAL4
    use mod_time, only: &
       TIME_CTIME
    use mod_calendar, only: &
       calendar_ss2yh
    use mod_grd, only: &
       GRD_gz
    use mod_runconf, only: &
       RUN_TYPE!,          &
!         NTAU_ISCCP,        &
!         NPRES_ISCCP,       &
!         LAND_TYPE,         &
!         OCEAN_TYPE
!    use mod_landvar_bc, only: &
!         KMAX_bc    => KMAX,    &
!         GLVNAME_bc => GLVNAME, &
!         I_ALL_bc   => I_ALL
!    use mod_landvar_matsiro, only: &
!         KMAX_mat    => KMAX,    &
!         GLVNAME_mat => GLVNAME, &
!         I_ALL_mat   => I_ALL
!    use mod_oceanvar_mixedlayer, only: &
!         KMAX_ocn  => KMAX, &
!         GOVNAME,           &
!         I_ALL_ocn => I_ALL
    implicit none

    character(len=ADM_NSYS)     :: hist3D_layername  != ''
    integer                     :: step_def          = 1
    character(len=ADM_NSYS)     :: ktype_def         != ''
    integer                     :: kstr_def          = 1
    integer                     :: kend_def          != ADM_vlayer
    integer                     :: kmax_def          != ADM_vlayer
    character(len=ADM_NSYS)     :: output_type_def   != 'SNAPSHOT'
    logical                     :: out_prelev_def    = .false.
    logical                     :: no_vintrpl        = .true.
    logical                     :: opt_wgrid_def     = .false.
    logical                     :: opt_lagintrpl_def = .true.
    logical                     :: doout_step0

    character(len=ADM_NSYS)     :: item
    character(len=ADM_MAXFNAME) :: file
    character(len=ADM_NSYS)     :: desc
    character(len=ADM_NSYS)     :: unit
    integer                     :: step
    character(len=ADM_NSYS)     :: ktype
    integer                     :: kstr
    integer                     :: kend
    integer                     :: kmax
    character(len=ADM_NSYS)     :: output_type
    logical                     :: out_prelev
    logical                     :: out_vintrpl
    logical                     :: opt_wgrid     
    logical                     :: opt_lagintrpl 

    namelist / NMHISD / &
         output_path,       &
         histall_fname,     &
         hist3D_layername,  &
         output_io_mode,    &
         direct_access,     &
         output_size,       &
         step,              &
         ktype,             &
         kstr,              &
         kend,              &
         kmax,              &
         output_type,       &
         out_prelev,        &
         no_vintrpl,        &
         opt_wgrid_def,     &
         opt_lagintrpl_def, &
         npreslev,          &
         pres_levs,         &
         check_flag,        &
         doout_step0

    namelist / NMHIST / &
         item,         &
         file,         &
         desc,         &
         unit,         &
         step,         &
         ktype,        &
         kstr,         &
         kend,         &
         kmax,         &
         output_type,  &
         out_prelev,   &
         out_vintrpl,  &
         opt_wgrid,    &
         opt_lagintrpl

    character(len=ADM_NSYS) :: lname

    integer :: idate(6)
    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    ! set default
    output_path       = ''
    histall_fname     = ''
    hist3D_layername  = ''
    output_io_mode    = 'LEGACY'
    ktype_def         = 'unknown'
    kend_def          = ADM_vlayer
    kmax_def          = ADM_vlayer
    output_type_def   = 'SNAPSHOT'
    pres_levs(:)      = CNST_PRE00

    ! nonsence prepare
    step        = step_def
    ktype       = ktype_def
    kstr        = kstr_def
    kend        = kend_def
    kmax        = kmax_def
    output_type = output_type_def
    out_prelev  = out_prelev_def

    doout_step0 = HIST_output_step0

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[history]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=NMHISD,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** NMHISD is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist NMHISD. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist NMHISD. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,NMHISD)

    ! nonsence restore
    step_def        = step
    ktype_def       = ktype
    kstr_def        = kstr
    kend_def        = kend
    kmax_def        = kmax
    output_type_def = output_type
    out_prelev_def  = out_prelev

    HIST_output_step0 = doout_step0

    if (      trim(output_io_mode) == 'ADVANCED' &
         .OR. trim(output_io_mode) == 'LEGACY'   ) then
       write(ADM_LOG_FID,*) '*** History output type:', trim(output_io_mode)
    else
       write(ADM_LOG_FID,*) 'xxx Invalid output_io_mode!', trim(output_io_mode)
       call ADM_proc_stop
    endif
    HIST_io_fname = trim(output_path)//trim(histall_fname)
    HIST_io_desc  = trim(RUN_TYPE)

    if ( output_size == 4 ) then
       HIST_dtype = FIO_REAL4
    elseif ( output_size == 8 ) then
       HIST_dtype = FIO_REAL8
    else
       write(*,*) 'output_size is not appropriate:',output_size
       call ADM_proc_stop
    endif


    ! listup history request
    rewind(ADM_CTL_FID)
    do n = 1, HIST_req_limit
       read(ADM_CTL_FID,nml=NMHIST,iostat=ierr)
       if ( ierr < 0 ) then
          exit
       elseif( ierr > 0 ) then
          write(*,          *) 'xxx Not appropriate names in namelist NMHIST. STOP.'
          write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist NMHIST. STOP.'
          call ADM_proc_stop
      endif
    enddo
    HIST_req_nmax = n - 1

    if    ( HIST_req_nmax > HIST_req_limit ) then
       write(ADM_LOG_FID,*) '*** request of history file is exceed! n >', HIST_req_limit
    elseif( HIST_req_nmax == 0 ) then
       write(ADM_LOG_FID,*) '*** No history file specified.'
       return
    else
       write(ADM_LOG_FID,*) '*** Number of requested history item : ', HIST_req_nmax
    endif

    allocate( item_save         (HIST_req_nmax) )
    allocate( file_save         (HIST_req_nmax) )
    allocate( desc_save         (HIST_req_nmax) )
    allocate( unit_save         (HIST_req_nmax) )
    allocate( step_save         (HIST_req_nmax) )
    allocate( ktype_save        (HIST_req_nmax) )
    allocate( kstr_save         (HIST_req_nmax) )
    allocate( kend_save         (HIST_req_nmax) )
    allocate( kmax_save         (HIST_req_nmax) )
    allocate( output_type_save  (HIST_req_nmax) )
    allocate( out_prelev_save   (HIST_req_nmax) )
    allocate( out_vintrpl_save  (HIST_req_nmax) )
    allocate( opt_wgrid_save    (HIST_req_nmax) )
    allocate( opt_lagintrpl_save(HIST_req_nmax) )
    item_save         (:) = ""
    file_save         (:) = ""
    desc_save         (:) = ""
    unit_save         (:) = ""
    step_save         (:) = 0
    ktype_save        (:) = ""
    kstr_save         (:) = -1
    kend_save         (:) = -1
    kmax_save         (:) = 0
    output_type_save  (:) = ""
    out_prelev_save   (:) = .false.
    out_vintrpl_save  (:) = .false.
    opt_wgrid_save    (:) = .false.
    opt_lagintrpl_save(:) = .false.

    allocate( lname_save        (HIST_req_nmax) )
    allocate( tmax_save         (HIST_req_nmax) )
    allocate( tstr_save         (HIST_req_nmax) )
    allocate( tend_save         (HIST_req_nmax) )
    allocate( month_old         (HIST_req_nmax) )
    allocate( l_region_save     (HIST_req_nmax) )
    allocate( flag_save         (HIST_req_nmax) )
    lname_save        (:) = ""
    tmax_save         (:) = 0
    tstr_save         (:) = 0.D0
    tend_save         (:) = 0.D0
    month_old         (:) = 0
    l_region_save     (:) = 0
    flag_save         (:) = .false.

    call calendar_ss2yh( idate, TIME_CTIME )

    rewind(ADM_CTL_FID)
    do n = 1, HIST_req_limit

       ! set default
       item          = ''
       file          = ''
       desc          = ''
       unit          = 'NIL'
       step          = step_def
       ktype         = ktype_def
       kstr          = kstr_def
       kend          = kend_def
       kmax          = -1
       output_type   = output_type_def
       out_prelev    = out_prelev_def
       if ( no_vintrpl ) then
          out_vintrpl  = .false.
       else
          out_vintrpl  = .true.
       endif
       opt_wgrid     = opt_wgrid_def
       opt_lagintrpl = opt_lagintrpl_def

       ! read namelist
       read(ADM_CTL_FID,nml=NMHIST,iostat=ierr)
       if( ierr /= 0 ) exit

       if ( item == '' ) then
          write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist NMHIST. STOP.'
          call ADM_proc_stop
       endif

       if( file == '' ) file = item

       ! set default layername
       if ( kmax == 1 ) then
          lname = "ZSSFC1"
       else
          lname = "LAYERNM"
       endif

       select case( trim(ktype) )
       case('3D')
          if ( out_prelev ) then
             kstr = 1
             kend = npreslev
          else
             kstr = ADM_kmin
             kend = ADM_kmax
          endif
          lname = hist3D_layername
       case('2D')
          kstr = 1
          kend = 1
          lname = "ZSSFC1"
!       case('ISCCP')
!          kstr = 1
!          kend = NTAU_ISCCP*NPRES_ISCCP
!       case('GL')
!          kstr =  1
!          kend = -1
!          if ( trim(LAND_TYPE) == 'BUCKET' ) then
!             do idx = 1, I_ALL_bc
!                if( trim(GLVNAME_bc(idx) ) == trim(item) ) kend = KMAX_bc(idx)
!             enddo
!          elseif( trim(LAND_TYPE) == 'MATSIRO' ) then
!             do idx = 1, I_ALL_mat
!                if( trim(GLVNAME_mat(idx)) == trim(item) ) kend = KMAX_mat(idx)
!             enddo
!          endif
!          if ( kend == -1 ) then
!             write(ADM_LOG_FID,*) 'xxx History item=', trim(item), ' is not in LAND module=', trim(LAND_TYPE),'. STOP.'
!             call ADM_proc_stop
!          endif
!       case('GO')
!          kstr =  1
!          kend = -1
!          if ( trim(OCEAN_TYPE) == 'MIXEDLAYER' ) then
!              do idx = 1, I_ALL_ocn
!                 if( trim(GOVNAME(idx)) == trim(item) ) kend = KMAX_ocn(idx)
!              enddo
!          endif
!          if ( kend == -1 ) then
!             write(ADM_LOG_FID,*) 'xxx History item=', trim(item), ' is not in OCEAN module=', trim(OCEAN_TYPE),'. STOP.'
!             call ADM_proc_stop
!          endif
       endselect

       ! check consistensy between kend and kmax
       if ( kmax > 0 ) then
          if ( kmax /= kend - kstr + 1 ) then
             kend = kstr + kmax - 1
          endif
       else
          kmax = kend - kstr + 1
       endif

       if ( out_prelev ) then
          if ( ktype /= '3D' ) then
             write(ADM_LOG_FID,*) '*** Only 3D vars can be output by pressure coordinates. item=', trim(item)
             out_prelev = .false.
          else
             calc_pressure = .true.
          endif
       endif

       if ( out_vintrpl ) then
          if ( ktype /= '3D' ) then
             out_vintrpl = .false.
          endif
       endif

       item_save         (n) = item
       file_save         (n) = file
       desc_save         (n) = desc
       unit_save         (n) = unit
       step_save         (n) = step
       ktype_save        (n) = ktype
       kstr_save         (n) = kstr
       kend_save         (n) = kend
       kmax_save         (n) = kmax
       output_type_save  (n) = output_type
       out_prelev_save   (n) = out_prelev
       out_vintrpl_save  (n) = out_vintrpl
       opt_wgrid_save    (n) = opt_wgrid
       opt_lagintrpl_save(n) = opt_lagintrpl

       lname_save        (n) = lname
       tmax_save         (n) = 0
       tstr_save         (n) = TIME_CTIME
       tend_save         (n) = 0.D0

       month_old         (n) = idate(2)
       l_region_save     (n) = 0
    enddo

    allocate( ksumstr(HIST_req_nmax) )
    allocate( ksumend(HIST_req_nmax) )

    ksum = 0
    do n = 1, HIST_req_nmax
       ksumstr(n) = ksum + 1
       ksumend(n) = ksum + kmax_save(n)
       ksum       = ksum + kmax_save(n)
    enddo

    allocate( tsum_save(HIST_req_nmax,ADM_lall) )
    tsum_save(:,:) = 0.D0

    ! k-merged history container
    allocate( v_save   (ADM_gall,   ksum,ADM_lall,   1) )
    allocate( v_save_pl(ADM_gall_pl,ksum,ADM_lall_pl,1) )
    v_save   (:,:,:,:) = 0.D0
    v_save_pl(:,:,:,:) = 0.D0

    allocate( zlev_save(ksum) )

    do n = 1, HIST_req_nmax
       select case( trim(ktype_save(n)) )
       case('3D')
          if ( out_prelev_save(n) ) then
             zlev_save( ksumstr(n):ksumend(n) ) = pres_levs(1:kmax_save(n))
          else
             zlev_save( ksumstr(n):ksumend(n) ) = GRD_gz(ADM_kmin:ADM_kmax)
          endif
       case('2D')
          zlev_save( ksumstr(n):ksumend(n) ) = 0.D0
!       case('ISCCP')
!          do k = 1, NTAU_ISCCP*NPRES_ISCCP
!             zlev_save( ksumstr(n) + k-1 ) = real(k,kind=8)
!          enddo
!       case('GL')
!       case('GO')
!          do k = 1, kmax_save(n)
!             zlev_save( ksumstr(n) + k-1 ) = real(k,kind=8)
!          enddo
       case default
          zlev_save( ksumstr(n):ksumend(n) ) = GRD_gz( ADM_kmin+kstr_save(n)-1:ADM_kmin+kend_save(n)-1 )
       endselect
    enddo

    allocate( pres_levs_ln(npreslev) )
    pres_levs_ln(1:npreslev) = log( pres_levs(1:npreslev) * 100 )  

    allocate( cnvpre_klev(ADM_gall,npreslev,ADM_lall) )
    allocate( cnvpre_fac1(ADM_gall,npreslev,ADM_lall) )
    allocate( cnvpre_fac2(ADM_gall,npreslev,ADM_lall) )

    return
  end subroutine history_setup

  !-----------------------------------------------------------------------------
  subroutine  history_in( item, gd, l_region )
    use mod_adm, only : &
       ADM_proc_stop,   &
       ADM_l_me,        &
       ADM_gall,        &
       ADM_gall_in,     &
       ADM_lall,        &
       ADM_kmin,        &
       ADM_IopJop_nmax, &
       ADM_IopJop,      &
       ADM_GIoJo
    use mod_cnst, only: &
       CNST_UNDEF
    use mod_time, only: &
       TIME_CSTEP, &
       TIME_DTL
    use mod_calendar, only: &
       calendar_ss2yh
    implicit none

    character(len=*), intent(in) :: item
    real(8),          intent(in) :: gd(:,:)
    integer,          intent(in), optional :: l_region

    character(len=ADM_NSYS) :: hitem
    integer                 :: ijdim_input
    integer                 :: kdim_input

    logical :: save_var
    integer :: kmax
    integer :: g, g2, k, k2, l, n
    !---------------------------------------------------------------------------

    hitem = trim(item)

    ijdim_input = size(gd,1)
    kdim_input  = size(gd,2)

    if (       ijdim_input /= ADM_gall_in &
         .AND. ijdim_input /= ADM_gall    ) then
       write(ADM_LOG_FID,*) '+++ Module[history]/Category[nhm share]'
       write(ADM_LOG_FID,*) 'xxx invalid dimension, item=', hitem, &
                            ', ijdim_input=', ijdim_input, &
                            ', ADM_gall_in=', ADM_gall_in, &
                            ', ADM_gall=',    ADM_gall
       call ADM_proc_stop
    endif

    if ( calc_pressure ) then
       call get_log_pres
    endif

    do n = 1, HIST_req_nmax

       save_var = .false.

       ! Item is required or not? ( Same name can be contained in item_save )
       if ( hitem == item_save(n) ) then

          flag_save(n) = .true.

          if ( ktype_save(n) == '3D' ) then ! trim HALO
             if ( kdim_input-2 /= kmax_save(n) ) then
                write(ADM_LOG_FID,*) '+++ Module[history]/Category[nhm share]'
                write(ADM_LOG_FID,*) '*** Size unmatch, item=', hitem, &
                                     ', kdim_input=', kdim_input, &
                                     ', kmax_save=',  kmax_save(n)
             endif
          else
             if ( kdim_input /= kmax_save(n) ) then
                write(ADM_LOG_FID,*) '+++ Module[history]/Category[nhm share]'
                write(ADM_LOG_FID,*) '*** Size unmatch, item=', hitem, &
                                     ', kdim_input=', kdim_input, &
                                     ', kmax_save=',  kmax_save(n)
             endif
          endif

          ! add data or not?
          if ( trim(output_type_save(n)) == 'SNAPSHOT' ) then
             if( mod(TIME_CSTEP+1,step_save(n)) == 0 ) save_var = .true.
          else
             save_var = .true.
          endif
       endif

       if ( save_var ) then

          if ( present(l_region) ) then
             l_region_save(n) = l_region
          elseif( ADM_l_me >= 1 .and. ADM_l_me <= ADM_lall ) then
             l_region_save(n) = ADM_l_me
          else
             l_region_save(n) = l_region_save(n) + 1
             if( l_region_save(n) > ADM_lall ) l_region_save(n) = 1 ! cyclic
          endif

          kmax = min( kmax_save(n), kdim_input )
          l    = l_region_save(n)

          if ( .NOT. out_prelev_save(n) ) then ! normal

             if ( ktype_save(n) == '3D' ) then ! trim HALO

                if (ijdim_input == ADM_gall_in) then
                   do k = 1, kmax
                   do g = 1, ADM_IopJop_nmax
                      k2 = ksumstr(n)-1 + k
                      g2 = ADM_IopJop(g,ADM_GIoJo)

                      v_save(g2,k2,l,1) = v_save(g2,k2,l,1) + gd(g,k+ADM_kmin-1) * TIME_DTL
                   enddo
                   enddo
                else ! ijdim_input == ADM_gall
                   do k = 1, kmax
                   do g = 1, ADM_gall
                      k2 = ksumstr(n)-1 + k

                      v_save(g,k2,l,1) = v_save(g,k2,l,1) + gd(g,k+ADM_kmin-1) * TIME_DTL
                   enddo
                   enddo
                endif

             else

                if (ijdim_input == ADM_gall_in) then
                   do k = 1, kmax
                   do g = 1, ADM_IopJop_nmax
                      k2 = ksumstr(n)-1 + k
                      g2 = ADM_IopJop(g,ADM_GIoJo)

                      v_save(g2,k2,l,1) = v_save(g2,k2,l,1) + gd(g,k) * TIME_DTL
                   enddo
                   enddo
                else ! ijdim_input == ADM_gall
                   do k = 1, kmax
                   do g = 1, ADM_gall
                      k2 = ksumstr(n)-1 + k

                      v_save(g,k2,l,1) = v_save(g,k2,l,1) + gd(g,k) * TIME_DTL
                   enddo
                   enddo
                endif

             endif

          else ! convert to pressure level

             if (ijdim_input == ADM_gall_in) then
                do k2 = 1, npreslev
                do g  = 1, ADM_IopJop_nmax
                   k  = cnvpre_klev(g2,k,l)
                   g2 = ADM_IopJop(g,ADM_GIoJo)

                   if ( k > ADM_kmin ) then
                      v_save(g2,k2,l,1) = ( cnvpre_fac1(g2,k2,l) * gd(g,k-1) &
                                          + cnvpre_fac2(g2,k2,l) * gd(g,k  ) ) * TIME_DTL
                   else
                      v_save(g2,k2,l,1) = CNST_UNDEF
                   endif
                enddo
                enddo
             else ! ijdim_input == ADM_gall
                do k = 1, kmax
                do g = 1, ADM_gall
                   k2 = ksumstr(n)-1 + k

                   if ( k > ADM_kmin ) then
                      v_save(g,k2,l,1) = ( cnvpre_fac1(g,k2,l) * gd(g,k-1) &
                                         + cnvpre_fac2(g,k2,l) * gd(g,k  ) ) * TIME_DTL
                   else
                      v_save(g,k2,l,1) = CNST_UNDEF
                   endif
                enddo
                enddo
             endif

          endif ! z*-level or p-level

          tsum_save(n,l) = tsum_save(n,l) + TIME_DTL

       endif ! save data ?
    enddo

    return
  end subroutine history_in

  !----------------------------------------------------------------------------
  subroutine history_out
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_kall,      &
       ADM_kmax,      &
       ADM_kmin
    use mod_time, only: &
       TIME_CSTEP, &
       TIME_CTIME
    use mod_calendar, only : &
       calendar_ss2yh, &
       Calendar_SS2CC
    use mod_comm, only : &
       COMM_data_transfer, &
       COMM_var
    use mod_fio, only: &
       FIO_output
    use mod_gtl, only: &
       GTL_max, &
       GTL_min, &
       GTL_output_var2, &
       GTL_output_var2_da
    use mod_vintrpl, only: &
       VINTRPL_z_level, &
       VINTRPL_z_level2

    real(8) :: tmp   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: tmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: val_max, val_min

    character(len=20)           :: HTIME
    character(len=ADM_NSYS)     :: item
    character(len=ADM_MAXFNAME) :: basename

    logical, save :: first = .true.

    integer :: idate(6)   
    logical :: out_var(HIST_req_limit)
    integer :: num_output
    integer :: g, k, l, n
    !---------------------------------------------------------------------------

    if ( first ) then
       call history_outlist
       first = .false.
    endif

    call calendar_ss2yh( idate, TIME_CTIME )

    ! count up output vars at this time
    out_var(:) = .false.
    num_output = 0
    do n = 1, HIST_req_nmax
       if ( flag_save(n) ) then
          if ( trim(output_type_save(n)) == 'MONTHLY_AVERAGE' ) then
             if ( idate(2) /= month_old(n) ) then
                out_var(n)   = .true.
                month_old(n) = idate(2)

                num_output = num_output + 1
             endif
          else
             if ( mod( TIME_CSTEP, step_save(n) ) == 0 ) then
                out_var(n) = .true.

                num_output = num_output + 1
             endif
          endif
       endif
    enddo

    ! At least one variable will output, do communication
    if ( num_output > 0 ) then 
       write(ADM_LOG_FID,*) '### HISTORY num_output = ', num_output
       call Calendar_SS2CC ( HTIME, TIME_CTIME )
       write(ADM_LOG_FID,*) '###         date-time  = ', HTIME

       call comm_var( v_save, v_save_pl, KSUM, 1, comm_type=2, NSval_fix=.true. )
    else
       return
    endif

    do n = 1, HIST_req_nmax

       if ( out_var(n) ) then

          tmax_save(n) = tmax_save(n) + 1
          tend_save(n) = TIME_CTIME

          do l = 1, ADM_lall
          do k = ksumstr(n), ksumend(n)
          do g = 1, ADM_gall
             v_save(g,k,l,1) = v_save(g,k,l,1) / tsum_save(n,l)
          enddo
          enddo
          enddo

          do l = 1, ADM_lall_pl
          do k = ksumstr(n), ksumend(n)
          do g = 1, ADM_gall_pl
             v_save_pl(g,k,l,1) = v_save_pl(g,k,l,1) / tsum_save(n,1)
          enddo
          enddo
          enddo

          item = item_save(n)
          val_max =  GTL_max( v_save   (:,:,:,1),          &
                              v_save_pl(:,:,:,1),          &
                              ksum, ksumstr(n), ksumend(n) )
          val_min =  GTL_min( v_save   (:,:,:,1),          &
                              v_save_pl(:,:,:,1),          &
                              ksum, ksumstr(n), ksumend(n) )

          if ( out_vintrpl_save(n) ) then
             do k = Adm_kmin, Adm_kmax
                tmp   (:,k,:) = v_save   (:,ksumstr(n)+k-2,:,1)
                tmp_pl(:,k,:) = v_save_pl(:,ksumstr(n)+k-2,:,1)
             enddo

             tmp   (:,ADM_kmin-1,:) = tmp   (:,ADM_kmin,:)
             tmp   (:,ADM_kmax+1,:) = tmp   (:,ADM_kmax,:)
             tmp_pl(:,ADM_kmin-1,:) = tmp_pl(:,ADM_kmin,:)
             tmp_pl(:,ADM_kmax+1,:) = tmp_pl(:,ADM_kmax,:)

             if ( opt_lagintrpl_save(n) ) then
                call VINTRPL_z_level ( tmp, tmp_pl, opt_wgrid_save(n) )
             else
                call VINTRPL_z_level2( tmp, tmp_pl, opt_wgrid_save(n) )
             endif

             do k = ADM_kmin, ADM_kmax
                v_save   (:,ksumstr(n)+k-2,:,1) = tmp   (:,k,:)
                v_save_pl(:,ksumstr(n)+k-2,:,1) = tmp_pl(:,k,:)
             enddo
          endif

          write(ADM_LOG_FID,'(A,A16,A,1PE24.17,A,E24.17)') ' [', item(1:16), '] max=', val_max, ', min=', val_min

          if ( trim(output_io_mode) == 'ADVANCED' ) then

             if ( trim(output_type_save(n)) == 'SNAPSHOT' ) then

                call FIO_output( v_save(:,:,:,1),                             &
                                 HIST_io_fname,    HIST_io_desc    , '',      &
                                 file_save(n),     desc_save(n), '',          &
                                 unit_save(n),     HIST_dtype,                &
                                 lname_save(n),    ksumstr(n),   ksumend(n),  &
                                 tmax_save(n),     tend_save(n), tend_save(n) )

             elseif(trim(output_type_save(n)) == 'AVERAGE') then

                call FIO_output( v_save(:,:,:,1),                             &
                                 HIST_io_fname,    HIST_io_desc    , '',      &
                                 file_save(n),     desc_save(n), '',          &
                                 unit_save(n),     HIST_dtype,                &
                                 lname_save(n),    ksumstr(n),   ksumend(n),  &
                                 tmax_save(n),     tstr_save(n), tend_save(n) )

             endif

          elseif( trim(output_io_mode) == 'LEGACY' ) then
             basename = trim(output_path)//file_save(n)

             if ( direct_access ) then
                call GTL_output_var2_da( basename, v_save(:,:,:,1),                        &
                                         ksumstr(n), ksumend(n), tmax_save(n), output_size )
             else
                call GTL_output_var2( basename, v_save(:,:,:,1), v_save_pl(:,:,:,1), &
                                      ksumstr(n), ksumend(n), output_size            )
             endif

             call history_timeinfo
          endif

          ! reset saved variable
          v_save   (:,ksumstr(n):ksumend(n),:,1) = 0.D0
          v_save_pl(:,ksumstr(n):ksumend(n),:,1) = 0.D0

          tstr_save(n) = TIME_CTIME
          tsum_save(n,:) = 0.D0
       endif
    enddo

  end subroutine history_out

  !-----------------------------------------------------------------------------
  subroutine history_outlist
    use mod_adm, only: &
       ADM_proc_stop
    implicit none

    character(len=ADM_NSYS)     :: item
    character(len=ADM_MAXFNAME) :: file
    character(len=ADM_NSYS)     :: unit
    character(len=ADM_NSYS)     :: ktype
    character(len=ADM_NSYS)     :: otype

    integer :: n
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '*** [HIST] Output item list '
    write(ADM_LOG_FID,*) '*** Total number of requested history item :', HIST_req_nmax
    write(ADM_LOG_FID,*) '============================================================================'
    write(ADM_LOG_FID,*) 'NAME            :Save name       :UNIT            :Avg.type        :interval'
    write(ADM_LOG_FID,*) '                :Vert.type       :# of layer      :p?  :z?  :zh? :lag.intrp?'
    write(ADM_LOG_FID,*) '============================================================================'

    do n = 1, HIST_req_nmax
       item  = item_save(n)
       file  = file_save(n)
       unit  = unit_save(n)
       ktype = ktype_save(n)
       otype = output_type_save(n)

       write(ADM_LOG_FID,'(1x,A16,A,A16,A,A16,A,A16,A,I8)')      item (1:16), &
                                                            ":", file (1:16), &
                                                            ":", unit (1:16), &
                                                            ":", otype(1:16), &
                                                            ":", step_save(n)

       write(ADM_LOG_FID,'(17x,A,A16,A,I016,A,L04,A,L04,A,L04,A,L04)') ":", ktype(1:16),           &
                                                                       ":", kmax_save(n),          &
                                                                       ":", out_prelev_save   (n), &
                                                                       ":", out_vintrpl_save  (n), &
                                                                       ":", opt_wgrid_save    (n), &
                                                                       ":", opt_lagintrpl_save(n)

       if ( .NOT. flag_save(n) ) then ! not stored yet or never
          write(ADM_LOG_FID,*) '+++ this variable is requested but not stored yet. check!'
          if ( check_flag ) then
             write(ADM_LOG_FID,*) 'xxx history check_flag is on. stop!'
             call ADM_proc_stop
          endif
       endif
    enddo

    write(ADM_LOG_FID,*) '============================================================================'
    write(ADM_LOG_FID,*)

    return
  end subroutine history_outlist

  !-----------------------------------------------------------------------------
  subroutine history_timeinfo
    use mod_misc, only: &
       MISC_get_available_fid
    use mod_adm, only: &
       ADM_prc_me,         &
       ADM_prc_run_master
    use mod_time, only: &
       TIME_DTL
    implicit none

    integer :: fid
    integer :: n, k
    !---------------------------------------------------------------------------

    if ( ADM_prc_me == ADM_prc_run_master ) then
       fid = MISC_get_available_fid()
       open( unit   = fid,                               &
             file   = trim(output_path)//'history.info', &
             form   = 'formatted',                       &
             status = 'replace'                          )

          do n = 1, HIST_req_nmax
             write(fid,'(I8,F16.2)') tmax_save(n), step_save(n)*TIME_DTL
             write(fid,'(I8)')       kmax_save(n)
             do k = 1, kmax_save(n)
                write(fid,'(F16.4)') zlev_save(ksumstr(n)+k-1)
             enddo
             write(fid,'(I8)')       1
             write(fid,'(A32)')      trim(file_save(n))
          enddo

       close(fid)
    endif

    return
  end subroutine history_timeinfo

  !-----------------------------------------------------------------------------
  subroutine get_log_pres  ! Y.Niwa add 071130
    !
    use mod_adm, only : &
         ADM_gall,  &
         ADM_kall,  &
         ADM_kmax,  &
         ADM_kmin,  &
         ADM_lall,  &
         ADM_GALL_PL, &
         ADM_LALL_PL, &
         ADM_KNONE
    use mod_cnvvar, only :  &
         cnvvar_p2d
    use mod_runconf, only : &
         nqmax => TRC_VMAX
    use mod_thrmdyn, only :    &
         thrmdyn_tempre,   &
         thrmdyn_qd
    use mod_vmtr, only :        &
         VMTR_GSGAM2
    use mod_sfcvar, only :    &
         sfcvar_get,          &
         I_PRE_SFC
    use mod_prgvar, only :     &
         prgvar_get
    !
    implicit none
    !
    integer :: l, nq, ij, kk, k
    real(8) :: rhog(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhoge(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogq(ADM_gall,ADM_kall,ADM_lall,nqmax)
    real(8) :: rhogq_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nqmax)
    !
    real(8) :: q(ADM_gall,ADM_kall,nqmax)
    real(8) :: rho(ADM_gall,ADM_kall)
    real(8) :: pre(ADM_gall,ADM_kall)
    real(8) :: tem(ADM_gall,ADM_kall)
    real(8) :: rgrho(ADM_gall,ADM_kall)
    real(8) :: ein(ADM_gall,ADM_kall)
    real(8) :: qd(ADM_gall,ADM_kall)
    real(8) :: v2d_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL)
    real(8) :: v3d_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d2_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d3_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d4_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d5_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d6_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: pre_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    integer, parameter :: num = 0 

    real(8) :: lpres_sfc(ADM_gall,ADM_lall)
    real(8) :: lpres(ADM_gall,ADM_kall)

    call prgvar_get(         &
         rhog,   v3d_pl,   &  !--- out
         rhogvx, v3d2_pl,  &  !--- out
         rhogvy, v3d3_pl,  &  !--- out
         rhogvz, v3d4_pl,  &  !--- out
         rhogw,  v3d5_pl,  &  !--- out
         rhoge,  v3d6_pl,  &  !--- out
         rhogq,  rhogq_pl, &  !--- out
         num )                !--- in

    call sfcvar_get(           &
         pre_sfc, v2d_pl,      &  !--- out
         vid = I_PRE_SFC       &  !--- in
         )
    lpres_sfc(:,:) = log(pre_sfc(:,ADM_KNONE,:))

    cnvpre_fac1(:,:,:) = 0.D0
    cnvpre_fac2(:,:,:) = 0.D0
    cnvpre_klev(:,:,:) = -1

    !--- calculation of tem, pre, and th
    do l=1, ADM_lall

       rgrho(:,:) = 1.0d0/rhog(:,:,l)
       ein(:,:) = rhoge(:,:,l) * rgrho(:,:)

       do nq=1, nqmax
          q(:,:,nq) = rhogq(:,:,l,nq) * rgrho(:,:)
       enddo

       call thrmdyn_qd( ADM_gall, qd, q )

       rho(:,:) = rhog(:,:,l)/VMTR_GSGAM2(:,:,l)
       call thrmdyn_tempre( &
            ADM_gall,       &
            tem,            &  !--- out
            pre,            &  !--- out
            ein,            &  !--- in
            rho,            &  !--- in
            qd,             &  !--- in
            q )                !--- in

       lpres(:,:) = log(pre(:,:))

       do kk=1, npreslev
          do ij=1, ADM_gall
             if(lpres_sfc(ij,l) < pres_levs_ln(kk)) then
                !cnvpre_klev(ij,kk,l) = -1
             else
                do k=ADM_kmin, ADM_kmax
                   if( pres_levs_ln(kk) > lpres(ij,k) ) exit
                enddo

                if( k == ADM_kmin )     k = ADM_kmin + 1 ! extrapolation
                if( k == ADM_kmax + 1 ) k = ADM_kmax     ! extrapolation

                cnvpre_klev(ij,kk,l) = k
                cnvpre_fac1(ij,kk,l) = ( lpres(ij,k) - pres_levs_ln(kk) ) / ( lpres(ij,k) - lpres(ij,k-1) )
                cnvpre_fac2(ij,kk,l) = ( pres_levs_ln(kk) - lpres(ij,k-1) ) / ( lpres(ij,k) - lpres(ij,k-1) )

             endif
          enddo
       enddo

    enddo

    return
  end subroutine get_log_pres

end module mod_history
!-------------------------------------------------------------------------------
