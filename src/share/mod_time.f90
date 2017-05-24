!-------------------------------------------------------------------------------
!> Module time management
!!
!! @par Description
!!          This module is for the time management
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_time
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: TIME_setup
  public :: TIME_report
  public :: TIME_advance

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: TIME_INTEG_TYPE = 'UNDEF'  ! Integration method in large steps
                                                    ! 'RK2'    ! Runge-Kutta 2nd
                                                    ! 'RK3'    ! Runge-Kutta 3rd
                                                    ! 'RK4'    ! Runge-Kutta 4th
                                                    ! 'TRCADV' ! Tracer advection only

  logical,  public :: TIME_SPLIT     = .true. ! Horizontally splitting?

  integer,  public :: TIME_LSTEP_MAX = 10     ! Max steps of large step
  integer,  public :: TIME_SSTEP_MAX          ! Max steps of small step

  real(DP), public :: TIME_DTL       = 5.0_DP ! Time interval for large step [sec]
  real(DP), public :: TIME_DTS                ! Time interval for small step [sec]
  !
  real(DP), public :: TIME_START              ! Start time [sec]
  real(DP), public :: TIME_END                ! End   time [sec]
  integer,  public :: TIME_NSTART             ! Time step at the start
  integer,  public :: TIME_NEND               ! Time step at the end

  real(DP), public :: TIME_CTIME              ! Current time [sec]
  integer,  public :: TIME_CSTEP              ! Current time step

  character(len=20), public :: TIME_HTIME     ! YYYY/MM/DD-HH:MM:SS

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: TIME_backward_sw = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup the temporal scheme and time management
  subroutine TIME_setup( &
       backward )
    use mod_process, only: &
       PRC_MPIstop
    use mod_calendar, only: &
       CALENDAR_yh2ss, &
       CALENDAR_ss2cc
    implicit none

    logical, intent(in), optional :: backward ! backward option (optional) [TM]

    character(len=H_SHORT) :: integ_type    !< integration method
    logical                :: split         !< time spliting flag
    real(DP)               :: dtl           !< delta t in large step
    integer                :: lstep_max     !< maximum number of large steps
    integer                :: sstep_max     !< division number in large step
    integer                :: start_date(6) !< start date
    integer                :: start_year    !< start year
    integer                :: start_month   !< start month
    integer                :: start_day     !< start day
    integer                :: start_hour    !< start hour
    integer                :: start_min     !< start min
    integer                :: start_sec     !< start sec

    namelist / TIMEPARAM / &
         integ_type,  &
         split,       &
         dtl,         &
         lstep_max,   &
         sstep_max,   &
         start_date,  &
         start_year,  &
         start_month, &
         start_day,   &
         start_hour,  &
         start_min,   &
         start_sec

    character(len=20) :: HTIME_start
    character(len=20) :: HTIME_end

    integer  :: ierr
    !---------------------------------------------------------------------------

    integ_type = TIME_integ_type
    split      = TIME_split
    dtl        = TIME_dtl
    lstep_max  = TIME_lstep_max
    sstep_max  = -999

    start_date(:) = -999
    start_year    = 0
    start_month   = 1
    start_day     = 1
    start_hour    = 0
    start_min     = 0
    start_sec     = 0

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[time]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=TIMEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** TIMEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist TIMEPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=TIMEPARAM)

    !--- rewrite
    TIME_integ_type = integ_type
    TIME_split      = split
    TIME_dtl        = dtl
    TIME_lstep_max  = lstep_max

    if ( sstep_max == -999 )  then
       if( IO_L ) write(IO_FID_LOG,*) 'TIME_integ_type is ', trim(TIME_integ_type)
       select case(TIME_integ_type)
       case('RK2')
          TIME_sstep_max = 4
       case('RK3')
          TIME_sstep_max = 6
        case('RK4')
          TIME_sstep_max = 8
        case('TRCADV')
          TIME_sstep_max = 0
       case default
          write(*,*) 'xxx Invalid TIME_INTEG_TYPE! STOP.'
       endselect
       if( IO_L ) write(IO_FID_LOG,*) 'TIME_sstep_max is automatically set to: ', TIME_sstep_max
    else
       TIME_sstep_max = sstep_max
    endif
    TIME_dts = TIME_dtl / max(real(TIME_sstep_max,kind=DP),1.0_DP)

    if ( start_date(1) == -999 ) start_date(1) = start_year
    if ( start_date(2) == -999 ) start_date(2) = start_month
    if ( start_date(3) == -999 ) start_date(3) = start_day
    if ( start_date(4) == -999 ) start_date(4) = start_hour
    if ( start_date(5) == -999 ) start_date(5) = start_min
    if ( start_date(6) == -999 ) start_date(6) = start_sec
    call CALENDAR_yh2ss( TIME_start, start_date )

    if ( present(backward) ) then
       TIME_backward_sw = backward ! [TM]
    else
       TIME_backward_sw = .false.
    endif

    !---< large step configuration >---

    if ( TIME_lstep_max < 0 ) then
       write(*,*) 'xxx TIME_lstep_max should be positive. STOP.'
       call PRC_MPIstop
    endif
    if ( TIME_sstep_max < 0 ) then
       write(*,*) 'xxx TIME_sstep_max should be positive. STOP.'
       call PRC_MPIstop
    endif
    if ( TIME_dtl < 0 ) then
       write(*,*) 'xxx TIME_dtl should be positive. STOP.'
       call PRC_MPIstop
    endif

    if ( .NOT. TIME_backward_sw ) then
       TIME_END = TIME_START + TIME_LSTEP_MAX * TIME_DTL
    else
       TIME_END = TIME_START - TIME_LSTEP_MAX * TIME_DTL ! [TM]
    endif

    TIME_NSTART = 0
    TIME_NEND   = TIME_NSTART + TIME_LSTEP_MAX

    TIME_CTIME  = TIME_START
    TIME_CSTEP  = TIME_NSTART

    !--- output the information for debug
    call CALENDAR_ss2cc( HTIME_start, TIME_START )
    call CALENDAR_ss2cc( HTIME_end,   TIME_END   )
    call CALENDAR_ss2cc( TIME_HTIME,  TIME_CTIME )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '====== Time management ======'
    if( IO_L ) write(IO_FID_LOG,*) '--- Time integration scheme (large step): ', trim(TIME_integ_type)
    if( IO_L ) write(IO_FID_LOG,*) '--- Backward integration?               : ', TIME_backward_sw
    if( IO_L ) write(IO_FID_LOG,*) '--- Time interval for large step        : ', TIME_DTL
    if( IO_L ) write(IO_FID_LOG,*) '--- Time interval for small step        : ', TIME_DTS
    if( IO_L ) write(IO_FID_LOG,*) '--- Max steps of large step             : ', TIME_LSTEP_MAX
    if( IO_L ) write(IO_FID_LOG,*) '--- Max steps of small step             : ', TIME_SSTEP_MAX
    if( IO_L ) write(IO_FID_LOG,*) '--- Start time (sec)                    : ', TIME_START
    if( IO_L ) write(IO_FID_LOG,*) '--- End time   (sec)                    : ', TIME_END
    if( IO_L ) write(IO_FID_LOG,*) '--- Start time (date)                   : ', HTIME_start
    if( IO_L ) write(IO_FID_LOG,*) '--- End time   (date)                   : ', HTIME_end
    if( IO_L ) write(IO_FID_LOG,*) '--- total integration time              : ', TIME_END - TIME_START
    if( IO_L ) write(IO_FID_LOG,*) '--- Time step at the start              : ', TIME_NSTART
    if( IO_L ) write(IO_FID_LOG,*) '--- Time step at the end                : ', TIME_NEND

    return
  end subroutine TIME_setup

  !-----------------------------------------------------------------------------
  subroutine TIME_report
    use mod_process, only: &
       PRC_IsMaster
    use mod_calendar, only: &
       CALENDAR_ss2cc
    implicit none
    !---------------------------------------------------------------------------

    call calendar_ss2cc( TIME_HTIME, TIME_CTIME )

    if( IO_L ) write(IO_FID_LOG,'(1x,3A,I8,A,I8,A)') &
                                '### TIME = ', TIME_HTIME, '( step = ', TIME_CSTEP, '/', TIME_LSTEP_MAX, ' )'
    if( PRC_IsMaster ) then
       write(*,'(1x,3A,I8,A,I8,A)') &
               '### TIME = ', TIME_HTIME, '( step = ', TIME_CSTEP, '/', TIME_LSTEP_MAX, ' )'
    endif

    return
  end subroutine TIME_report

  !-----------------------------------------------------------------------------
  subroutine TIME_advance
    implicit none
    !---------------------------------------------------------------------------

    ! time advance
    if ( .NOT. TIME_backward_sw ) then
       TIME_CTIME = TIME_CTIME + TIME_DTL
    else
       TIME_CTIME = TIME_CTIME - TIME_DTL ! [TM]
    endif
    TIME_CSTEP = TIME_CSTEP + 1

    call TIME_report

    return
  end subroutine TIME_advance

end module mod_time
