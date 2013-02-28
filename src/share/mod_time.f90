!-------------------------------------------------------------------------------
!>
!! Time management module
!!
!! @par Description
!!         This module is for the time management.
!! @author  H.Tomita
!! @par History
!! @li      2004-02-17 (H.Tomita) Imported from igdc-4.33
!! @li      2004-05-31 (      )   Calculation of "num_of_iteration_[sl]step" are moved to mod[onestep].
!<
module mod_time
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TIME_setup
  public :: TIME_report
  public :: TIME_advance

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !--- Integration method in large steps
  character(ADM_NSYS), public, save :: TIME_INTEG_TYPE = 'RK2' ! Runge-Kutta 2nd
  !                                                    = 'RK3'   Runge-Kutta 3rd
  !
  !--- Horizontally spliting or not
  logical, public, save :: TIME_SPLIT = .TRUE.
  !
  !--- Time interval for large step
  real(8), public, save :: TIME_DTL = 5.0D0
  !
  !--- Time interval for small step ( calculated in sub[TIME_setup] )
  real(8), public, save :: TIME_DTS
  !
  !--- Start time
  real(8), public, save :: TIME_START = 0.0D0
  !
  !--- End time ( calculated in sub[TIME_setup] )
  real(8), public, save :: TIME_END
  !
  !--- Max steps of large step
  integer, public, save :: TIME_LSTEP_MAX = 10
  !
  !--- Max steps of small step
  integer, public, save :: TIME_SSTEP_MAX = 10
  !
  !--- Number of initial steps
  integer, public, save :: TIME_NUM_OF_INITIAL_STEPS = 1
  !
  !--- Current time step ( calculated in sub[TIME_setup] )
  integer, public, save :: TIME_CSTEP
  !
  !--- Current time ( calculated in sub[TIME_setup] )
  real(8), public, save :: TIME_CTIME
  !
  !--- Time step at the start ( calculated in sub[TIME_setup] )
  integer, public, save :: TIME_NSTART
  !
  !--- Time step at the end ( calculated in sub[TIME_setup] )
  integer, public, save :: TIME_NEND
  !
  !--- Time filter for leapfrog ( if alpha < 0: no use of filter )
  real(8), public, save :: TIME_ALPHA_TIME_FILTER = 0.05D0

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
  !>
  !> Setup the temporal scheme and time management
  !>
  subroutine TIME_setup
    use mod_adm, only: &
       ADM_CTL_FID, &
       ADM_proc_stop
    use mod_calendar, only: &
       calendar_yh2ss, &
       calendar_ss2cc
    implicit none

    character(len=ADM_NSYS) :: integ_type

    integer :: start_year
    integer :: start_month
    integer :: start_day
    integer :: start_hour
    integer :: start_min
    integer :: start_sec

    real(8) :: dtl
    integer :: lstep_max
    integer :: sstep_max
    integer :: num_of_initial_steps
    real(8) :: alpha_time_filter
    logical :: split

    namelist / TIMEPARAM /     &
         integ_type,           &  !--- integration method
         dtl,                  &  !--- delta t in large step
         lstep_max,            &  !--- maximum number of large steps
         sstep_max,            &  !--- division number in large step
         num_of_initial_steps, &  !--- initial step number
         alpha_time_filter,    &  !--- time filter for leapfrog
         split,                &  !--- time spliting flag
         start_year,           &  !--- start year
         start_month,          &  !--- start month
         start_day,            &  !--- start day
         start_hour,           &  !--- start hour
         start_min,            &  !--- start min
         start_sec                !--- start sec

    integer           :: start_date(6)
    character(len=20) :: HTIME_start
    character(len=20) :: HTIME_end

    integer :: ierr
    !---------------------------------------------------------------------------

    integ_type           = TIME_integ_type
    dtl                  = TIME_dtl
    lstep_max            = TIME_lstep_max
    sstep_max            = TIME_sstep_max
    num_of_initial_steps = TIME_num_of_initial_steps
    alpha_time_filter    = TIME_alpha_time_filter
    split                = TIME_split

    start_year =  0
    start_month = 1
    start_day   = 1
    start_hour  = 0
    start_min   = 0
    start_sec   = 0

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[time]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=TIMEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** TIMEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist TIMEPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist TIMEPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,TIMEPARAM)

    !--- rewrite
    TIME_integ_type = integ_type
    TIME_dtl                  = dtl
    TIME_lstep_max            = lstep_max
    TIME_sstep_max            = sstep_max
    TIME_num_of_initial_steps = num_of_initial_steps
    TIME_alpha_time_filter    = alpha_time_filter
    TIME_split                = split

    start_date(1) = start_year
    start_date(2) = start_month
    start_date(3) = start_day
    start_date(4) = start_hour
    start_date(5) = start_min
    start_date(6) = start_sec
    call calendar_yh2ss( TIME_start, start_date )

    if ( TIME_LSTEP_MAX > 0 ) then
       TIME_END = TIME_START + TIME_LSTEP_MAX * TIME_DTL
    endif

    call calendar_ss2cc ( HTIME_start, TIME_START )
    call calendar_ss2cc ( HTIME_end,   TIME_END   )

    !------ set the step number at the intial step
    TIME_NSTART = 0

    !------ set the step number at the final step
    TIME_NEND = TIME_NSTART + TIME_LSTEP_MAX

    !------ set the current step
    TIME_CSTEP = TIME_NSTART

    !------ intitialize the current time
    TIME_CTIME = TIME_START

    !--- < small step configuration > ---
    if ( TIME_INTEG_TYPE == 'LEAPFROG' ) then
       TIME_DTS = 2.0D0 * TIME_DTL / dble(TIME_SSTEP_MAX)
    else
       TIME_DTS = TIME_DTL / dble(TIME_SSTEP_MAX)
    endif

    !------ Fix TIME_NUM_OF_INITIAL_STEPS for consistency with
    !------ initial step method.
    !------ If TIME_INTEG_TYPE=='LEAPFROG',
    !------ the initial conditions for two steps are required.
    !
    if ( TIME_NUM_OF_INITIAL_STEPS < 0 ) then
       if ( TIME_INTEG_TYPE == 'LEAPFROG' ) then
          TIME_NUM_OF_INITIAL_STEPS = 1
       else
          TIME_NUM_OF_INITIAL_STEPS = 0
       endif
    endif

    !--- output the information for debug
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '====== Time management ======'
    write(ADM_LOG_FID,*) '--- Time integration scheme (large step): ', trim(TIME_INTEG_TYPE)
    write(ADM_LOG_FID,*) '--- Time interval for large step        : ', TIME_DTL
    write(ADM_LOG_FID,*) '--- Time interval for small step        : ', TIME_DTS
    write(ADM_LOG_FID,*) '--- Max steps of large step             : ', TIME_LSTEP_MAX
    write(ADM_LOG_FID,*) '--- Max steps of small step             : ', TIME_SSTEP_MAX
    write(ADM_LOG_FID,*) '--- Start time (sec)                    : ', TIME_START
    write(ADM_LOG_FID,*) '--- End time   (sec)                    : ', TIME_END
    write(ADM_LOG_FID,*) '--- Start time (date)                   : ', HTIME_start
    write(ADM_LOG_FID,*) '--- End time   (date)                   : ', HTIME_end
    write(ADM_LOG_FID,*) '--- total integration time              : ', TIME_END - TIME_START
    write(ADM_LOG_FID,*) '--- Time step at the start              : ', TIME_NSTART
    write(ADM_LOG_FID,*) '--- Time step at the end                : ', TIME_NEND
    write(ADM_LOG_FID,*) '--- Number of initial steps             : ', TIME_NUM_OF_INITIAL_STEPS
    write(ADM_LOG_FID,*) '--- Time filter for leapfrog            : ', TIME_ALPHA_TIME_FILTER

    return
  end subroutine TIME_setup

  !-----------------------------------------------------------------------------
  subroutine TIME_report
    use mod_adm, only: &
       ADM_prc_run_master, &
       ADM_prc_me
    use mod_calendar, only: &
       calendar_ss2cc
    implicit none

    character(len=20) :: HTIME
    !---------------------------------------------------------------------------

    call calendar_ss2cc ( HTIME, TIME_CTIME )

    write(ADM_LOG_FID,*) '### TIME =', HTIME,'( step = ', TIME_CSTEP, ' )' 
    if( ADM_prc_me == ADM_prc_run_master ) then
       write(*,*) '### TIME = ', HTIME,'( step = ', TIME_CSTEP, ' )' 
    endif

    return
  end subroutine TIME_report

  !-----------------------------------------------------------------------------
  subroutine TIME_advance
    implicit none
    !---------------------------------------------------------------------------

    ! time advance
    TIME_CTIME = TIME_CTIME + TIME_DTL
    TIME_CSTEP = TIME_CSTEP + 1

    call TIME_report

    return
  end subroutine TIME_advance

end module mod_time
!-------------------------------------------------------------------------------
