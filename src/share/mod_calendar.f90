!-------------------------------------------------------------------------------
!>
!! Calendar module
!!
!! @par Description
!!         This module provides the subroutines for calendar.
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_calendar
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
  public :: CALENDAR_setup
  public :: CALENDAR_ss2yh
  public :: CALENDAR_yh2ss
  public :: CALENDAR_ss2cc

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  public :: CALENDAR_ss2ds
  public :: CALENDAR_ds2ss
  public :: CALENDAR_rs2hm
  public :: CALENDAR_hm2rs
  public :: CALENDAR_dd2ym
  public :: CALENDAR_ym2dd

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: oauto = .true.    ! automatic setting of DOY?
                                        ! automatic means,
                                        ! yr =    0- 999 : 360day
                                        ! yr = 1000-1899 : 365day
                                        ! yr = 1900-     : gregorian

  logical, private :: ogrego = .true.   ! use gregorian calendar?

  logical, private :: oideal = .false.  ! use ideal calendar?
  integer, private :: imonyr = 12       ! # of months  for a year  (for ideal case)
  integer, private :: idaymo = 30       ! # of days    for a month (for ideal case)

  integer, private :: monday(12,2)      ! # of days    for a month
  data monday / 31,28,31,30,31,30,31,31,30,31,30,31, & ! normal year
                31,29,31,30,31,30,31,31,30,31,30,31  / ! leap   year

  integer, private :: ihrday = 24       ! # of hours   for a days
  integer, private :: iminhr = 60       ! # of minutes for a hour
  integer, private :: isecmn = 60       ! # of seconds for a minute

  logical, private :: operpt = .false.  ! use perpetual date?
  integer, private :: iyrpp  = 0        ! perpetual year
  integer, private :: imonpp = 3        ! perpetual month
  integer, private :: idaypp = 21       ! perpetual day

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine CALENDAR_setup
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / NM_CALENDAR / &
       oauto,  &
       ogrego, &
       oideal, &
       imonyr, &
       idaymo, &
       ihrday, &
       iminhr, &
       isecmn, &
       operpt, &
       iyrpp,  &
       imonpp, &
       idaypp

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[calendar]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=NM_CALENDAR,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** NM_CALENDAR is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist NM_CALENDAR. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=NM_CALENDAR)

    return
  end subroutine CALENDAR_setup

  !-----------------------------------------------------------------------------
  !> seconds -> date
  subroutine CALENDAR_ss2yh( &
       idate, &
       dsec   )
    implicit none

    integer,  intent(out) :: idate(6) ! date [yymmddhhmmss]
    real(DP), intent(in)  :: dsec     ! time [sec]

    integer  :: idays ! # of days
    real(DP) :: rsec  ! # of seconds in a day
    !---------------------------------------------------------------------------

    call CALENDAR_ss2ds( idays, rsec, & ! [OUT]
                         dsec         ) ! [IN]

    call CALENDAR_dd2ym( idate(1), idate(2), idate(3), & ! [OUT]
                         idays                         ) ! [IN]

    call CALENDAR_rs2hm( idate(4), idate(5), idate(6), & ! [OUT]
                         rsec                          ) ! [IN]

    return
  end subroutine CALENDAR_ss2yh

  !-----------------------------------------------------------------------------
  !> date -> seconds
  subroutine CALENDAR_yh2ss( &
       dsec, &
       idate )
    implicit none

    real(DP), intent(out) :: dsec     ! time [sec]
    integer,  intent(in)  :: idate(6) ! date [yymmddhhmmss]

    integer  :: idays ! # of days
    real(DP) :: rsec  ! # of seconds in a day
    !---------------------------------------------------------------------------

    call CALENDAR_ym2dd( idays,                       & ! [OUT]
                         idate(1), idate(2), idate(3) ) ! [IN]

    call CALENDAR_hm2rs( rsec,                        & ! [OUT]
                         idate(4), idate(5), idate(6) ) ! [IN]

    call CALENDAR_ds2ss( dsec,       & ! [OUT]
                         idays, rsec ) ! [IN]

    return
  end subroutine CALENDAR_yh2ss

  !-----------------------------------------------------------------------------
  !> seconds -> ddss
  subroutine CALENDAR_ss2ds( &
       idays, &
       rsec,  &
       dsec   )
    implicit none

    integer,  intent(out) :: idays
    real(DP), intent(out) :: rsec
    real(DP), intent(in)  :: dsec

    integer  :: isecdy ! # of seconds for a day
    !---------------------------------------------------------------------------

    isecdy = isecmn * iminhr * ihrday
    idays  = int( dsec / real(isecdy,kind=DP) ) + 1
    rsec   = dsec - real(idays-1,kind=DP) * real(isecdy,kind=DP)

    if ( nint(rsec) >= isecdy ) then
       idays = idays + 1
       rsec  = rsec  - real(isecdy,kind=DP)
    endif

    return
  end subroutine CALENDAR_ss2ds

  !-----------------------------------------------------------------------------
  !> ddss -> seconds
  subroutine CALENDAR_ds2ss( &
       dsec,  &
       idays, &
       rsec   )
    implicit none

    real(DP), intent(out) :: dsec
    integer,  intent(in)  :: idays
    real(DP), intent(in)  :: rsec

    integer :: isecdy ! # of seconds for a day
    !---------------------------------------------------------------------------

    isecdy = isecmn * iminhr * ihrday
    dsec   = real(idays-1,kind=DP) * real(isecdy,kind=DP) + real(rsec,kind=DP)

    return
  end subroutine CALENDAR_ds2ss

  !-----------------------------------------------------------------------------
  !> seconds -> hhmmss
  subroutine CALENDAR_rs2hm( &
       ihour, &
       imin,  &
       isec,  &
       rsec   )
    implicit none

    integer,  intent(out) :: ihour
    integer,  intent(out) :: imin
    integer,  intent(out) :: isec
    real(DP), intent(in)  :: rsec

    integer :: isechr ! # of seconds for a hour
    !---------------------------------------------------------------------------

    isechr = isecmn * iminhr
    ihour  = int ( rsec / real(isechr,kind=DP) )
    imin   = int ( ( rsec - real(ihour*isechr,kind=DP) ) / real(isecmn,kind=DP) )
    isec   = nint( rsec - real(ihour*isechr,kind=DP) - real(imin*isecmn,kind=DP) )

    if ( isec >= isecmn ) then
       imin  = imin + 1
       isec  = isec - isecmn
    endif

    if ( imin == iminhr ) then
       ihour = ihour + 1
       imin  = imin  - iminhr
    endif

    return
  end subroutine CALENDAR_rs2hm

  !-----------------------------------------------------------------------------
  !> hhmmss -> seconds
  subroutine CALENDAR_hm2rs( &
       rsec,  &
       ihour, &
       imin,  &
       isec   )
    implicit none

    real(DP), intent(out) :: rsec
    integer,  intent(in)  :: ihour
    integer,  intent(in)  :: imin
    integer,  intent(in)  :: isec
    !---------------------------------------------------------------------------

    rsec = real( ihour*isecmn*iminhr + imin*isecmn + isec, kind=DP )

    return
  end subroutine CALENDAR_hm2rs

  !-----------------------------------------------------------------------------
  !> day -> yymmdd
  subroutine CALENDAR_dd2ym( &
       iyear,  &
       imonth, &
       iday,   &
       idays   )
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: imonth
    integer, intent(out) :: iday
    integer, intent(in)  :: idays

    integer :: jy
    integer :: j4, j100, j400
    integer :: idays0
    integer :: ileap

    integer :: idy      ! # of days in a year
    integer :: m, id
    !---------------------------------------------------------------------------

    if ( operpt ) then ! return perpetual date
       iyear  = iyrpp
       imonth = imonpp
       iday   = idaypp
       return
    endif

    if ( oauto ) then ! switch
       if ( idays >= 693961 ) then ! 1900*365 + 1900/4-19+5
          ogrego = .true.
          oideal = .false.
       else
          ogrego = .false.
          if ( idays >= 1000*365 ) then
             oideal = .false.
          else
             oideal = .true.
             imonyr = 12
             idaymo = 30
          endif
       endif
    endif

    if ( ogrego ) then
       jy  = int( real(idays,kind=DP) / 365.24_DP )
       do
          j4     = (jy+3)   / 4
          j100   = (jy+99)  / 100
          j400   = (jy+399) / 400
          idays0 = jy * 365 + j4 - j100 + j400
          if ( idays <= idays0 ) then
             jy = jy -1
             if( jy < 0 ) exit
          endif
       enddo
       iyear = jy
       idy   = idays - idays0

       if ( checkleap(iyear) ) then
          ileap  = 2
       else
          ileap  = 1
       endif
    else
       if ( oideal ) then
          iyear = idays / (idaymo*imonyr)
          idy   = idays - iyear*(idaymo*imonyr)
          ileap = -1
       else
          iyear = idays / 365
          idy   = idays - iyear*365
          ileap = 1
       endif
    endif

    if ( oideal ) then
       imonth = idy / idaymo + 1
       iday   = idays - iyear*(idaymo*imonyr) - (imonth-1)*idaymo
    else
       id = 0
       do m = 1, 12
          id = id + monday(m,ileap)
          if ( idy <= id ) then
             imonth = m
             iday   = idy - id + monday(m,ileap)
             exit
          endif
       enddo
    endif

    return
  end subroutine CALENDAR_dd2ym

  !-----------------------------------------------------------------------------
  !> yymmdd -> day
  subroutine CALENDAR_ym2dd( &
       idays,  &
       iyear,  &
       imonth, &
       iday    )
    implicit none

    integer, intent(out) :: idays
    integer, intent(in)  :: iyear
    integer, intent(in)  :: imonth
    integer, intent(in)  :: iday

    integer :: jyear, jmonth
    integer :: j4, j100, j400
    integer :: idays0
    integer :: ileap

    integer :: m, id
    !---------------------------------------------------------------------------

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
          oideal = .false.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             imonyr = 12
             idaymo = 30
          endif
       endif
    endif

    if ( .NOT. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1) / 12
          jmonth = mod(imonth-1,12) + 1
       else
          jyear  = iyear - (-imonth) / 12 - 1
          jmonth = 12 - mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       j4     = (jyear+3)   / 4
       j100   = (jyear+99)  / 100
       j400   = (jyear+399) / 400
       idays0 = jyear*365 + j4 - j100 + j400

       if ( checkleap(jyear) ) then
          ileap = 2
       else
          ileap = 1
       endif
    else
       if ( oideal ) then
          idays0 = iyear * (idaymo*imonyr)
          ileap  = -1
       else
          idays0 = jyear * 365
          ileap  = 1
       endif
    endif

    if ( oideal ) then
       idays = idays0 + (imonth-1) * idaymo + iday
    else
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       enddo
       idays = idays0 + id + iday
    endif

    return
  end subroutine CALENDAR_ym2dd

  !-----------------------------------------------------------------------------
  !> seconds -> character
  subroutine CALENDAR_ss2cc( &
       htime, &
       dsec   )
    implicit none

    character(len=*), intent(out) :: htime
    real(DP),         intent(in)  :: dsec

    integer :: idate(6)
    integer :: i
    !---------------------------------------------------------------------------

    call CALENDAR_ss2yh( idate, & ! [OUT]
                         dsec   ) ! [IN]

    write(htime,'(I4.4,"/",I2.2,"/",I2.2,"-",I2.2,":",I2.2,":",I2.2 )') (idate(i),i=1,6)

    return
  end subroutine CALENDAR_ss2cc

  !-----------------------------------------------------------------------------
  !> Check leap year
  !> @return checkleap
  function checkleap( iyear )
    implicit none

    integer, intent(in) :: iyear     !< current year
    logical             :: checkleap

    integer :: check4, check100, check400
    !---------------------------------------------------------------------------

    check4   = mod(iyear,4  )
    check100 = mod(iyear,100)
    check400 = mod(iyear,400)

    checkleap = .false.
    if( check4   == 0 ) checkleap = .true.
    if( check100 == 0 ) checkleap = .false.
    if( check400 == 0 ) checkleap = .true.

  end function checkleap

end module mod_calendar
