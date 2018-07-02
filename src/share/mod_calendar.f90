!-------------------------------------------------------------------------------
!> Module CALENDAR
!!
!! @par Description
!!         Calendar calculation
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
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: CALENDAR_setup
  public :: CALENDAR_ss2cc

  public :: CALENDAR_ss2yh
  public :: CALENDAR_yh2ss
  public :: CALENDAR_ss2ds
  public :: CALENDAR_ds2ss
  public :: CALENDAR_ym2dd
  public :: CALENDAR_dd2ym
  public :: CALENDAR_rs2hm
  public :: CALENDAR_hm2rs

  public :: CALENDAR_ym2yd ! 17/07/31 C.Kodama [add]
  public :: CALENDAR_ss2yd
  public :: CALENDAR_ss2ym
  public :: CALENDAR_xx2ss
  public :: CALENDAR_ointvl
  public :: CALENDAR_ssaft

  public :: CALENDAR_perpr
  public :: CALENDAR_daymo
  public :: CALENDAR_dayyr
  public :: CALENDAR_secdy

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  public :: CALENDAR_leapyr
  public :: CALENDAR_dgaus

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: oauto  = .true.  ! automatic calendar?
                                       ! year 0   - 999 : 360day
                                       ! year 1000-1899 : 365day
                                       ! year 1900-     : gregorian

  logical, private :: ogrego = .true.  ! use gregorian calendar?

  logical, private :: oideal = .false. ! ideal calender (n day per month)?

  integer, private :: imonyr = 12      ! 1 year   = x months in the ideal case
  integer, private :: idaymo = 30      ! 1 month  = x days   in the ideal case

  integer, private :: monday(12,2)     ! number of days in a month (normal year/leap year)

  data monday / 31,28,31,30,31,30,31,31,30,31,30,31, &
                31,29,31,30,31,30,31,31,30,31,30,31  /

  integer, private :: ihrday = 24      ! 1 day    = x hours
  integer, private :: iminhr = 60      ! 1 hour   = x minutes
  integer, private :: isecmn = 60      ! 1 minute = x sec.

  logical, private :: operpt = .false. ! perpetual?
  integer, private :: iyrpp  = 0       ! perpetual date(year)
  integer, private :: imonpp = 3       ! perpetual date(month)
  integer, private :: idaypp = 21      ! perpetual date(day)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine CALENDAR_setup
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist /CALENDARPARAM/ &
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

    logical :: dummy
    namelist /nm_calendar/ dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[calendar]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=CALENDARPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** CALENDARPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist CALENDARPARAM. STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist CALENDARPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=CALENDARPARAM)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF, nm_calendar, iostat=ierr)
    if ( ierr >= 0 ) then
       write(*         ,*) 'xxx namelist nm_calendar is changed to CALENDARPARAM. STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx namelist nm_calendar is changed to CALENDARPARAM. STOP.'
       call PRC_MPIstop
    endif

    return
  end subroutine CALENDAR_setup

  !-----------------------------------------------------------------------------
  !> second to character
  subroutine CALENDAR_ss2cc( &
       htime,      &
       dsec,       &
       number_only )
    implicit none

    character(len=*), intent(out) :: htime
    real(DP),         intent(in)  :: dsec
    logical,          intent(in), optional :: number_only

    integer :: itime(6)
    !---------------------------------------------------------------------------

    call CALENDAR_ss2yh( itime, dsec )

    if ( present(number_only) ) then
       if ( number_only ) then
          write(htime,'(I4.4,5(I2.2))') itime(1), itime(2), itime(3), &
                                        itime(4), itime(5), itime(6)
          return
       endif
    endif

    write(htime,'(I4.4,5(A,I2.2))') itime(1), '/', itime(2), '/', itime(3), &
                               '-', itime(4), ':', itime(5), ':', itime(6)

    return
  end subroutine CALENDAR_ss2cc

  !-----------------------------------------------------------------------------
  !> second -> date
  subroutine CALENDAR_ss2yh( &
       idate, &
       dsec   )
    implicit none

    integer,  intent(out) :: idate(6) ! yymmddhhmmss
    real(DP), intent(in)  :: dsec     ! time(second)

    integer  :: idays ! serial number of day
    real(DP) :: rsec  ! seconds in a day
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
  !> date -> second
  subroutine CALENDAR_yh2ss( &
       dsec, &
       idate )
    implicit none

    real(DP), intent(out) :: dsec     ! time(second)
    integer,  intent(in)  :: idate(6) ! yymmddhhmmss

    integer  :: idays ! serial number of day
    real(DP) :: rsec  ! seconds in a day
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
  !> seconds -> day,seconds
  subroutine CALENDAR_ss2ds( &
       idays, &
       rsec,  &
       dsec   )
    implicit none

    integer,  intent(out) :: idays
    real(DP), intent(out) :: rsec
    real(DP), intent(in)  :: dsec

    integer  :: isecdy
    !---------------------------------------------------------------------------

    isecdy = isecmn * iminhr * ihrday
    idays  = int( dsec / real(isecdy,kind=DP) ) + 1
    rsec   = dsec - real(idays-1,kind=DP) * real(isecdy,kind=DP)

    if ( nint(rsec) >= isecdy ) then
       idays = idays + 1
       rsec  = rsec - real(isecdy,kind=DP)
    endif

    return
  end subroutine CALENDAR_ss2ds

  !-----------------------------------------------------------------------------
  !> day,seconds -> seconds
  subroutine CALENDAR_ds2ss( &
       dsec,  &
       idays, &
       rsec   )
    implicit none

    real(DP), intent(out) :: dsec
    integer,  intent(in)  :: idays
    real(DP), intent(in)  :: rsec

    integer  :: isecdy
    !---------------------------------------------------------------------------

    isecdy = isecmn * iminhr * ihrday
    dsec   = rsec + real(idays-1,kind=DP) * real(isecdy,kind=DP)

    return
  end subroutine CALENDAR_ds2ss

  !-----------------------------------------------------------------------------
  !> day -> year,month,day
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

    integer  :: jyear, jy4, jcent, jcent4, ileap
    integer  :: jdays, id, idayyr
    integer  :: m
    !---------------------------------------------------------------------------

    if ( oauto ) then
       if ( idays >= 693961 ) then ! 1900*365 + 1900/4 - 19 + 5
          ogrego = .true.
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
       jyear = int( real(idays,kind=DP) / 365.24_DP ) ! first guess
       do
          jy4    = ( jyear + 3   ) / 4
          jcent  = ( jyear + 99  ) / 100
          jcent4 = ( jyear + 399 ) / 400

          id = jyear * 365 + jy4 - jcent + jcent4
          if ( idays <= id ) then
             jyear = jyear - 1
          else
             exit
          endif
       enddo
       iyear = jyear
       jdays = idays - id

       if ( CALENDAR_leapyr(iyear) ) then
          ileap  = 2
       else
          ileap  = 1
       endif
    elseif( .NOT. oideal ) then
       iyear = idays / 365
       jdays = idays - iyear*365
       ileap = 1
    endif

    if ( ogrego .OR. .NOT. oideal ) then
       id = 0
       do m = 1, 12
          id = id + monday(m,ileap)
          if ( jdays <= id ) then
             imonth = m
             iday   = jdays + monday(m,ileap) - id
             exit
          endif
       enddo
    else
       idayyr = idaymo * imonyr
       iyear  = ( idays - 1 ) / idayyr
       imonth = ( idays - iyear*idayyr - 1 ) / idaymo + 1
       iday   = idays - iyear*idayyr - (imonth-1)*idaymo
    endif

    return
  end subroutine CALENDAR_dd2ym

  !-----------------------------------------------------------------------------
  !> year,month,day -> day
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

    integer  :: jyear, jy4, jcent, jcent4, ileap
    integer  :: jmonth, m
    integer  :: idayyr
    !---------------------------------------------------------------------------

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
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

    if ( ogrego .OR. .NOT. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1) / 12
          jmonth = mod(imonth-1,12) + 1
       else
          jyear  = iyear - (-imonth) / 12 - 1
          jmonth = 12 - mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       jy4    = ( jyear + 3   ) / 4
       jcent  = ( jyear + 99  ) / 100
       jcent4 = ( jyear + 399 ) / 400

       idays = iday + jyear * 365 + jy4 - jcent + jcent4

       if ( CALENDAR_leapyr(jyear) ) then
          ileap  = 2
       else
          ileap  = 1
       endif
    elseif( .NOT. oideal ) then
       idays = iday + jyear * 365
       ileap = 1
    endif

    if ( ogrego .OR. .NOT. oideal ) then
       do m = 1, jmonth-1
          idays = idays + monday(m,ileap)
       enddo
    else
       idayyr = idaymo * imonyr
       idays  = iday
       idays  = idays + iyear * idayyr
       idays  = idays + (imonth-1) * idaymo
    endif

    return
  end subroutine CALENDAR_ym2dd

  !-----------------------------------------------------------------------------
  !> second -> hour,minute,second
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

    real(DP) :: isecs
    !---------------------------------------------------------------------------

    ihour = int( rsec / real(iminhr*isecmn,kind=DP) )
    isecs = rsec - real(ihour*iminhr*isecmn,kind=DP)

    imin  = int( isecs / real(isecmn,kind=DP) )
    isecs = isecs - real(imin*isecmn,kind=DP)

    isec  = nint( isecs )

    if ( isec >= isecmn ) then
       imin = imin + 1
       isec = isec - isecmn
    endif

    if ( imin == iminhr ) then
       ihour = ihour + 1
       imin  = imin  - iminhr
    endif

    return
  end subroutine CALENDAR_rs2hm

  !-----------------------------------------------------------------------------
  !> hour,minute,second -> second
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

    integer :: isecs
    !---------------------------------------------------------------------------

    isecs = ihour*iminhr*isecmn + imin*isecmn + isec
    rsec  = real(isecs,kind=DP)

    return
  end subroutine CALENDAR_hm2rs

  !-----------------------------------------------------------------------------
  !> year,month,day -> day of year
  subroutine CALENDAR_ym2yd( &
       idays,  &
       iyear,  &
       imonth, &
       iday    )
    implicit none

    integer, intent(out) :: idays
    integer, intent(in)  :: iyear
    integer, intent(in)  :: imonth
    integer, intent(in)  :: iday

    integer  :: jyear, ileap
    integer  :: jmonth, m
    !---------------------------------------------------------------------------

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
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

    if ( ogrego .OR. .NOT. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1) / 12
          jmonth = mod(imonth-1,12) + 1
       else
          jyear  = iyear - (-imonth) / 12 - 1
          jmonth = 12 - mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       if ( CALENDAR_leapyr(jyear) ) then
          ileap = 2
       else
          ileap = 1
       endif
    elseif( .NOT. oideal ) then
       ileap = 1
    endif

    if ( ogrego .OR. .NOT. oideal ) then
       idays = iday
       do m = 1, jmonth-1
          idays = idays + monday(m,ileap)
       enddo
    else
       idays = iday + (imonth-1) * idaymo
    endif

    return
  end subroutine CALENDAR_ym2yd

  !-----------------------------------------------------------------------------
  !> second -> year, day of year
  subroutine CALENDAR_ss2yd( &
       iyear, &
       idays, &
       dsec   )
    implicit none

    integer,  intent(out) :: iyear
    integer,  intent(out) :: idays
    real(DP), intent(in)  :: dsec

    integer :: imonth, iday
    !---------------------------------------------------------------------------

    call CALENDAR_ss2ym( iyear, imonth, iday, & ! [OUT]
                         dsec                 ) ! [IN]

    call CALENDAR_ym2yd( idays,               & ! [OUT]
                         iyear, imonth, iday  ) ! [IN]

    return
  end subroutine CALENDAR_ss2yd

  !-----------------------------------------------------------------------------
  !> second -> year, month, day
  subroutine CALENDAR_ss2ym( &
       iyear,  &
       imonth, &
       iday,   &
       dsec    )
    implicit none

    integer,  intent(out) :: iyear
    integer,  intent(out) :: imonth
    integer,  intent(out) :: iday
    real(DP), intent(in)  :: dsec

    integer  :: idays ! serial number of day
    real(DP) :: rsec  ! seconds in a day
    !---------------------------------------------------------------------------

    call CALENDAR_ss2ds( idays, rsec, & ! [OUT]
                         dsec         ) ! [IN]

    call CALENDAR_dd2ym( iyear, imonth, iday, & ! [OUT]
                         idays                ) ! [IN]

    return
  end subroutine CALENDAR_ss2ym

  !-----------------------------------------------------------------------------
  subroutine CALENDAR_xx2ss( &
       ddsec, &
       rtdur, &
       hunit, &
       dsec   )
    implicit none

    real(DP),         intent(out) :: ddsec
    real(DP),         intent(in)  :: rtdur
    character(len=*), intent(in)  :: hunit
    real(DP),         intent(in)  :: dsec

    character(len=10) :: hunitx

    integer :: iyear, imonth, iday, ndaymo, ndayyr
    !---------------------------------------------------------------------------

    hunitx = trim(hunit)

    if    ( hunitx(1:1) == 's'  .OR. hunitx(1:1) == 'S'  ) then

       ddsec = rtdur

    elseif( hunitx(1:2) == 'mi' .OR. hunitx(1:2) == 'MI' ) then

       ddsec = rtdur * real(isecmn,kind=DP)

    elseif( hunitx(1:1) == 'h'  .OR. hunitx(1:1) == 'H'  ) then

       ddsec = rtdur * real(iminhr*isecmn,kind=DP)

    elseif( hunitx(1:1) == 'd' .OR. hunitx(1:1) == 'D' ) then

       ddsec = rtdur * real(ihrday*iminhr*isecmn,kind=DP)

    elseif( hunitx(1:2) == 'mo' .OR. hunitx(1:2) == 'MO' ) then

       call CALENDAR_ss2ym( iyear, imonth, iday, & ! [OUT]
                            dsec                 ) ! [IN]

       call CALENDAR_daymo( ndaymo,       & ! [OUT]
                            iyear, imonth ) ! [IN]

       ddsec = rtdur * real(ndaymo,kind=DP) * real(ihrday*iminhr*isecmn,kind=DP)

    elseif( hunitx(1:1) == 'y'  .OR. hunitx(1:1) == 'Y'  ) then

       call CALENDAR_ss2ym( iyear, imonth, iday, & ! [OUT]
                            dsec                 ) ! [IN]

       call CALENDAR_dayyr( ndayyr, & ! [OUT]
                            iyear   ) ! [IN]

       ddsec = rtdur * real(ndayyr,kind=DP) * real(ihrday*iminhr*isecmn,kind=DP)

    else

       write(*,*) 'xxx [CALENDAR_xx2ss] unknown unit: ', hunit, '. Assumed as second.'
       ddsec = rtdur

    endif

    return
  end subroutine CALENDAR_xx2ss

  !-----------------------------------------------------------------------------
  subroutine CALENDAR_ssaft( &
       dsec_after, &
       dsec,       &
       raftr,      &
       hunit       )
    implicit none

    real(DP),         intent(out) :: dsec_after
    real(DP),         intent(in)  :: dsec
    real(DP),         intent(in)  :: raftr
    character(len=*), intent(in)  :: hunit

    real(DP) :: dt
    !---------------------------------------------------------------------------

    call CALENDAR_xx2ss( dt, raftr, hunit, dsec )
    dsec_after = dsec + dt

    return
  end subroutine CALENDAR_ssaft

  !-----------------------------------------------------------------------------
  !> time step passed ?
  function CALENDAR_ointvl( &
       dtime,  &
       dtprev, &
       dtorgn, &
       rintv,  &
       htunit  )
    implicit none

    real(DP),         intent(in) :: dtime
    real(DP),         intent(in) :: dtprev
    real(DP),         intent(in) :: dtorgn
    real(DP),         intent(in) :: rintv
    character(len=*), intent(in) :: htunit
    logical                      :: CALENDAR_ointvl

    character(len=5) :: hunit
    integer          :: iyear,  imon,  iday
    integer          :: iyearp, imonp, idayp
    real(DP)         :: dt

    real(DP)         :: ry, rmo
    integer          :: ndayyr, ndaymo
    !---------------------------------------------------------------------------

    hunit = trim(htunit)

    if ( dtime == dtprev ) then
       CALENDAR_ointvl = .true.
       return
    endif

    CALENDAR_ointvl = .false.

    call CALENDAR_ss2ym( iyear,  imon,  iday,  dtime  )
    call CALENDAR_ss2ym( iyearp, imonp, idayp, dtprev )
    call CALENDAR_xx2ss( dt, rintv, hunit, dtime )

    if ( dtime >= dtorgn ) then
       if      ( hunit(1:1) == 'y' .OR. hunit(1:1) == 'Y' ) then

          call CALENDAR_dayyr( ndayyr, iyear )

          ry = real( iyear-iyearp, kind=DP )                           &
             + real( imon -imonp,  kind=DP ) / real( imonyr, kind=DP ) &
             + real( iday -idayp,  kind=DP ) / real( ndayyr, kind=DP )

          if ( ry >= rintv ) then
             CALENDAR_ointvl = .true.
          endif

       elseif( hunit(1:2) == 'mo' .OR. hunit(1:2) == 'MO' ) then

          call CALENDAR_daymo( ndaymo, iyear, imon )

          rmo = real( iyear-iyearp-1, kind=DP ) * real( imonyr, kind=DP ) &
              + real( imon -imonp,    kind=DP )                           &
              + real( iday -idayp,    kind=DP ) / real( ndaymo, kind=DP )

          if ( rmo >= rintv ) then
             CALENDAR_ointvl = .true.
          endif

       elseif( CALENDAR_dgaus((dtime-dtorgn)/dt) > CALENDAR_dgaus((dtprev-dtorgn)/dt) ) then

          CALENDAR_ointvl = .true.

       endif
    endif

    return
  end function CALENDAR_ointvl

  !-----------------------------------------------------------------------------
  !> return perpetual date
  subroutine CALENDAR_perpr( &
       iyear , &
       imonth, &
       iday,   &
       ooperp  )
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: imonth
    integer, intent(out) :: iday
    logical, intent(out) :: ooperp
    !---------------------------------------------------------------------------

    ooperp = operpt
    iyear  = iyrpp
    imonth = imonpp
    iday   = idaypp

    return
  end subroutine CALENDAR_perpr

  !-----------------------------------------------------------------------------
  !> return # of days in a month
  subroutine CALENDAR_daymo(&
       ndaymo, &
       iyear,  &
       imonth  )
    implicit none

    integer, intent(out) :: ndaymo
    integer, intent(in)  :: iyear
    integer, intent(in)  :: imonth
    !---------------------------------------------------------------------------

    if ( oauto ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego ) then
       if ( CALENDAR_leapyr(iyear) ) then
          ndaymo = monday(imonth,2)
       else
          ndaymo = monday(imonth,1)
       endif
    elseif( .NOT. oideal ) then
       ndaymo = monday(imonth,1)
    else
       ndaymo = idaymo
    endif

    return
  end subroutine CALENDAR_daymo

  !-----------------------------------------------------------------------------
  !> return # of days in an year
  subroutine CALENDAR_dayyr(&
       ndayyr, &
       iyear   )
    implicit none

    integer, intent(out) :: ndayyr
    integer, intent(in)  :: iyear
    !---------------------------------------------------------------------------

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego ) then
       if ( CALENDAR_leapyr(iyear) ) then
          ndayyr = 366
       else
          ndayyr = 365
       endif
    elseif( .NOT. oideal ) then
       ndayyr = 365
    else
       ndayyr = idaymo * imonyr
    endif

    return
  end subroutine CALENDAR_dayyr

  !-----------------------------------------------------------------------------
  !> return # of seconds in a day
  subroutine CALENDAR_secdy( &
       nsecdy )
    implicit none

    integer, intent(out) :: nsecdy
    !---------------------------------------------------------------------------

    nsecdy = ihrday * iminhr * isecmn

    return
  end subroutine CALENDAR_secdy

  !-----------------------------------------------------------------------------
  !> check leap year
  function CALENDAR_leapyr( &
       iyear )
    implicit none

    integer, intent(in) :: iyear
    logical             :: CALENDAR_leapyr

    integer :: iy, iycen, icent
    !---------------------------------------------------------------------------

    iy     = mod(iyear,4)
    iycen  = mod(iyear,100)
    icent  = mod(iyear/100,4)

    if ( iy == 0 .AND. ( iycen /= 0 .OR. icent == 0 ) ) then
       CALENDAR_leapyr = .true.
    else
       CALENDAR_leapyr = .false.
    endif

    return
  end function CALENDAR_leapyr

  !-----------------------------------------------------------------------------
  !> dicard gaussian
  function CALENDAR_dgaus( &
       dx                )
    implicit none

    real(DP), intent(in) :: dx
    real(DP)             :: CALENDAR_dgaus
    !---------------------------------------------------------------------------

    CALENDAR_dgaus = aint(dx) + aint( dx - aint(dx) + 1.0_DP ) - 1.0_DP

  end function CALENDAR_dgaus

end module mod_calendar
