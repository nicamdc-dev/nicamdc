!-------------------------------------------------------------------------------
!> module RANDOM
!!
!! @par Description
!!         Random number generation module
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_random
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
  public :: RANDOM_setup
  public :: RANDOM_get

  interface RANDOM_get
     module procedure RANDOM_get_SP
     module procedure RANDOM_get_DP
  end interface RANDOM_get

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: RANDOM_reset

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: RANDOM_FIX = .false.

  integer, private, allocatable :: RANDOM_seedvar(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine RANDOM_setup
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / RANDOMPARAM / &
       RANDOM_FIX

    integer  :: nseeds, ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[random]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=RANDOMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** RANDOMPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist RANDOMPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=RANDOMPARAM)

    call random_seed
    call random_seed(size=nseeds)

    allocate( RANDOM_seedvar(nseeds))

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Array size for random seed:', nseeds
    if ( RANDOM_FIX ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** random seed is fixed.'
    endif

    return
  end subroutine RANDOM_setup

  !-----------------------------------------------------------------------------
  !> Reset random seed
  subroutine RANDOM_reset
    use mod_process, only: &
       PRC_myrank
    implicit none

    integer  :: time1(8)
    real(DP) :: time2
    !---------------------------------------------------------------------------

    if ( RANDOM_FIX ) then
       time1(1) = 2011
       time1(2) = 12
       time1(3) = 5
       time1(4) = 10
       time1(5) = 20
       time1(6) = 41
       time2    = 0.0_RP
    else
       call date_and_time(values=time1)
       call cpu_time(time2)
    endif

    RANDOM_seedvar(:) = ( time1(1) - 1970 ) * 32140800 &
                      + time1(2) * 2678400 &
                      + time1(3) * 86400 &
                      + time1(4) * 60 &
                      + time1(5) * 3600 &
                      + time1(6) * 60 &
                      + time1(7) &
                      + int(time2*1.E6_DP) + PRC_myrank

    call random_seed(put=RANDOM_seedvar)

    return
  end subroutine RANDOM_reset

  !-----------------------------------------------------------------------------
  !> Get random number
  subroutine RANDOM_get_SP( var )
    implicit none

    real(SP), intent(out) :: var(:,:,:)
    !---------------------------------------------------------------------------

    call RANDOM_reset
    call random_number(var)

    return
  end subroutine RANDOM_get_SP

  !-----------------------------------------------------------------------------
  !> Get random number
  subroutine RANDOM_get_DP( var )
    implicit none

    real(DP), intent(out) :: var(:,:,:)
    !---------------------------------------------------------------------------

    call RANDOM_reset
    call random_number(var)

    return
  end subroutine RANDOM_get_DP

end module mod_random
