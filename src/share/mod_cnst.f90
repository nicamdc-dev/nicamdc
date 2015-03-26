!-------------------------------------------------------------------------------
!>
!! Constant definition module
!!
!! @par Description
!!         This module defines constant numbers.
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita) Imported from igdc-4.33
!! @li      2009-09-28 (      )   add CNST_undef2 for 2 byte integer
!!
!<
module mod_cnst
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
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
  public :: CNST_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  REAL(RP), public :: CNST_ERADIUS = 6.37122D+6 ! Radius of the Earth [m]
  REAL(RP), public :: CNST_EOHM    = 7.292D-5   ! Angular velocity of the Earth [/s]
  REAL(RP), public :: CNST_EGRAV   = 9.80616D0  ! Gravitational accerlaration of the Earth [m/s2]

  REAL(RP), public :: CNST_RAIR    =  287.0D0   ! Gas constant of air
  REAL(RP), public :: CNST_RVAP    =  461.5D0   ! Gas constant of vapor

  REAL(RP), public :: CNST_CP      = 1004.5D0   ! Specific heat of air (consant pressure)
  REAL(RP), public :: CNST_CV                   ! Specific heat of air (consant volume)

  REAL(RP), public :: CNST_CPV     = 1846.0D0   ! Specific heat of vapor (consant pressure)
  REAL(RP), public :: CNST_CVV                  ! Specific heat of vapor (consant volume)
  REAL(RP), public :: CNST_CL      = 4218.0D0   ! Specific heat of water
  REAL(RP), public :: CNST_CI      = 2006.0D0   ! Specific heat of ice

  REAL(RP), public :: CNST_SOUND                ! Speed of sound at 0C, dry

  !------ cp/cv
  REAL(RP), public :: CNST_GAMMA
  !<----- calculated in sub[CNST_setup].
  !
  !------ R/cp
  REAL(RP), public :: CNST_KAPPA
  !<----- calculated in sub[CNST_setup].
  !
  !------ dry lapse rate [K/m]
  REAL(RP), public :: CNST_LAPS
  !<----- calculated in sub[CNST_setup].
  !
  !------ molecular weight ( water/air )
  REAL(RP), public :: CNST_EPSV
  !<----- calculated in sub[CNST_setup].
  !
  !------ 1/epsv-1
  REAL(RP), public :: CNST_EPSVT
  !<----- calculated in sub[CNST_setup].
  !
  !------ Density of water
  REAL(RP), public :: CNST_DWATR = 1000.D0
  !
  !------ Saturate pressure of water vapor at 0C
  REAL(RP), public :: CNST_PSAT0 = 610.7D0
  !<----- unit : [Pa]
  !
  !------ Latent heat of vaporizaion at 0C
!  REAL(RP), public :: CNST_LH0   = 2.5008D+6 [mod] 20120704 H.Yashiro
  REAL(RP), public :: CNST_LH0   = 2.501D+6
  !
  !------ Latent heat of vaporizaion at 0K
  REAL(RP), public :: CNST_LH00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of sublimation at 0C
!  REAL(RP), public :: CNST_LHS0  = 2.8342D+6 [mod] 20120704 H.Yashiro
  REAL(RP), public :: CNST_LHS0  = 2.834D+6
  !
  !------ Latent heat of sublimation at 0K
  REAL(RP), public :: CNST_LHS00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of fusion at 0C
  REAL(RP), public :: CNST_LHF0
  !<----- calculated in sub[CNST_setup].
  !
  !------ latent heat of fusion at 0K
  REAL(RP), public :: CNST_LHF00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of melting
  REAL(RP), public :: CNST_EMELT = 3.40D+5
  !
  !------ Melting temperature of water
  REAL(RP), public :: CNST_TMELT = 273.15D0
  !
  !------ Freeze point of sea
  REAL(RP), public :: CNST_TFRZS  = 271.35D0
  !
  !------ Wet-bulb temp. rain/snow
  REAL(RP), public :: CNST_TQICE = 273.15D0
  !
  !------ Stefan-Boltzman constant
  REAL(RP), public :: CNST_STB   = 5.67D-8
  !
  !------ Karman constant
  REAL(RP), public :: CNST_KARMAN = 0.4D0
  !
  !------ Surface pressure
  REAL(RP), public :: CNST_PRES0    = 101325.0D0
  !
  !------ Surface temperature
  REAL(RP), public :: CNST_TEMS0    = 300.0D0
  !
  !------ Standard pressure
  REAL(RP), public :: CNST_PRE00    = 1.0D+5
  !
  !------ Standard temperature
  REAL(RP), public :: CNST_TEM00    = 273.15D0
  !
  !------ Standard density
  REAL(RP), public :: CNST_RHO00
  !<----- calculated in sub[CNST_setup].
  !
  !====== Misc. constants ======
  !
  !------ Definition of PI
  REAL(RP), public :: CNST_PI = 3.14159265358979323846D0

  REAL(RP), public :: CNST_D2R

  !------ Allowable minimum value
  REAL(RP), public, parameter :: CNST_EPS_ZERO = 1.0D-99
  !
  !------ Allowable maximum value
  REAL(RP), public, parameter :: CNST_MAX_REAL = 1.0D+99
  !
  !------ Missing value
  REAL(RP), public, parameter :: CNST_VMISS    = 0.0D0
  !
  !------ Undefined value
  real(RP), public :: CNST_UNDEF
  !
  !------ Undefined value
  REAL(RP), public, parameter :: CNST_UNDEF8   = -99.9D+33
  !
  !------ Undefined value
  real(4), public, parameter :: CNST_UNDEF4   = -99.9E+33
  !
  !------ Undefined value
  integer(4), public, parameter :: CNST_UNDEF2   = -32768

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine CNST_setup
  !>
  subroutine CNST_setup
    use mod_adm, only: &
       ADM_LOG_FID, &
       ADM_CTL_FID, &
       ADM_proc_stop
    implicit none

    REAL(RP) :: earth_radius               ! Earth radius
    REAL(RP) :: earth_angvel               ! Anguler velocity of the earth
    REAL(RP) :: small_planet_factor = 1.D0 ! small planet factor
    REAL(RP) :: earth_gravity              ! Gravitational accelaration
    REAL(RP) :: gas_cnst                   ! Gas constant of dry air
    REAL(RP) :: gas_cnst_vap               ! Gas constant of water vapour
    REAL(RP) :: specific_heat_pre          ! Specific heat of air( const pre )
    REAL(RP) :: specific_heat_pre_vap      ! Specific heat of water vapour ( const pre )
    REAL(RP) :: latent_heat_vap            ! latent heat of vaporization LH0 ( 0 deg )
    REAL(RP) :: latent_heat_sub            ! latent heat of sublimation LHS0 ( 0 deg )

    namelist / CNSTPARAM / &
       earth_radius,          &
       earth_angvel,          &
       small_planet_factor,   &
       earth_gravity,         &
       gas_cnst,              &
       gas_cnst_vap,          &
       specific_heat_pre,     &
       specific_heat_pre_vap, &
       latent_heat_vap,       &
       latent_heat_sub

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- initialization of controled parameters
    earth_radius          = CNST_ERADIUS
    earth_angvel          = CNST_EOHM
    earth_gravity         = CNST_EGRAV
    gas_cnst              = CNST_RAIR
    gas_cnst_vap          = CNST_RVAP
    specific_heat_pre     = CNST_CP
    specific_heat_pre_vap = CNST_CPV
    latent_heat_vap       = CNST_LH0
    latent_heat_sub       = CNST_LHS0

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[cnst]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=CNSTPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** CNSTPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist CNSTPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist CNSTPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=CNSTPARAM)

    CNST_ERADIUS = earth_radius / small_planet_factor
    CNST_EOHM    = earth_angvel * small_planet_factor
    CNST_EGRAV   = earth_gravity
    CNST_RAIR    = gas_cnst
    CNST_RVAP    = gas_cnst_vap
    CNST_CP      = specific_heat_pre
    CNST_CPV     = specific_heat_pre_vap
    CNST_LH0     = latent_heat_vap
    CNST_LHS0    = latent_heat_sub

    !--- calculate other parameters
    if    ( RP == SP ) then
       CNST_UNDEF = real(CNST_UNDEF4,kind=RP)
    elseif( RP == DP ) then
       CNST_UNDEF = real(CNST_UNDEF8,kind=RP)
    else
       write(*,*) 'xxx unsupported precision: ', RP
       call ADM_proc_stop
    endif

    CNST_PI    = 4.D0 * atan( 1.D0 )
    CNST_D2R   = CNST_PI / 180.D0

    CNST_CV    = CNST_CP - CNST_RAIR
    CNST_SOUND = sqrt( CNST_CP * CNST_RAIR / ( CNST_CP - CNST_RAIR ) * CNST_TEM00 )
    CNST_GAMMA = CNST_CP / CNST_CV
    CNST_KAPPA = CNST_RAIR / CNST_CP
    CNST_RHO00 = CNST_PRE00 / CNST_RAIR / CNST_TEM00
    CNST_LAPS  = CNST_EGRAV / CNST_CP

    CNST_CVV   = CNST_CPV - CNST_RVAP
    CNST_EPSV  = CNST_RAIR / CNST_RVAP
    CNST_EPSVT = 1.0D0/CNST_EPSV - 1.0D0

    CNST_LH00  = CNST_LH0  - ( CNST_CPV - CNST_CL ) * CNST_TEM00
    CNST_LHS00 = CNST_LHS0 - ( CNST_CPV - CNST_CI ) * CNST_TEM00
    CNST_LHF0  = CNST_LHS0 - CNST_LH0
    CNST_LHF00 = CNST_LHF0 - ( CNST_CL - CNST_CI ) * CNST_TEM00 ! bugfix: CNST_LHS0 -> CNST_LHF0

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '====== Physical Constants ======'
    write(ADM_LOG_FID,*) '--- The circluar constant PI                 : ', CNST_PI
    write(ADM_LOG_FID,*) '--- Allowable minimum value                  : ', CNST_EPS_ZERO
    write(ADM_LOG_FID,*) '--- Allowable maximum value                  : ', CNST_MAX_REAL
    write(ADM_LOG_FID,*) '--- Missing value                            : ', CNST_VMISS
    write(ADM_LOG_FID,*) '--- Radius of the Earth                      : ', CNST_ERADIUS
    write(ADM_LOG_FID,*) '--- Angular velocity of the Earth            : ', CNST_EOHM
    write(ADM_LOG_FID,*) '--- Gravitational accerlaration of the Earth : ', CNST_EGRAV
    write(ADM_LOG_FID,*) '--- Gas constant of air                      : ', CNST_RAIR
    write(ADM_LOG_FID,*) '--- Gas constant of vapor                    : ', CNST_RVAP
    write(ADM_LOG_FID,*) '--- Specific heat of air (consant pressure)  : ', CNST_CP
    write(ADM_LOG_FID,*) '--- Specific heat of air (consant volume)    : ', CNST_CV
    write(ADM_LOG_FID,*) '--- Speed of sound at 0C                     : ', CNST_SOUND
    write(ADM_LOG_FID,*) '--- Rair/Cp                                  : ', CNST_KAPPA
    write(ADM_LOG_FID,*) '--- Surface pressure                         : ', CNST_PRES0
    write(ADM_LOG_FID,*) '--- Standard pressure                        : ', CNST_PRE00
    write(ADM_LOG_FID,*) '--- Standard temperature                     : ', CNST_TEM00
    write(ADM_LOG_FID,*) '--- Standard density                         : ', CNST_RHO00
    write(ADM_LOG_FID,*) '--- Specific heat of vapor (consant pressure): ', CNST_CPV
    write(ADM_LOG_FID,*) '--- Specific heat of vapor (consant volume)  : ', CNST_CVV
    write(ADM_LOG_FID,*) '--- Specific heat of water                   : ', CNST_CL
    write(ADM_LOG_FID,*) '--- Specific heat of ice                     : ', CNST_CI
    write(ADM_LOG_FID,*) '--- Mocular weight                           : ', CNST_EPSV
    write(ADM_LOG_FID,*) '--- 1/epsv-1                                 : ', CNST_EPSVT
    write(ADM_LOG_FID,*) '--- Density of water                         : ', CNST_DWATR
    write(ADM_LOG_FID,*) '--- Saturate pressure of water vapor at 0C   : ', CNST_PSAT0
    write(ADM_LOG_FID,*) '--- Latent heat of vaporizaion at 0C         : ', CNST_LH0
    write(ADM_LOG_FID,*) '--- Latent heat of vaporizaion at 0K         : ', CNST_LH00
    write(ADM_LOG_FID,*) '--- Latent heat of sublimation at 0C         : ', CNST_LHS0
    write(ADM_LOG_FID,*) '--- Latent heat of sublimation at 0K         : ', CNST_LHS00
    write(ADM_LOG_FID,*) '--- Latent heat of fusion at 0C              : ', CNST_LHF0
    write(ADM_LOG_FID,*) '--- Latent heat of fusion at 0K              : ', CNST_LHF00
    write(ADM_LOG_FID,*) '--- Latent heat of melt                      : ', CNST_EMELT
    write(ADM_LOG_FID,*) '--- Melting temperature of water             : ', CNST_TMELT
    write(ADM_LOG_FID,*) '--- Freeze point of sea                      : ', CNST_TFRZS
    write(ADM_LOG_FID,*) '--- Wet-bulb temperature rain/snow           : ', CNST_TQICE
    write(ADM_LOG_FID,*) '--- Stefan-Boltzman constant                 : ', CNST_STB
    write(ADM_LOG_FID,*) '--- Karman constant                          : ', CNST_KARMAN
    write(ADM_LOG_FID,*) '--- Cp/Cv                                    : ', CNST_GAMMA

    return
  end subroutine CNST_setup

end module mod_cnst
!-------------------------------------------------------------------------------
