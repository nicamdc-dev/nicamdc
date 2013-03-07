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
  !======  Earth parameters ======
  !
  !------ Radius of the Earth
!  real(8), public, save :: CNST_ERADIUS  = 6.37122D+6 [mod] 20120704 H.Yashiro
  real(8), public, save :: CNST_ERADIUS = 6371000.D0
  !<----- unit : [m]
  !
  !------ Angular velocity of the Earth
!  real(8), public, save :: CNST_EOHM     = 7.292D-5 [mod] 20120704 H.Yashiro
  real(8), public, save :: CNST_EOHM = 7.292115D-5
  !<----- unit : [/s]
  !
  !------ Gravitational accerlaration of the Earth
!  real(8), public, save :: CNST_EGRAV    = 9.80616D0 [mod] 20120704 H.Yashiro
  real(8), public, save :: CNST_EGRAV    = 9.79764D0
  !<----- unit : [m/s^2]
  !
  !====== Gas parameters ======
  !
  !------ Gas constant of air
  real(8), public, save :: CNST_RAIR     = 287.04D0
  !
  !------ Gas constant of vapor
  real(8), public, save :: CNST_RVAP     = 461.50D0
  !
  !------ Specific heat of air (consant pressure)
  real(8), public, save :: CNST_CP       = 1004.6D0
  !
  !------ Specific heat of vapor (consant pressure)
!  real(8), public, save :: CNST_CPV      = 1850.0D0 [mod] 20120704 H.Yashiro
  real(8), public, save :: CNST_CPV      = 1846.0D0
  !
  !------ Specific heat of water
  real(8), public, save :: CNST_CL   = 4218.0D0
  !
  !------ Specific heat of ice
  real(8), public, save :: CNST_CI   = 2006.0D0
  !
  !------ Specific heat of air (consant volume)
  real(8), public, save :: CNST_CV
  !<----- calculated in sub[CNST_setup].
  !
  !------ Specific heat of vapor (consant volume)
  real(8), public, save :: CNST_CVV
  !<----- calculated in sub[CNST_setup].
  !
  !------ cp/cv
  real(8), public, save :: CNST_GAMMA
  !<----- calculated in sub[CNST_setup].
  !
  !------ R/cp
  real(8), public, save :: CNST_KAPPA
  !<----- calculated in sub[CNST_setup].
  !
  !------ dry lapse rate [K/m]
  real(8), public, save :: CNST_LAPS
  !<----- calculated in sub[CNST_setup].
  !
  !------ molecular weight ( water/air )
  real(8), public, save :: CNST_EPSV
  !<----- calculated in sub[CNST_setup].
  !
  !------ 1/epsv-1
  real(8), public, save :: CNST_EPSVT
  !<----- calculated in sub[CNST_setup].
  !
  !------ Density of water
  real(8), public, save :: CNST_DWATR = 1000.D0
  !
  !------ Saturate pressure of water vapor at 0C
  real(8), public, save :: CNST_PSAT0 = 610.7D0
  !<----- unit : [Pa]
  !
  !------ Latent heat of vaporizaion at 0C
!  real(8), public, save :: CNST_LH0   = 2.5008D+6 [mod] 20120704 H.Yashiro
  real(8), public, save :: CNST_LH0   = 2.501D+6
  !
  !------ Latent heat of vaporizaion at 0K
  real(8), public, save :: CNST_LH00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of sublimation at 0C
!  real(8), public, save :: CNST_LHS0  = 2.8342D+6 [mod] 20120704 H.Yashiro
  real(8), public, save :: CNST_LHS0  = 2.834D+6
  !
  !------ Latent heat of sublimation at 0K
  real(8), public, save :: CNST_LHS00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of fusion at 0C
  real(8), public, save :: CNST_LHF0
  !<----- calculated in sub[CNST_setup].
  !
  !------ latent heat of fusion at 0K
  real(8), public, save :: CNST_LHF00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of melting
  real(8), public, save :: CNST_EMELT = 3.40D+5
  !
  !------ Melting temperature of water
  real(8), public, save :: CNST_TMELT = 273.15D0
  !
  !------ Freeze point of sea
  real(8), public, save :: CNST_TFRZS  = 271.35D0
  !
  !------ Wet-bulb temp. rain/snow
  real(8), public, save :: CNST_TQICE = 273.15D0
  !
  !------ Stefan-Boltzman constant
  real(8), public, save :: CNST_STB   = 5.67D-8
  !
  !------ Karman constant
  real(8), public, save :: CNST_KARMAN = 0.4D0
  !
  !------ Surface pressure
  real(8), public, save :: CNST_PRES0    = 101325.0D0
  !
  !------ Surface temperature
  real(8), public, save :: CNST_TEMS0    = 300.0D0
  !
  !------ Standard pressure
  real(8), public, save :: CNST_PRE00    = 1.0D+5
  !
  !------ Standard temperature
  real(8), public, save :: CNST_TEM00    = 273.15D0
  !
  !------ Standard density
  real(8), public, save :: CNST_RHO00
  !<----- calculated in sub[CNST_setup].
  !
  !====== Misc. constants ======
  !
  !------ Definition of PI
  real(8), public, save :: CNST_PI = 3.14159265358979323846D0
  !
  !------ Allowable minimum value
  real(8), public, parameter :: CNST_EPS_ZERO = 1.0D-99
  !
  !------ Allowable maximum value
  real(8), public, parameter :: CNST_MAX_REAL = 1.0D+99
  !
  !------ Missing value
  real(8), public, parameter :: CNST_VMISS    = 0.0D0
  !
  !------ Undefined value
  real(8), public, parameter :: CNST_UNDEF    = -99.9D+33
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

    !--- Parameters to be controled
    real(8) :: earth_radius      !--- Earth radius
    real(8) :: earth_angvel      !--- Anguler velocity of the earth
    real(8) :: earth_gravity     !--- Gravitational accelaration
    real(8) :: gas_cnst          !--- Gas constant of dry air
    real(8) :: gas_cnst_vap      !--- Gas constant of water vapour
    real(8) :: specific_heat_pre !--- Specific heat of air( const pre )
    real(8) :: specific_heat_pre_vap !--- Specific heat of water vapour ( const pre )
    real(8) :: latent_heat_vap   !--- latent heat of vaporization LH0 ( 0 deg )
    real(8) :: latent_heat_sub   !--- latent heat of sublimation LHS0 ( 0 deg )

    namelist / CNSTPARAM / &
         earth_radius,     & !--- earth radius
         earth_angvel,     & !--- anguler velocity of the earth
         earth_gravity,    & !--- gravitational accelaration
         gas_cnst,         & !--- gas constant of dry air ( Rd )
         gas_cnst_vap,     & !--- gas constant of water vapour
         specific_heat_pre,& !--- specific heat of air ( Cp )
         specific_heat_pre_vap,&  !--- specific heat of water vapour
         latent_heat_vap,  &
         latent_heat_sub

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- initialization of controled parameters
    earth_radius = CNST_ERADIUS
    earth_angvel = CNST_EOHM
    earth_gravity   = CNST_EGRAV
    gas_cnst    = CNST_RAIR
    gas_cnst_vap = CNST_RVAP
    specific_heat_pre      = CNST_CP
    specific_heat_pre_vap  = CNST_CPV
    latent_heat_vap = CNST_LH0
    latent_heat_sub = CNST_LHS0

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
    write(ADM_LOG_FID,CNSTPARAM)

    CNST_ERADIUS = earth_radius
    CNST_EOHM    = earth_angvel
    CNST_EGRAV   = earth_gravity
    CNST_RAIR    = gas_cnst
    CNST_RVAP    = gas_cnst_vap
    CNST_CP      = specific_heat_pre
    CNST_CPV     = specific_heat_pre_vap
    CNST_LH0     = latent_heat_vap
    CNST_LHS0    = latent_heat_sub

    !--- calculate other parameters
    CNST_PI    = 4.D0 * atan( 1.D0 )

    CNST_CV    = CNST_CP - CNST_RAIR
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
