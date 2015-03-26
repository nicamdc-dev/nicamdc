!-------------------------------------------------------------------------------
!
!+  Module of Updating UVW for Tracer Advection Test
!
!-------------------------------------------------------------------------------
module mod_af_trcadv
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module is for the Dyn Core Test Initialization.
  !
  !
  !++ Current Corresponding Author : R.Yoshida
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      12-10-19   artificial updating  R.Yoshida
  !
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: test11_velocity
  public :: test12_velocity

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: Sp_Unit_East
  private :: Sp_Unit_North

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  REAL(RP), private, parameter :: a  = 6371220.D0        ! Earth's Radius [m]
  REAL(RP), private, parameter :: Rd = 287.D0            ! Ideal gas const dry air [J/kg*K]
  REAL(RP), private, parameter :: g  = 9.80616D0         ! Gravity [m/s2]
  REAL(RP), private, parameter :: cp = 1004.5D0          ! Specific heat capacity [J/kg*K]
  REAL(RP), private, parameter :: pi = 4.D0 * atan(1.D0) ! pi

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> original by H.Miura, 20130612 R.Yoshida
  subroutine test11_velocity( &
       time, &
       lon,  &
       lat,  &
       zf,   &
       zh,   &
       vx,   &
       vy,   &
       vz,   &
       w     )
    implicit none

    REAL(RP),intent(in)  :: time
    REAL(RP),intent(in)  :: lon
    REAL(RP),intent(in)  :: lat
    REAL(RP),intent(in)  :: zf
    REAL(RP),intent(in)  :: zh
    REAL(RP),intent(out) :: vx
    REAL(RP),intent(out) :: vy
    REAL(RP),intent(out) :: vz
    REAL(RP),intent(out) :: w

    REAL(RP), parameter :: tau     = 12.D0 * 86400.D0 ! period of motion 12 days
    REAL(RP), parameter :: u0      = 2.D0*pi*a/tau    ! 2 pi a / 12 days
    REAL(RP), parameter :: k0      = 10.D0*a/tau      ! Velocity Magnitude
    REAL(RP), parameter :: omega0  = 23000.D0*pi/tau  ! Velocity Magnitude
    REAL(RP), parameter :: T0      = 300.D0           ! temperature
    REAL(RP), parameter :: H       = Rd * T0 / g      ! scale height
    REAL(RP), parameter :: p0      = 100000.D0        ! reference pressure (Pa)

    REAL(RP) :: u ! Zonal wind      [m/s]
    REAL(RP) :: v ! Meridional wind [m/s]
    REAL(RP) :: dlon, lonp, bs, height, p, ptop, s, ud
    REAL(RP) :: east(3), nrth(3)
    !---------------------------------------------------------------------------

    dlon = 2.D0 * pi * time / tau
    lonp = lon - dlon
    bs   = 0.2D0

    ! full level
    height = zf
    p      = p0 * exp(-zf/H)
    ptop   = p0 * exp(-12000.D0/H)

    s = 1.D0 + exp( (ptop-p0) / (bs*ptop) ) &
             - exp( (p-p0)    / (bs*ptop) ) &
             - exp( (ptop-p)  / (bs*ptop) )

    ud = (omega0*a) / (bs*ptop) * cos(lonp) * cos(lat)**2 * cos(dlon) &
       * ( -exp( (p-p0)/(bs*ptop) ) + exp( (ptop-p)/(bs*ptop) ) )

    u = k0 * sin(2.D0*lat)  * sin(lonp)**2 * cos(0.5D0*dlon) + u0 * cos(lat) + ud
    v = k0 * sin(2.D0*lonp) * cos(lat)     * cos(0.5D0*dlon)

    east = Sp_Unit_East (lon)
    nrth = Sp_Unit_North(lon,lat)

    vx = east(1) * u + nrth(1) * v
    vy = east(2) * u + nrth(2) * v
    vz = east(3) * u + nrth(3) * v

    height = zh
    p      = p0 * exp(-zh/H)
    ptop   = p0 * exp(-12000.D0/H)

    s = 1.0 + exp( (ptop-p0) / (bs*ptop) ) &
            - exp( (p-p0)    / (bs*ptop) ) &
            - exp( (ptop-p)  / (bs*ptop) )

    w = -(Rd*T0) / (g*p) * omega0 * sin(lonp) * cos(lat) * cos(dlon) * s

    return
  end subroutine test11_velocity

  !-----------------------------------------------------------------------------
  !> original by H.Miura, 20130612 R.Yoshida
  subroutine test12_velocity( &
       time, &
       lon,  &
       lat,  &
       zf,   &
       zh,   &
       vx,   &
       vy,   &
       vz,   &
       w     )
    implicit none

    REAL(RP),intent(in)  :: time
    REAL(RP),intent(in)  :: lon
    REAL(RP),intent(in)  :: lat
    REAL(RP),intent(in)  :: zf
    REAL(RP),intent(in)  :: zh
    REAL(RP),intent(out) :: vx
    REAL(RP),intent(out) :: vy
    REAL(RP),intent(out) :: vz
    REAL(RP),intent(out) :: w

    REAL(RP), parameter :: tau  = 1.D0 * 86400.D0 ! period of motion 1 day (in s)
    REAL(RP), parameter :: u0   = 40.D0           ! Zonal velocity magnitude (m/s)
    REAL(RP), parameter :: w0   = 0.15D0          ! Vertical velocity magnitude (m/s), changed in v5
    REAL(RP), parameter :: T0   = 300.D0          ! temperature
    REAL(RP), parameter :: H    = Rd * T0 / g     ! scale height
    REAL(RP), parameter :: K    = 5.D0            ! number of Hadley-like cells
    REAL(RP), parameter :: z1   = 2000.D0         ! position of lower tracer bound (m), changed in v5
    REAL(RP), parameter :: z2   = 5000.D0         ! position of upper tracer bound (m), changed in v5
    REAL(RP), parameter :: ztop = 12000.D0        ! model top (m)
    REAL(RP), parameter :: p0   = 100000.D0       ! reference pressure (Pa)

    REAL(RP) :: u ! Zonal wind      [m/s]
    REAL(RP) :: v ! Meridional wind [m/s]
    REAL(RP) :: height, p, rho, t, rho0
    REAL(RP) :: east(3), nrth(3)
    !---------------------------------------------------------------------------

    t    = T0
    rho0 = p0 / (Rd*t)

    height = zf
    p      = p0 * exp(-zf/H)
    rho    = p / (Rd*t)

    u = u0 * cos(lat)
    v = -(rho0/rho) * a*w0*pi / (K*ztop) * cos(lat) * sin(K*lat) &
      * cos(pi*height/ztop) * cos(pi*time/tau)

    east = Sp_Unit_East (lon)
    nrth = Sp_Unit_North(lon,lat)

    vx = east(1) * u + nrth(1) * v
    vy = east(2) * u + nrth(2) * v
    vz = east(3) * u + nrth(3) * v

    height = zh
    p      = p0 * exp(-zh/H)
    rho    = p / (Rd*t)

    w = (rho0/rho) * (w0/K) * (-2.D0*sin(K*lat)*sin(lat)+K*cos(lat)*cos(K*lat)) &
      * sin(pi*height/ztop) * cos(pi*time/tau)

    return
  end subroutine test12_velocity

  !-----------------------------------------------------------------------------
  function Sp_Unit_East( lon ) result( unit_east )
    implicit none

    REAL(RP), intent(in) :: lon ! [rad]
    REAL(RP)             :: unit_east(3)
    !---------------------------------------------------------------------------

    unit_east(1) = -sin(lon) ! x-direction
    unit_east(2) =  cos(lon) ! y-direction
    unit_east(3) = 0.D0      ! z-direction

    return
  end function Sp_Unit_East

  !-----------------------------------------------------------------------------
  function Sp_Unit_North( lon, lat ) result( unit_north )
    implicit none

    REAL(RP), intent(in) :: lon, lat ! [rad]
    REAL(RP)             :: unit_north(3)
    !---------------------------------------------------------------------------

    unit_north(1) = -sin(lat) * cos(lon) ! x-direction
    unit_north(2) = -sin(lat) * sin(lon) ! y-direction
    unit_north(3) =  cos(lat)            ! z-direction

    return
  end function Sp_Unit_North

end module mod_af_trcadv
