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
!  use mod_adm, only: &
!     ADM_LOG_FID,  &
!     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! physical parameters configurations
    real(8), parameter :: &
         a	= 6371220.0d0,	& ! Earth's Radius (m)
         Rd 	= 287.0d0,	& ! Ideal gas const dry air (J kg^-1 K^1)
         g	= 9.80616d0,	& ! Gravity (m s^2)
         cp	= 1004.5d0,	& ! Specific heat capacity (J kg^-1 K^1)
         pi	= 4.d0*atan(1.d0) ! pi

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: test11_velocity
  public :: test12_velocity
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: Sp_Unit_East
  private :: Sp_Unit_North
  !-----------------------------------------------------------------------------
contains
  !
  !-----------------------------------------------------------------------------
  subroutine test11_velocity (time,lon,lat,zf,zh,vx,vy,vz,w)
    !----- [add; original by H.Miura] 20130612 R.Yoshida
    implicit none
    !
    real(8),intent(in) :: time
    real(8),intent(in) :: zf,zh    ! height (real altitude)
    real(8),intent(in) :: lon
    real(8),intent(in) :: lat
    real(8),intent(out) :: vx
    real(8),intent(out) :: vy
    real(8),intent(out) :: vz
    real(8),intent(out) :: w
    !
    ! <constant>
    real(8), parameter :: 	&
         tau     = 12.d0 * 86400.d0,	& ! period of motion 12 days
         u0      = (2.d0*pi*a)/tau,	& ! 2 pi a / 12 days
         k0	   = (10.d0*a)/tau,	& ! Velocity Magnitude
         omega0  = (23000.d0*pi)/tau,	& ! Velocity Magnitude
         T0      = 300.d0,		& ! temperature
         H       = Rd * T0 / g,	& ! scale height
         RR      = 1.d0/2.d0,		& ! horizontal half width divided by 'a'
         ZZ      = 1000.d0,		& ! vertical half width
         z0      = 5000.d0,		& ! center point in z
         lambda0 = 5.d0*pi/6.d0,	& ! center point in longitudes
         lambda1 = 7.d0*pi/6.d0,	& ! center point in longitudes
         phi0    = 0.d0,		& ! center point in latitudes
         phi1    = 0.d0        
    real(8), parameter ::   p0 = 100000.d0 ! reference pressure (Pa)
    ! <work>
    real(8) :: &
         u, &	! Zonal wind (m s^-1)
         v    ! Meridional wind (m s^-1)
    real(8) :: lonp,height,p,ptop,bs,s,ud,east(3),nrth(3)
    !-----

      ! full level
      height = zf
      p = p0 * exp(-zf/H)

      ptop    = p0*exp(-12000.d0/H)

      lonp = lon - 2.d0*pi*time/tau

      bs = 0.2
      s = 1.0 + exp( (ptop-p0)/(bs*ptop) ) &
              - exp( (p-p0)/(bs*ptop)) &
              - exp( (ptop-p)/(bs*ptop))

      ud = (omega0*a)/(bs*ptop) * cos(lonp) * (cos(lat)**2.0) * cos(2.0*pi*time/tau) * &
           ( - exp( (p-p0)/(bs*ptop)) + exp( (ptop-p)/(bs*ptop))  )

      u = k0*sin(lonp)*sin(lonp)*sin(2.d0*lat)*cos(pi*time/tau) + u0*cos(lat) + ud
      v = k0*sin(2.d0*lonp)*cos(lat)*cos(pi*time/tau)

      east=Sp_Unit_East (lon)
      nrth=Sp_Unit_North (lon,lat) 
      
      vx=east(1)*u+nrth(1)*v
      vy=east(2)*u+nrth(2)*v
      vz=east(3)*u+nrth(3)*v

      ! half level
      height = zh
      p = p0 * exp(-zh/H)

      ptop    = p0*exp(-12000.d0/H)

      lonp = lon - 2.d0*pi*time/tau

      bs = 0.2
      s = 1.0 + exp( (ptop-p0)/(bs*ptop) ) &
              - exp( (p-p0)/(bs*ptop)) &
              - exp( (ptop-p)/(bs*ptop))

      w = -((Rd*T0)/(g*p))*omega0*sin(lonp)*cos(lat)*cos(2.0*pi*time/tau)*s

      return
  end subroutine test11_velocity
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine test12_velocity (time,lon,lat,zf,zh,vx,vy,vz,w)
    !----- [add; original by H.Miura] 20130612 R.Yoshida
    implicit none
    !
    real(8),intent(in) :: time
    real(8),intent(in) :: zf,zh    ! height (real altitude)
    real(8),intent(in) :: lon
    real(8),intent(in) :: lat
    real(8),intent(out) :: vx
    real(8),intent(out) :: vy
    real(8),intent(out) :: vz
    real(8),intent(out) :: w
    !
    ! <constant>
    real(8), parameter :: &
         tau     = 1.d0 * 86400.d0, &	! period of motion 1 day (in s)
         u0      = 40.d0,           &	! Zonal velocity magnitude (m/s)
         w0	   = 0.15d0,          &	! Vertical velocity magnitude (m/s), changed in v5
         T0      = 300.d0,          &	! temperature (K)
         H       = Rd * T0 / g,     &	! scale height
         K       = 5.d0,            &	! number of Hadley-like cells
         z1      = 2000.d0,    & ! position of lower tracer bound (m), changed in v5
         z2      = 5000.d0,    & ! position of upper tracer bound (m), changed in v5
         z0      = 0.5d0*(z1+z2),   &	! midpoint (m)       
         ztop    = 12000.d0		! model top (m)
    real(8), parameter ::   p0 = 100000.d0 ! reference pressure (Pa)
    ! <work>
    real(8) :: &
         u, &	! Zonal wind (m s^-1)
         v    ! Meridional wind (m s^-1)
    real(8) :: lonp,height,p,ptop,bs,s,ud,east(3),nrth(3),rho,rho0,t,phis
    !-----

      ! full level
      height = zf
      p = p0 * exp(-zf/H)

      t = T0
      phis = 0.d0

      rho = p/(Rd*t)
      rho0 = p0/(Rd*t)

      u = u0*cos(lat)
      v = -(rho0/rho) * (a*w0*pi)/(K*ztop) *cos(lat)*sin(K*lat)&
                         *cos(pi*height/ztop)*cos(pi*time/tau)

      east=Sp_Unit_East (lon)
      nrth=Sp_Unit_North (lon,lat) 
      
      vx=east(1)*u+nrth(1)*v
      vy=east(2)*u+nrth(2)*v
      vz=east(3)*u+nrth(3)*v

      ! half level
      height = zh
      p = p0 * exp(-zf/H)

      t = T0
      phis = 0.d0

      rho = p/(Rd*t)
      rho0 = p0/(Rd*t)

      w = (rho0/rho) *(w0/K)*(-2.d0*sin(K*lat)*sin(lat) + K*cos(lat)*cos(K*lat)) &
		*sin(pi*height/ztop)*cos(pi*time/tau)

      return
  end subroutine test12_velocity
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  function Sp_Unit_East( lon ) result( unit_east )
    ! imported from prg_mkinit_ncep.f90 (original written by H.Miura)
    ! ------
    ! Compute local eastward unit vector (unit_east)
    ! in the Cartesian coordinate system
    ! given longitude (lon) of a position.
    !
    implicit none
    real(8), intent(in) :: lon  ! [ rad ]
    real(8) :: unit_east(3)
    !
    unit_east(1) = - sin( lon )    ! --- x-direction
    unit_east(2) =   cos( lon )    ! --- y-direction
    unit_east(3) = 0.0_8           ! --- z-direction
    return
  end function Sp_Unit_East
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  function Sp_Unit_North( lon, lat ) result( unit_north )
    ! imported from prg_mkinit_ncep.f90 (original written by H.Miura)
    ! ------
    ! Compute local northward unit vector (unit_north) 
    ! in the Cartesian coordinate system
    ! given longitude (lon) and latitude (lat) of a position.
    !
    implicit none
    real(8), intent(in) :: lon, lat  ! [ rad ]
    real(8) :: unit_north(3)
    !
    unit_north(1) = - sin( lat ) * cos( lon )    ! --- x-direction
    unit_north(2) = - sin( lat ) * sin( lon )    ! --- y-direction
    unit_north(3) =   cos( lat )                 ! --- z-direction
    return
  end function Sp_Unit_North
  !-----------------------------------------------------------------------------
  !
end module mod_af_trcadv
!-------------------------------------------------------------------------------
