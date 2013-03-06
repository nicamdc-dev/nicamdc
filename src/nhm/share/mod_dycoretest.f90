!-------------------------------------------------------------------------------
!
!+  Module of Dynamical Core Test initial condition
!
!-------------------------------------------------------------------------------
module mod_dycoretest
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
  !      0.00      12-10-19   Imported from mod_restart.f90 of NICAM
  !
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! physical parameters configurations
  real(8), private, parameter :: R = 287.04d0             ! ideal gas constant of dry air
  real(8), private, parameter :: Cp = 1004.D0             ! heat capacity of const. pressure
  real(8), private, parameter :: kai = R / Cp             ! temporal value
  real(8), private, parameter :: g = 9.80616d0            ! gravity accelaration [ms^-2]
  real(8), private, parameter :: pi = 3.14159265358979323846D0
  real(8), private, parameter :: d2r = pi/180.D0          ! Degree to Radian
  real(8), private, parameter :: r2d = 180.D0/pi          ! Radian to Degree
  real(8), private, parameter :: eps = 1.D-14            ! minimum value
  real(8), private, parameter :: eps_r = 1.D-9           ! minimum value (reduce version)
  real(8), private, parameter :: zero = 0.D0              ! zero
  !
  ! test configurations
  integer, private, parameter :: PRCS_D = 8
  ! for Held and Suarez
  real(8), private, save :: deltaT = 60.D0
  real(8), private, save :: deltaTh = 10.D0
  ! for Jablonowski
  real(8), private, save :: clat = 40.D0             ! perturbation center: latitude [deg]
  real(8), private, save :: clon = 20.D0             ! perturbation center: longitude [deg]
  real(8), private, save :: etaS = 1.D0              ! surface eta level
  real(8), private, save :: etaT = 0.2d0              ! threashold of vertical profile
  real(8), private, save :: eta0 = 0.252d0            ! threashold of vertical profile
  real(8), private, save :: t0 = 288.D0              ! [K]
  real(8), private, save :: delT = 4.8d+5             ! [K]
  real(8), private, save :: ganma = 0.005d0           ! [Km^-1]
  real(8), private, save :: u0 = 35.D0               ! [ms^-1]
  real(8), private, save :: uP = 1.D0                ! [ms^-1]
  real(8), private, save :: p0 = 1.D+5               ! [Pa]
  real(8), private, save :: ps = 1.D+5               ! [Pa]
  real(8), private, save :: a = 6.371229d+6           ! [m] mean radis of the earth
  real(8), private, save :: omega = 7.29212d-5        ! [s^-1] rotation of the earth
  logical, private, save :: deep_atm = .false.        ! deep atmosphere setting
  logical, private, parameter :: message = .false.
  integer, private, parameter :: itrmax = 100                  ! # of iteration maximum

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: dycore_input
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
  private :: hs_init
  private :: jbw_init
  private :: tracer_init
  private :: Sp_Unit_East
  private :: Sp_Unit_North
  private :: sphere_xyz_to_lon
  private :: sphere_xyz_to_lat
  private :: eta_vert_coord_NW
  private :: steady_state
  private :: geo2prs
  private :: perturbation
  private :: conv_vxvyvz

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine dycore_input( DIAG_var )
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_gall,      &
       ADM_kall,      &
       ADM_lall
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(out) :: DIAG_var(ADM_gall,ADM_kall,ADM_lall,6+TRC_VMAX)

    character(len=ADM_NSYS) :: init_type

    namelist / DYCORETESTPARAM / &
         init_type

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[dycoretest]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=DYCORETESTPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** DYCORETESTPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist DYCORETESTPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist DYCORETESTPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,DYCORETESTPARAM)

    write(ADM_LOG_FID,*) '*** type: ', trim(init_type)
    if ( trim(init_type) == "Heldsuarez" ) then
       call hs_init ( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )
    elseif( trim(init_type) == "Jablonowski" ) then
       call jbw_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )
    elseif( trim(init_type) == "Traceradvection" ) then ! tentative test
       call tracer_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )
    else
       write(ADM_LOG_FID,*) 'xxx Invalid input_io_mode. STOP.'
       call ADM_proc_stop
    endif

    return
  end subroutine dycore_input

  !-----------------------------------------------------------------------------
  subroutine hs_init( &
     ijdim,   &
     kdim,    &
     lall,    &
     DIAG_var )
    use mod_misc, only: &
       MISC_get_latlon
    use mod_adm, only: &
       ADM_KNONE
    use mod_grd, only: &
       GRD_x,          &
       GRD_XDIR,       &
       GRD_YDIR,       &
       GRD_ZDIR,       &
       GRD_vz
    use mod_cnst, only :  &
       CNST_RAIR, &
         CNST_EGRAV
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: lall
    real(8), intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(8), parameter :: p_surf = 1.D+5         ! [Pa]
    real(8), parameter :: LASPdry = 9.771958146487295D-3
    real(8) :: t_1stlev, t_surf, t_atthelev, t_guess
    real(8) :: p_1stlev(2)
    real(8) :: dZ, dT, f, df
    real(8) :: lat, lon, e

    integer :: n, k, l, itr, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE
    e = exp(1.D0)

    DIAG_var(:,:,:,:) = 0.D0

    do l = 1, lall
    do n = 1, ijdim
       k = 1
       call MISC_get_latlon( lat, lon,               &
                             GRD_x(n,K0,l,GRD_XDIR), &
                             GRD_x(n,K0,l,GRD_YDIR), &
                             GRD_x(n,K0,l,GRD_ZDIR)  )
       t_surf = 315.D0 - deltaT*sin(lat)**2.D0
       p_1stlev(:) = p_surf
       dZ = GRD_vz(n,k,l,1) - 0.D0

       ! Newton-Lapson
       p_1stlev(1) = 0.D0 ! to do the first step
       do itr = 1, itrmax
          if( abs(p_1stlev(1)-p_1stlev(2)) <= eps_r ) exit

          p_1stlev(1) = p_1stlev(2)
          t_1stlev = ( 315.D0 - deltaT*sin(lat)**2 - deltaTh*log(p_1stlev(1)/p0) &
                                            * cos(lat)**2 ) * (p_1stlev(1)/p0)**kai
          f  = ( log( p_surf ) - log( p_1stlev(1) ) ) / ( dZ ) &
             + CNST_EGRAV / CNST_RAIR * 2.D0 / ( t_surf + t_1stlev )
          df = 1.D0 / dZ / p_1stlev(1)

          p_1stlev(2) = p_1stlev(1) + f/df
       enddo
       if ( itr > itrmax ) then
          write(*,*) 'xxx iteration not converged!', &
                     p_1stlev(1)-p_1stlev(2), p_1stlev(2), t_1stlev, t_surf
       endif

       DIAG_var(n,k,l,1) = p_1stlev(2)
       t_atthelev = ( 315.D0 - deltaT*sin(lat)**2 - deltaTh*log(DIAG_var(n,k,l,1)/p0) &
                                         * cos(lat)**2 ) * (DIAG_var(n,k,l,1)/p0)**kai
       DIAG_var(n,k,l,2) = max(200.D0, t_atthelev)

       ! guess pressure field upper k=1
       do k=2, kdim
          dZ = GRD_vz(n,k,l,1) - GRD_vz(n,k-1,l,1)
          t_guess = DIAG_var(n,k-1,l,2) - LASPdry * dZ
          dT = ( t_guess + DIAG_var(n,k-1,l,2) )/2.D0
          DIAG_var(n,k,l,1) = DIAG_var(n,k-1,l,1) * e**( -1.D0 * dZ * (g/(R*dT)) )
          t_atthelev = ( 315.D0 - deltaT*sin(lat)**2 - deltaTh*log(DIAG_var(n,k,l,1)/p0) &
                                                   * cos(lat)**2 ) * (DIAG_var(n,k,l,1)/p0)**kai
          DIAG_var(n,k,l,2) = max(200.D0, t_atthelev)
       enddo

    enddo
    enddo
 
    return
  end subroutine hs_init  

  !-----------------------------------------------------------------------------
  subroutine jbw_init( &
     ijdim,   &
     kdim,    &
     lall,    &
     DIAG_var )
    use mod_misc, only: &
       MISC_get_latlon
    use mod_adm, only: &
       ADM_KNONE
    use mod_grd, only: &
       GRD_vz,         &
       GRD_x,          &
       GRD_x_pl,       &
       GRD_XDIR,       &
       GRD_YDIR,       &
       GRD_ZDIR,       &
       GRD_Z,          &
       GRD_ZH
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: lall
    real(8), intent(out) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    ! work paramters
    real(PRCS_D) :: lat, lon                 ! latitude, longitude on Icosahedral grid
    real(PRCS_D) :: eta(kdim,2), geo(kdim)   ! eta & geopotential in ICO-grid field
    real(PRCS_D) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(PRCS_D) :: wix(kdim),   wiy(kdim)   ! wind components in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)

    logical :: signal ! if ture, continue iteration
    integer :: n, l, k, itr, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    DIAG_var(:,:,:,:) = 0.D0

    do l = 1, lall
    do n = 1, ijdim
       z_local(1) = GRD_vz(n,2,l,GRD_ZH)
       do k = 2, kdim
          z_local(k) = GRD_vz(n,k,l,GRD_Z)
       enddo

       call MISC_get_latlon( lat, lon,               &
                             GRD_x(n,K0,l,GRD_XDIR), &
                             GRD_x(n,K0,l,GRD_YDIR), &
                             GRD_x(n,K0,l,GRD_ZDIR)  )

       signal = .true.
       ! iteration -----
       do itr = 1, itrmax

          if( itr == 1 ) then
             eta(:,:) = 1.D-7 ! This initial value is recommended by Jablonowsky.
          else
             call eta_vert_coord_NW( kdim, itr, z_local, tmp, geo, eta, signal )
          endif

          call steady_state( kdim, lat, eta, wix, wiy, tmp, geo )

          if( .NOT. signal ) exit
       enddo

       if ( itr > itrmax ) then
          write(ADM_LOG_FID,*) 'ETA ITERATION ERROR: NOT CONVERGED', n, l
          stop
       endif

       call geo2prs     ( kdim, tmp, geo, prs )
       call perturbation( kdim, lat, lon, tmp, prs, wix, wiy )
       call conv_vxvyvz ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k=1, kdim
          DIAG_var(n,k,l,1) = prs(k)
          DIAG_var(n,k,l,2) = tmp(k)
          DIAG_var(n,k,l,3) = vx_local(k)
          DIAG_var(n,k,l,4) = vy_local(k)
          DIAG_var(n,k,l,5) = vz_local(k)
       enddo

    enddo
    enddo

    !call surface_height( lall, ijdim, a, omega, g )

    return
  end subroutine jbw_init  

  !-----------------------------------------------------------------------------
  subroutine tracer_init( &
     ijdim,   &
     kdim,    &
     lall,    &
     DIAG_var )
    use mod_misc, only: &
       MISC_get_latlon, &
       MISC_get_distance
    use mod_adm, only: &
       ADM_KNONE
    use mod_cnst, only: &
       CNST_PI, &
       CNST_ERADIUS
    use mod_grd, only: &
       GRD_vz,         &
       GRD_x,          &
       GRD_x_pl,       &
       GRD_XDIR,       &
       GRD_YDIR,       &
       GRD_ZDIR,       &
       GRD_Z,          &
       GRD_ZH
    use mod_gmtr, only: &
       GMTR_lon, &
       GMTR_lat
    use mod_runconf, only: &
       TRC_vmax, &
       NCHEM_STR
    use mod_chemvar, only: &
       chemvar_getid
    implicit none

    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: lall
    real(8), intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(8) :: lon_center   = 120.D0 ! [degree]
    real(8) :: lat_center   =  30.D0 ! [degree]
    real(8) :: alt_center   =  10.D3 ! [m]
    real(8) :: bubble_radius_h = 2000.D3 ! [m]
    real(8) :: bubble_radius_v =    5.D3 ! [m]
    real(8) :: bubble_conc     =  100.D0

    ! work paramters
    real(8) :: dist, dist_h
    real(8) :: rlon_center
    real(8) :: rlat_center

    real(8) :: lat, lon                 ! latitude, longitude on Icosahedral grid
    real(8) :: eta(kdim,2), geo(kdim)   ! eta & geopotential in ICO-grid field
    real(8) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)   ! wind components in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)

    integer :: I_passive1
    logical :: signal ! if ture, continue iteration
    integer :: n, l, k, itr, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    DIAG_var(:,:,:,:) = 0.D0

    I_passive1 = 6 + chemvar_getid( "passive1" ) + NCHEM_STR - 1

    rlon_center = lon_center / 180.D0 * CNST_PI
    rlat_center = lat_center / 180.D0 * CNST_PI

    do l = 1, lall
    do n = 1, ijdim
       z_local(1) = GRD_vz(n,2,l,GRD_ZH)
       do k = 2, kdim
          z_local(k) = GRD_vz(n,k,l,GRD_Z)
       enddo

       call MISC_get_latlon( lat, lon,               &
                             GRD_x(n,K0,l,GRD_XDIR), &
                             GRD_x(n,K0,l,GRD_YDIR), &
                             GRD_x(n,K0,l,GRD_ZDIR)  )

       signal = .true.
       ! iteration -----
       do itr = 1, itrmax

          if( itr == 1 ) then
             eta(:,:) = 1.D-7 ! This initial value is recommended by Jablonowsky.
          else
             call eta_vert_coord_NW( kdim, itr, z_local, tmp, geo, eta, signal )
          endif

          call steady_state( kdim, lat, eta, wix, wiy, tmp, geo )

          if( .NOT. signal ) exit
       enddo

       if ( itr > itrmax ) then
          write(ADM_LOG_FID,*) 'ETA ITERATION ERROR: NOT CONVERGED', n, l
          stop
       endif

       call geo2prs     ( kdim, tmp, geo, prs )
       call conv_vxvyvz ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k=1, kdim
          DIAG_var(n,k,l,1) = prs(k)
          DIAG_var(n,k,l,2) = tmp(k)
          DIAG_var(n,k,l,3) = vx_local(k)
          DIAG_var(n,k,l,4) = vy_local(k)
          DIAG_var(n,k,l,5) = vz_local(k)
       enddo

       call MISC_get_distance( CNST_ERADIUS,  &
                               rlon_center,   &
                               rlat_center,   &
                               GMTR_lon(n,l), &
                               GMTR_lat(n,l), &
                               dist_h         )

       do k = 1, kdim
          dist = ( dist_h                          /bubble_radius_h )**2 &
               + ( (GRD_vz(n,k,l,GRD_Z)-alt_center)/bubble_radius_v )**2

          DIAG_var(n,k,l,I_passive1) = bubble_conc * cos( 0.5D0*PI*sqrt( min(dist,1.D0) ) )**2
       enddo
    enddo
    enddo

    return
  end subroutine tracer_init  

  !-----------------------------------------------------------------------------
  ! eta vertical coordinate by Newton Method
  subroutine eta_vert_coord_NW( &
      kdim,   &  !--- IN : # of z dimension
      itr,    &  !--- IN : iteration number
      z,      &  !--- IN : z-height vertical coordinate
      tmp,    &  !--- IN : guessed temperature
      geo,    &  !--- IN : guessed geopotential
      eta,    &  !--- INOUT : eta level vertical coordinate
      signal  )  !--- INOUT : iteration signal
    !
    implicit none
    integer, intent(in) :: itr
    integer, intent(in) :: kdim
    real(PRCS_D), intent(in) :: z(kdim)
    real(PRCS_D), intent(in) :: geo(kdim), tmp(kdim)
    real(PRCS_D), intent(inout) :: eta(kdim,2)
    logical, intent(inout) :: signal
    integer :: k
    real(PRCS_D) :: diffmax, diff(kdim)
    real(PRCS_D) :: F(kdim), Feta(kdim)
    !
    do k=1, kdim
       F(k) = -g*z(k) + geo(k)
       Feta(k) = -1.D0 * ( R/eta(k,1) ) * tmp(k)
       eta(k,2) = eta(k,1) - ( F(k)/Feta(k) )
       if(eta(k,2) > 1.D0) eta(k,2) = 1.D0     ! not allow over 1.0 for eta
       if(eta(k,2) < 0.D0) eta(k,2) = 1.D-20
       diff(k) = abs( eta(k,2) - eta(k,1) )
    enddo
    !
    eta(:,1) = eta(:,2)
    diffmax = maxval(diff)
    if (message) write (*, '("| Eta  ",I4,": -- MAX: ",F23.20,3X,"MIN: ",F23.20)') itr, maxval(eta(:,1)), minval(eta(:,1))
    if (message) write (*, '("| Diff ",I4,": -- MAX: ",F23.20,3X,"MIN: ",F23.20)') itr, diffmax, minval(diff)
    !
    if(diffmax < eps) then
       signal = .false.
    else
       if (message) write (ADM_LOG_FID,*) "| ----- Iterating ", itr
    endif
    !
    return
  end subroutine eta_vert_coord_NW
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  ! calculation of steady state
  subroutine steady_state( &
      kdim,   &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      eta,  &  !--- IN : eta level vertical coordinate
      wix,  &  !--- INOUT : zonal wind component
      wiy,  &  !--- INOUT : meridional wind component
      tmp,  &  !--- INOUT : mean temperature
      geo   )  !--- INOUT : mean geopotential height
    !
    implicit none
    integer :: k
    integer, intent(in) :: kdim
    real(PRCS_D), intent(in) :: lat
    real(PRCS_D), intent(in) :: eta(kdim,2)
    real(PRCS_D), intent(inout) :: wix(kdim)
    real(PRCS_D), intent(inout) :: wiy(kdim)
    real(PRCS_D), intent(inout) :: tmp(kdim)
    real(PRCS_D), intent(inout) :: geo(kdim)
    real(PRCS_D) :: eta_v
    real(PRCS_D) :: work1, work2
    !
    ! ---------- horizontal mean
    work1 = pi/2.D0
    work2 = R*ganma/g
    do k=1, kdim
       eta_v = (eta(k,1) - eta0)*(work1)
       wix(k) = u0 * (cos(eta_v))**1.5d0 * (sin(2.D0*lat))**2.D0
       !
      !if( etaS >= eta(k,1) .and. eta(k,1) >= etaT ) then  ! not allow over 1.0 for eta
        if( eta(k,1) >= etaT ) then
           tmp(k) = t0 * eta(k,1)**work2
           geo(k) = t0*g/ganma * ( 1.D0 - eta(k,1)**work2 )
        elseif( eta(k,1) < etaT ) then
           tmp(k) = t0 * eta(k,1)**work2 + delT*(etaT - eta(k,1))**5.D0
           !
           geo(k) = t0*g/ganma * ( 1.D0 - eta(k,1)**work2 ) - R * delT *                                     &
                   ( ( log(eta(k,1)/etaT) + 137.D0/60.D0 )*etaT**5.D0 - 5.D0*(etaT**4.D0)*eta(k,1)       &
                     + 5.D0*(etaT**3.D0)*(eta(k,1)**2.D0) - (10.D0/3.D0)*(etaT**2.D0)*(eta(k,1)**3.D0) &
                     + (5.D0/4.D0)*etaT*(eta(k,1)**4.D0) - (1.D0/5.D0)*(eta(k,1)**5.D0)                 &
                   )
        else
           write (ADM_LOG_FID,'(A)') "|-- ETA BOUNDARY ERROR: [steady state calc.]"
           write (ADM_LOG_FID,'("|-- (",I3,")  eta: ",F10.4)') k, eta(k,1)
           stop
        endif
      !else
      !   write (ADM_LOG_FID,'(A)') "|-- OVER 1.0 for eta: [steady state calc.]"
      !   stop
      !endif
       !
    enddo
    !
    ! ---------- meridional distribution for temeperature and geopotential
    work1 = pi/2.D0
    work2 = 3.D0/4.D0 * ( pi*u0 / R )
    do k=1, kdim
       eta_v = (eta(k,1) - eta0)*(work1)
       tmp(k) = tmp(k)                                           &
                    + work2*eta(k,1) * sin(eta_v) * (cos(eta_v))**0.5d0  &
                    * ( ( -2.D0 * (sin(lat))**6.D0 * (cos(lat)**2.D0 + 1.D0/3.D0) + 10.D0/63.D0 )   &
                         * 2.D0*u0*(cos(eta_v))**1.5d0                   &
                        + ( 8.D0/5.D0 * (cos(lat))**3.D0 * ((sin(lat))**2.D0 + 2.D0/3.D0) - pi/4.D0 ) &
                         * a*omega                                       &
                      )
       geo(k) = geo(k)                                           &
                    + u0*(cos(eta_v))**1.5d0  &
                    * ( ( -2.D0 * (sin(lat))**6.D0 * (cos(lat)**2.D0 + 1.D0/3.D0) + 10.D0/63.D0 )   &
                         * u0*(cos(eta_v))**1.5d0                        &
                        + ( 8.D0/5.D0 * (cos(lat))**3.D0 * ((sin(lat))**2.D0 + 2.D0/3.D0) - pi/4.D0 ) &
                         * a*omega                                       &
                      )
    enddo
    !
    wiy(:) = 0.D0
    !
    return
  end subroutine steady_state
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  ! convert geopotential height to pressure
  subroutine geo2prs( &
      kdim,   &  !--- IN : # of z dimension
      tmp,  &  !--- IN : temperature
      geo,  &  !--- IN : geopotential height at full height
      prs   )  !--- INOUT : pressure
    !
    implicit none
    integer :: k
    integer, intent(in) :: kdim
    real(PRCS_D), intent(in) :: tmp(kdim)
    real(PRCS_D), intent(in) :: geo(kdim)
    real(PRCS_D), intent(inout) :: prs(kdim)
    real(PRCS_D) :: dZ, dT, e
    !
    e = exp(1.D0)
    prs(1) = p0
    ! guess pressure field upper k=0
    do k=2, kdim
       dZ = ( geo(k) - geo(k-1) )/g
       dT = ( tmp(k) + tmp(k-1) )/2.D0
       prs(k) = prs(k-1) * e**( -1.D0 * dZ * (g/(R*dT)) )
    enddo
    !
    if (message) then
       write (*,'(A)') "| ----- Pressure (Final Guess) -----"
       do k=1, kdim
          write(ADM_LOG_FID, '("| K(",I3,") -- ",F20.13)') k, prs(k)
       enddo
    endif
    !
    return
  end subroutine geo2prs
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  ! setting perturbation
  subroutine perturbation( &
      kdim,   &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      lon,  &  !--- IN : longitude information
      tmp,  &  !--- INOUT : temperature
      prs,  &  !--- INOUT : pressure
      wix,  &  !--- INOUT : zonal wind component
      wiy   )  !--- INOUT : meridional wind component
    !
    implicit none
    integer :: k
    integer, intent(in) :: kdim
    real(PRCS_D), intent(in) :: lat, lon
    real(PRCS_D), intent(inout) :: tmp(kdim)
    real(PRCS_D), intent(inout) :: prs(kdim)
    real(PRCS_D), intent(inout) :: wix(kdim)
    real(PRCS_D), intent(inout) :: wiy(kdim)
    real(PRCS_D) :: r, rr, rbyrr, cla, clo
    !
    cla = clat * d2r
    clo = clon * d2r
    r = a * acos( sin(cla)*sin(lat) + cos(cla)*cos(lat)*cos(lon-clo) )
    rr = a / 10.D0
    rbyrr = r/rr
    !
    do k=1, kdim
       !tmp(k) = tmp(k)
       !prs(k) = prs(k)
       wix(k) = wix(k) + uP * exp( -1.D0*rbyrr**2.D0 )
       !wiy(k) = wiy(k)
    enddo
    !
    if (message) then
       write (ADM_LOG_FID,'(A)') "|-----------------------------------------------------"
       write (ADM_LOG_FID, '("| T: -- MAX: ",F9.3,3X,"MIN: ",F9.3)') maxval(tmp), minval(tmp)
       write (ADM_LOG_FID, '("| P: -- MAX: ",F9.2,3X,"MIN: ",F9.2)') maxval(prs), minval(prs)
       write (ADM_LOG_FID, '("| U: -- MAX: ",F9.5,3X,"MIN: ",F9.5)') maxval(wix), minval(wix)
       write (ADM_LOG_FID, '("| V: -- MAX: ",F9.5,3X,"MIN: ",F9.5)') maxval(wiy), minval(wiy)
    endif
    !
    return
  end subroutine perturbation
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine conv_vxvyvz( &
      kdim, &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      lon,  &  !--- IN : longitude information
      wix,  &  !--- IN : zonal wind component on latlon
      wiy,  &  !--- IN : meridional wind component on latlon
      vx1d, &  !--- INOUT : horizontal-x component on absolute system for horizontal wind
      vy1d, &  !--- INOUT : horizontal-y component on absolute system for horizontal wind
      vz1d  )  !--- INOUT : vertical component on absolute system for horizontal wind
    !
    implicit none
    integer, intent(in) :: kdim
    real(PRCS_D), intent(in) :: lat, lon
    real(PRCS_D), intent(in) :: wix(kdim)
    real(PRCS_D), intent(in) :: wiy(kdim)
    real(PRCS_D), intent(inout) :: vx1d(kdim)
    real(PRCS_D), intent(inout) :: vy1d(kdim)
    real(PRCS_D), intent(inout) :: vz1d(kdim)
    !
    integer :: k
    real(PRCS_D) :: unit_east(3), unit_north(3)
    !
    ! imported from NICAM/nhm/mkinit/prg_mkinit_ncep.f90 (original written by H.Miura)
    ! *** compute vx, vy, vz as 1-dimensional variables
    do k=1, kdim
       unit_east  = Sp_Unit_East( lon )
       unit_north = Sp_Unit_North( lon, lat )
       !
       vx1d(k) = unit_east(1) * wix(k) + unit_north(1) * wiy(k)
       vy1d(k) = unit_east(2) * wix(k) + unit_north(2) * wiy(k)
       vz1d(k) = unit_east(3) * wix(k) + unit_north(3) * wiy(k)
    enddo
    !
    return
  end subroutine conv_vxvyvz
  !-----------------------------------------------------------------------------
  !
  !  "subroutine surface_height" was moved to share/mod_grd.f90
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
    real(PRCS_D), intent(in) :: lon  ! [ rad ]
    real(PRCS_D) :: unit_east(3)
    !
    unit_east(1) = - sin( lon )    ! --- x-direction
    unit_east(2) =   cos( lon )    ! --- y-direction
    unit_east(3) = 0.0_PRCS_D      ! --- z-direction
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
    real(PRCS_D), intent(in) :: lon, lat  ! [ rad ]
    real(PRCS_D) :: unit_north(3)
    !
    unit_north(1) = - sin( lat ) * cos( lon )    ! --- x-direction
    unit_north(2) = - sin( lat ) * sin( lon )    ! --- y-direction
    unit_north(3) =   cos( lat )                 ! --- z-direction
    return
  end function Sp_Unit_North
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  function sphere_xyz_to_lon (xyz) result (lon)  !(x,y,z) to longitude [-pi,pi].
    real(PRCS_D),intent(in) :: xyz(3) !< (x,y,z)
    real(PRCS_D) :: lon !< longitude [rad]
    real(PRCS_D) :: proj
    !
    proj=sqrt (xyz(1)*xyz(1) + xyz(2)*xyz(2))
    if (proj<eps) then
       lon=0.0_PRCS_D       !# pole points
    else
       lon=atan2 (xyz(2),xyz(1))
    end if
    return
  end function sphere_xyz_to_lon
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  function sphere_xyz_to_lat (xyz) result (lat)  !(x,y,z) to latitude [-pi/2,pi/2].
    real(PRCS_D),intent(in) :: xyz(3) !< (x,y,z)
    real(PRCS_D) :: lat !< latitude [rad]
    real(PRCS_D) :: proj
    !
    proj=sqrt (xyz(1)*xyz(1) + xyz(2)*xyz(2))
    if (proj<eps) then
       lat=sign (0.5_PRCS_D*pi,xyz(3))       !# pole points
    else
       lat=atan (xyz(3)/proj)
    end if
    return
  end function sphere_xyz_to_lat
  !-----------------------------------------------------------------------------
  !
end module mod_dycoretest
!-------------------------------------------------------------------------------
