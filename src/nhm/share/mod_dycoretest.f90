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
  real(8), private, parameter :: Rd = 287.04d0            ! ideal gas constant of dry air
  real(8), private, parameter :: Cp = 1004.D0             ! heat capacity of const. pressure
  real(8), private, parameter :: kai = R / Cp             ! temporal value
  real(8), private, parameter :: g = 9.80616d0            ! gravity accelaration [ms^-2]
  real(8), private, parameter :: pi = 3.14159265358979323846D0
  real(8), private, parameter :: d2r = pi/180.D0          ! Degree to Radian
  real(8), private, parameter :: r2d = 180.D0/pi          ! Radian to Degree
  real(8), private, parameter :: eps = 1.D-8              ! minimum value
  real(8), private, parameter :: zero = 0.D0              ! zero

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
  private :: mountwave_init
  private :: gravwave_init
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

    character(len=ADM_NSYS) :: init_type, test_case

    namelist / DYCORETESTPARAM / &
         init_type, test_case

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
       call jbw_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )
    elseif( trim(init_type) == "Traceradvection" ) then ! tentative test
       call tracer_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )
    elseif( trim(init_type) == "Mountainwave" ) then
       call mountwave_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )
    elseif( trim(init_type) == "Gravitywave" ) then
       call gravwave_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )
    else
       write(ADM_LOG_FID,*) 'xxx Invalid input_io_mode. STOP.'
       call ADM_proc_stop
    endif

    if ( trim(init_type) == "Jablonowski"     .or. &
         trim(init_type) == "Traceradvection" .or. &
         trim(init_type) == "Mountainwave"         ) then
       write(ADM_LOG_FID,*) '*** test case: ', trim(test_case)
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
       ADM_KNONE, &
       ADM_kmin,  &
       ADM_kmax
    use mod_cnst, only :  &
       GRAV  => CNST_EGRAV, &
       RAIR  => CNST_RAIR,  &
       KAPPA => CNST_KAPPA,  &
       PRE00 => CNST_PRE00
    use mod_grd, only: &
       GRD_x,    &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_vz,   &
       GRD_Z
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: lall
    real(8), intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(8) :: pre(kdim), tem(kdim), dz(kdim)
    real(8) :: pre_sfc, tem_sfc
    real(8) :: pre_save

    real(8) :: dT, f, df
    real(8) :: lat, lon

    integer :: n, k, l, itr
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.D0

    do l = 1, lall
    do n = 1, ijdim

       dz(ADM_kmin) = GRD_vz(n,ADM_kmin,l,GRD_Z)
       do k = ADM_kmin+1, ADM_kmax+1
          dz(k) = GRD_vz(n,k,l,GRD_Z) - GRD_vz(n,k-1,l,GRD_Z)
       enddo

       call MISC_get_latlon( lat, lon,               &
                             GRD_x(n,ADM_KNONE,l,GRD_XDIR), &
                             GRD_x(n,ADM_KNONE,l,GRD_YDIR), &
                             GRD_x(n,ADM_KNONE,l,GRD_ZDIR)  )

       pre_sfc = PRE00
!       tem_sfc = 300.D0
       tem_sfc = 315.D0 - deltaT*sin(lat)**2

       !---< from ground surface to lowermost atmosphere >---
       k = ADM_kmin
       ! first guess
       pre(k) = pre_sfc
       tem(k) = tem_sfc

       ! Newton-Lapson
       do itr = 1, itrmax
          pre_save = pre(k) ! save

          f  = log(pre(k)/pre_sfc) / dz(k) + GRAV / ( RAIR * 0.5D0 * (tem(k)+tem_sfc) )
          df = 1.D0 / (pre(k)*dz(k))

          pre(k) = pre(k) - f / df
!          tem(k) = 300.D0 * ( pre(k)/PRE00 )**KAPPA
          tem(k) = ( 315.D0 - deltaT*sin(lat)**2 - deltaTh*log(pre(k)/PRE00)*cos(lat)**2 ) &
                 * ( pre(k)/PRE00 )**KAPPA
          tem(k) = max( 200.D0, tem(k) )

          if( abs(pre_save-pre(k)) <= eps ) exit
       enddo

       if ( itr > itrmax ) then
          write(ADM_LOG_FID,*) 'xxx iteration not converged!', k, pre_save-pre(k), pre(k), pre_sfc, tem(k), tem_sfc
          write(*,          *) 'xxx iteration not converged!', k, pre_save-pre(k), pre(k), pre_sfc, tem(k), tem_sfc
          stop
       endif

       !---< from lowermost to uppermost atmosphere >---
       do k = ADM_kmin+1, ADM_kmax+1

          ! first guess
          pre(k) = pre(k-1)
          tem(k) = 300.D0 * ( pre(k)/PRE00 )**KAPPA
          tem(k) = max( 200.D0, tem(k) )

          ! Newton-Lapson
          do itr = 1, itrmax
             pre_save = pre(k) ! save

             f  = log(pre(k)/pre(k-1)) / dz(k) + GRAV / ( RAIR * 0.5D0 * (tem(k)+tem(k-1)) )
             df = 1.D0 / (pre(k)*dz(k))

             pre(k) = pre(k) - f / df
!             tem(k) = 300.D0 * ( pre(k)/PRE00 )**KAPPA
             tem(k) = ( 315.D0 - deltaT*sin(lat)**2 - deltaTh*log(pre(k)/PRE00)*cos(lat)**2 ) &
                    * ( pre(k)/PRE00 )**KAPPA
             tem(k) = max( 200.D0, tem(k) )

             if( abs(pre_save-pre(k)) <= eps ) exit

          enddo

          if ( itr > itrmax ) then
             write(ADM_LOG_FID,*) 'xxx iteration not converged!', k, pre_save-pre(k), pre(k), pre(k-1), tem(k), tem(k-1)
             write(*,          *) 'xxx iteration not converged!', k, pre_save-pre(k), pre(k), pre(k-1), tem(k), tem(k-1)
             stop
          endif
       enddo

       DIAG_var(n,ADM_kmin-1,l,1) = pre_sfc ! tentative
       DIAG_var(n,ADM_kmin-1,l,2) = tem_sfc ! tentative
       do k = ADM_kmin, ADM_kmax+1
          DIAG_var(n,k,l,1) = pre(k)
          DIAG_var(n,k,l,2) = tem(k)
       enddo

    enddo
    enddo

    return
  end subroutine hs_init  

  !-----------------------------------------------------------------------------
  subroutine jbw_init( &
     ijdim,      &
     kdim,       &
     lall,       &
     test_case,  &
     DIAG_var    )
    use mod_misc, only: &
       MISC_get_latlon
    use mod_adm, only: &
       ADM_KNONE,      &
       ADM_NSYS
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
    character(len=ADM_NSYS), intent(in) :: test_case
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
    logical :: pertb  ! if ture, with perturbation
    integer :: n, l, k, itr, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    DIAG_var(:,:,:,:) = 0.D0

    select case( trim(test_case) )
    case ('1', '4-1')  ! with perturbation
       write(ADM_LOG_FID,*) "Jablonowski Initialize - case 1: with perturbation"
       pertb = .true.
    case ('2', '4-2')  ! without perturbation
       write(ADM_LOG_FID,*) "Jablonowski Initialize - case 2: without perturbation"
       pertb = .false.
    case default
       write(ADM_LOG_FID,*) "Unknown test_case: '"//trim(test_case)//"' specified."
       write(ADM_LOG_FID,*) "Force changed to case 1 (with perturbation)"
       pertb = .true.
    end select

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
       call perturbation( kdim, lat, lon, tmp, prs, wix, wiy, pertb )
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
  subroutine mountwave_init( &
     ijdim,      &
     kdim,       &
     lall,       &
     test_case,  &
     DIAG_var    )
    use mod_misc, only: &
       MISC_get_latlon, &
       MISC_get_distance
    use mod_adm, only: &
       ADM_KNONE,      &
       ADM_NSYS,       &
       ADM_proc_stop
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
    character(len=ADM_NSYS), intent(in) :: test_case
    real(8), intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    ! work paramters
    real(8) :: lat, lon                 ! latitude, longitude on Icosahedral grid
    real(8) :: eta(kdim,2), geo(kdim)   ! eta & geopotential in ICO-grid field
    real(8) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)   ! wind components in ICO-grid field
    real(8) :: wiz(kdim)                ! vertical wind components in ICO-grid field
    real(8) :: q(kdim),     rho(kdim)   ! tracer and rho in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)

    integer :: I_passive1
    integer :: n, l, k, K0

    integer :: shear
    logical :: fault = .true.
    logical :: hybrid_eta = .false.
    real(8) :: hyam, hybm, phis, ps
    integer, parameter :: zcoords = 1
    !---------------------------------------------------------------------------

    hyam = 0.0d0
    hybm = 0.0d0
    fault = .false.
    hybrid_eta = .false.

    K0 = ADM_KNONE

    DIAG_var(:,:,:,:) = 0.D0

    I_passive1 = 6 + chemvar_getid( "passive1" ) + NCHEM_STR - 1

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

       select case( trim(test_case) )
       ! DCMIP: TEST CASE 2-0 - Steady-State Atmosphere at Rest in the Presence of Orography
       case ('0', '2-0')
          do k=1, kdim
             call test2_steady_state_mountain (lon, lat, prs(k), z_local(k), zcoords, &
                      hybrid_eta, hyam, hybm, wix(k), wiy(k), wiz(k), tmp(k), phis, &
                      ps, rho(k), q(k) )
          enddo
       ! DCMIP: TEST CASE 2-1 - Non-hydrostatic Mountain Waves over a Schaer-type Mountain
       case ('1', '2-1')
          shear = 0   ! test case: 2-1 (constant u)
          do k=1, kdim
             call test2_schaer_mountain (lon, lat, prs(k), z_local(k), zcoords, &
                      hybrid_eta, hyam, hybm, shear, wix(k), wiy(k), wiz(k), tmp(k), phis, &
                      ps, rho(k), q(k) )
          enddo
       ! DCMIP: TEST CASE 2-2 - Non-hydrostatic Mountain Waves over a Schaer-type Mountain
       case ('2', '2-2')
          shear = 1   ! test case: 2-2 (sheared u)
          do k=1, kdim
             call test2_schaer_mountain (lon, lat, prs(k), z_local(k), zcoords, &
                      hybrid_eta, hyam, hybm, shear, wix(k), wiy(k), wiz(k), tmp(k), phis, &
                      ps, rho(k), q(k) )
          enddo
       case default
          write(ADM_LOG_FID,*) "Unknown test_case: '"//trim(test_case)//"' specified."
          call ADM_proc_stop
       end select

       call conv_vxvyvz ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k=1, kdim
          DIAG_var(n,k,l,1) = prs(k)
          DIAG_var(n,k,l,2) = tmp(k)
          DIAG_var(n,k,l,3) = vx_local(k)
          DIAG_var(n,k,l,4) = vy_local(k)
          DIAG_var(n,k,l,5) = vz_local(k)
          DIAG_var(n,k,l,6) = wiz(k)
          DIAG_var(n,k,l,I_passive1) = q(k)
       enddo
    enddo
    enddo

    return
  end subroutine mountwave_init

  !-----------------------------------------------------------------------------
  subroutine gravwave_init( &
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

    ! work paramters
    real(8) :: lat, lon                 ! latitude, longitude on Icosahedral grid
    real(8) :: eta(kdim,2), geo(kdim)   ! eta & geopotential in ICO-grid field
    real(8) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)   ! wind components in ICO-grid field
    real(8) :: wiz(kdim)                ! vertical wind components in ICO-grid field
    real(8) :: q(kdim),     rho(kdim)   ! tracer and rho in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)

    integer :: I_passive1
    integer :: n, l, k, K0

    integer :: shear
    logical :: fault = .true.
    logical :: hybrid_eta = .false.
    real(8) :: hyam, hybm, phis, ps
    integer, parameter :: zcoords = 1
    !---------------------------------------------------------------------------

    hyam = 0.0d0
    hybm = 0.0d0
    fault = .false.
    hybrid_eta = .false.

    K0 = ADM_KNONE

    DIAG_var(:,:,:,:) = 0.D0

    I_passive1 = 6 + chemvar_getid( "passive1" ) + NCHEM_STR - 1

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

       do k=1, kdim
          call test3_gravity_wave (lon, lat, prs(k), z_local(k), zcoords, &
                   wix(k), wiy(k), wiz(k), tmp(k), phis, &
                   ps, rho(k), q(k) )
       enddo

       call conv_vxvyvz ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k=1, kdim
          DIAG_var(n,k,l,1) = prs(k)
          DIAG_var(n,k,l,2) = tmp(k)
          DIAG_var(n,k,l,3) = vx_local(k)
          DIAG_var(n,k,l,4) = vy_local(k)
          DIAG_var(n,k,l,5) = vz_local(k)
          DIAG_var(n,k,l,6) = wiz(k)
          DIAG_var(n,k,l,I_passive1) = q(k)
       enddo
    enddo
    enddo

    return
  end subroutine gravwave_init

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
  ! setting perturbation
  subroutine perturbation( &
      kdim,   &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      lon,  &  !--- IN : longitude information
      tmp,  &  !--- INOUT : temperature
      prs,  &  !--- INOUT : pressure
      wix,  &  !--- INOUT : zonal wind component
      wiy,  &  !--- INOUT : meridional wind component
      with  )  !--- IN : perturbation switch
    !
    implicit none
    integer :: k
    integer, intent(in) :: kdim
    logical, intent(in) :: with
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
    if ( with ) then
       do k=1, kdim
          !tmp(k) = tmp(k)
          !prs(k) = prs(k)
          wix(k) = wix(k) + uP * exp( -1.D0*rbyrr**2.D0 )
          !wiy(k) = wiy(k)
       enddo
    endif
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


!=========================================================================
! Test 2-0:  Steady-State Atmosphere at Rest in the Presence of Orography
!=========================================================================
SUBROUTINE test2_steady_state_mountain (lon,lat,p,z,zcoords,hybrid_eta,hyam,hybm,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z, &		! Height (m)
				hyam, &		! A coefficient for hybrid-eta coordinate, at model level midpoint
				hybm		! B coefficient for hybrid-eta coordinate, at model level midpoint

	logical, intent(in)  :: hybrid_eta      ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used
                                                ! if set to .true., then the pressure will be computed via the 
                                                !    hybrid coefficients hyam and hybm, they need to be initialized
                                                ! if set to .false. (for pressure-based models): the pressure is already pre-computed
                                                !    and is an input value for this routine 
                                                ! for height-based models: pressure will always be computed based on the height and
                                                !    hybrid_eta is not used

	real(8), intent(inout) :: p		! Pressure  (Pa)
				
	integer,  intent(in) :: zcoords 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q 		! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we compute or use p
        !
	! In hybrid-eta coords: p = hyam p0 + hybm ps
        !
        ! The grid-point based initial data are computed in this routine. 

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: 	T0      = 300.d0,		&	! temperature (K)
                                gamma   = 0.0065d0,             &       ! temperature lapse rate (K/m)      
                            	lambdam = 3.d0*pi/2.d0,		&	! mountain longitude center point (radians)   
                            	phim    = 0.d0,			&	! mountain latitude center point (radians)    
                            	h0      = 2000.d0,		&	! peak height of the mountain range (m)
                            	Rm      = 3.d0*pi/4.d0,		&	! mountain radius (radians)
                            	zetam   = pi/16.d0,		&	! mountain oscillation half-width (radians) 
                            	ztop    = 12000.d0			! model top (m)         
                            
      real(8) :: height							! Model level heights (m)
      real(8) :: r							! Great circle distance (radians)
      real(8) :: zs							! Surface elevation (m)
      real(8) :: exponent                                               ! exponent: g/(Rd * gamma)
      real(8) :: exponent_rev                                           ! reversed exponent


!-----------------------------------------------------------------------
!    compute exponents 
!-----------------------------------------------------------------------
     exponent     = g/(Rd*gamma)
     exponent_rev = 1.d0/exponent

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------
      
	r = acos( sin(phim)*sin(lat) + cos(phim)*cos(lat)*cos(lon - lambdam) )

	if (r .lt. Rm) then

		zs = (h0/2.d0)*(1.d0+cos(pi*r/Rm))*cos(pi*r/zetam)**2.d0   ! mountain height

	else

		zs = 0.d0

	endif

	phis = g*zs

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = p0 * (1.d0 - gamma/T0*zs)**exponent


!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------

	! Height and pressure are aligned (p = p0 * (1.d0 - gamma/T0*z)**exponent)

	if (zcoords .eq. 1) then

		height = z
		p = p0 * (1.d0 - gamma/T0*z)**exponent

	else

                if (hybrid_eta) p = hyam*p0 + hybm*ps              ! compute the pressure based on the surface pressure and hybrid coefficients
		height = T0/gamma * (1.d0 - (p/p0)**exponent_rev)  ! compute the height at this pressure

	endif

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE ZERO (STATE AT REST)
!-----------------------------------------------------------------------

	! Zonal Velocity

	u = 0.d0

	! Meridional Velocity

	v = 0.d0

        ! Vertical Velocity

        w = 0.d0

!-----------------------------------------------------------------------
!    TEMPERATURE WITH CONSTANT LAPSE RATE
!-----------------------------------------------------------------------

	t = T0 - gamma*height

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE test2_steady_state_mountain



!=====================================================================================
! Tests 2-1 and 2-2:  Non-hydrostatic Mountain Waves over a Schaer-type Mountain
!=====================================================================================

SUBROUTINE test2_schaer_mountain (lon,lat,p,z,zcoords,hybrid_eta,hyam,hybm,shear,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z,   &		! Height (m)
				hyam, &		! A coefficient for hybrid-eta coordinate, at model level midpoint
				hybm		! B coefficient for hybrid-eta coordinate, at model level midpoint

	logical, intent(in)  :: hybrid_eta      ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used
                                                ! if set to .true., then the pressure will be computed via the 
                                                !    hybrid coefficients hyam and hybm, they need to be initialized
                                                ! if set to .false. (for pressure-based models): the pressure is already pre-computed
                                                !    and is an input value for this routine 
                                                ! for height-based models: pressure will always be computed based on the height and
                                                !    hybrid_eta is not used

	real(8), intent(inout) :: p		! Pressure  (Pa)
				

	integer,  intent(in) :: zcoords, &	! 0 or 1 see below
				shear	 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q 		! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we either compute or use p

	! if shear = 1, then we use shear flow
	! if shear = 0, then we use constant u

!-----------------------------------------------------------------------
!     test case parameters
!-----------------------------------------------------------------------
        real(8), parameter :: a = 6.371229d+6                         ! [m] mean radis of the earth
	real(8), parameter :: 	X       = 500.d0,		&	! Reduced Earth reduction factor
			    	Om      = 0.d0,			&	! Rotation Rate of Earth
                            	as      = a/X,			&	! New Radius of small Earth     
			    	ueq     = 20.d0,		&	! Reference Velocity 
                            	Teq     = 300.d0,		&	! Temperature at Equator    
			    	Peq     = 100000.d0,		&	! Reference PS at Equator
                            	ztop    = 30000.d0,		&	! Model Top       
				lambdac = pi/4.d0, 		& 	! Lon of Schar Mountain Center
				phic    = 0.d0,	 		& 	! Lat of Schar Mountain Center
				h0      = 250.d0, 		& 	! Height of Mountain
				d       = 5000.d0, 		& 	! Mountain Half-Width
				xi      = 4000.d0, 		& 	! Mountain Wavelength
				cs      = 0.00025d0 		 	! Wind Shear (shear=1)
                            
      real(8) :: height							! Model level heights
      real(8) :: sin_tmp, cos_tmp					! Calculation of great circle distance
      real(8) :: r							! Great circle distance
      real(8) :: zs							! Surface height
      real(8) :: c							! Shear

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------

	sin_tmp = sin(lat) * sin(phic)
	cos_tmp = cos(lat) * cos(phic)
	
	! great circle distance with 'a/X'  

	r  = as * ACOS (sin_tmp + cos_tmp*cos(lon-lambdac))     
	zs   = h0 * exp(-(r**2)/(d**2))*(cos(pi*r/xi)**2)
	phis = g*zs

!-----------------------------------------------------------------------
!    SHEAR FLOW OR CONSTANT FLOW
!-----------------------------------------------------------------------

	if (shear .eq. 1) then

		c = cs

	else

		c = 0.d0

	endif

!-----------------------------------------------------------------------
!    TEMPERATURE 
!-----------------------------------------------------------------------

	t = Teq *(1.d0 - (c*ueq*ueq/(g))*(sin(lat)**2) )

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = peq*exp( -(ueq*ueq/(2.d0*Rd*Teq))*(sin(lat)**2) - phis/(Rd*t)    )

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE 
!-----------------------------------------------------------------------

	if (zcoords .eq. 1) then

		height = z
		p = peq*exp( -(ueq*ueq/(2.d0*Rd*Teq))*(sin(lat)**2) - g*height/(Rd*t)    )

	else
                if (hybrid_eta) p = hyam*p0 + hybm*ps ! compute the pressure based on the surface pressure and hybrid coefficients
		height = (Rd*t/(g))*log(peq/p) - (t*ueq*ueq/(2.d0*Teq*g))*(sin(lat)**2)

	endif

!-----------------------------------------------------------------------
!    THE VELOCITIES
!-----------------------------------------------------------------------

	! Zonal Velocity

	u = ueq * cos(lat) * sqrt( (2.d0*Teq/(t))*c*height + t/(Teq) )

	! Meridional Velocity

	v = 0.d0

	! Vertical Velocity = Vertical Pressure Velocity = 0

	w = 0.d0

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE test2_schaer_mountain

!==========
! Test 3-1
!==========
SUBROUTINE test3_gravity_wave (lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z		! Height (m)

	real(8), intent(inout) :: p		! Pressure  (Pa)
				

	integer,  intent(in) :: zcoords 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q 		! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output z
	! if zcoords = 0, then we use p

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
        real(8), parameter :: a = 6.371229d+6                         ! [m] mean radis of the earth
	real(8), parameter :: 	X       = 125.d0,		&	! Reduced Earth reduction factor
			    	Om      = 0.d0,			&	! Rotation Rate of Earth
                            	as      = a/X,			&	! New Radius of small Earth     
			    	u0      = 20.d0,		&	! Reference Velocity 
                            	Teq     = 300.d0,		&	! Temperature at Equator    
			    	Peq     = 100000.d0,		&	! Reference PS at Equator
                            	ztop    = 10000.d0,		&	! Model Top       
				lambdac = 2.d0*pi/3.d0, 	& 	! Lon of Pert Center
				d       = 5000.d0, 		& 	! Width for Pert
				phic    = 0.d0,	 		& 	! Lat of Pert Center
				delta_theta = 1.d0, 		& 	! Max Amplitude of Pert
				Lz      = 20000.d0, 		& 	! Vertical Wavelength of Pert
				N       = 0.01d0, 		& 	! Brunt-Vaisala frequency
				N2      = N*N, 		 	&	! Brunt-Vaisala frequency Squared
				bigG    = (g*g)/(N2*cp)			! Constant
                            
      real(8) :: height							! Model level height
      real(8) :: sin_tmp, cos_tmp					! Calculation of great circle distance
      real(8) :: r, s							! Shape of perturbation
      real(8) :: TS 							! Surface temperature
      real(8) :: t_mean, t_pert						! Mean and pert parts of temperature
      real(8) :: theta_pert						! Pot-temp perturbation

!-----------------------------------------------------------------------
!    THE VELOCITIES
!-----------------------------------------------------------------------

	! Zonal Velocity 

	u = u0 * cos(lat)

	! Meridional Velocity

	v = 0.d0

	! Vertical Velocity = Vertical Pressure Velocity = 0

	w = 0.d0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------

	phis = 0.d0

!-----------------------------------------------------------------------
!    SURFACE TEMPERATURE 
!-----------------------------------------------------------------------

	TS = bigG + (Teq-bigG)*exp( -(u0*N2/(4.d0*g*g))*(u0+2.d0*om*as)*(cos(2.d0*lat)-1.d0)    ) 

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = peq*exp( (u0/(4.0*bigG*Rd))*(u0+2.0*Om*as)*(cos(2.0*lat)-1.0)  ) &
		* (TS/Teq)**(cp/Rd)

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE AND MEAN TEMPERATURE
!-----------------------------------------------------------------------

	if (zcoords .eq. 1) then

		height = z
		p = ps*( (bigG/TS)*exp(-N2*height/g)+1.d0 - (bigG/TS)  )**(cp/Rd)

	else

		height = (-g/N2)*log( (TS/bigG)*( (p/ps)**(Rd/cp) - 1.d0  ) + 1.d0 )

	endif

	t_mean = bigG*(1.d0 - exp(N2*height/g))+ TS*exp(N2*height/g)

!-----------------------------------------------------------------------
!    rho (density), unperturbed using the background temperature t_mean
!-----------------------------------------------------------------------

!***********
!       change in version 3: density is now initialized with unperturbed background temperature,
!                            temperature perturbation is added afterwards
!***********
	rho = p/(Rd*t_mean)

!-----------------------------------------------------------------------
!    POTENTIAL TEMPERATURE PERTURBATION, 
!    here: converted to temperature and added to the temperature field
!    models with a prognostic potential temperature field can utilize
!    the potential temperature perturbation theta_pert directly and add it
!    to the background theta field (not included here)
!-----------------------------------------------------------------------

	sin_tmp = sin(lat) * sin(phic)
	cos_tmp = cos(lat) * cos(phic)

		! great circle distance with 'a/X' 

	r  = as * ACOS (sin_tmp + cos_tmp*cos(lon-lambdac)) 

	s = (d**2)/(d**2 + r**2)

	theta_pert = delta_theta*s*sin(2.d0*pi*height/Lz)

	t_pert = theta_pert*(p/p0)**(Rd/cp)

	t = t_mean + t_pert

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE test3_gravity_wave

  !
end module mod_dycoretest
!-------------------------------------------------------------------------------
