!-------------------------------------------------------------------------------
!
!+  Module of Dynamical Core Test initial condition
!
!-------------------------------------------------------------------------------
module mod_ideal_init
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
  !      0.01      13-06-12   Test cases in DCMIP2012 were imported
  !      0.02      14-01-28   Test case in Tomita and Satoh 2004 was imported
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
  use dcmip_initial_conditions_test_1_2_3, only: &
     test2_steady_state_mountain, &
     test2_schaer_mountain,       &
     test3_gravity_wave
  use mod_cnst, only: &
     pi    => CNST_PI,      &
     a     => CNST_ERADIUS, &
     omega => CNST_EOHM,    &
     g     => CNST_EGRAV,   &
     Rd    => CNST_RAIR,    &
     Cp    => CNST_CP,      &
     KAPPA => CNST_KAPPA,   &
     PRE00 => CNST_PRE00
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: dycore_input
  public :: tracer_input

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=ADM_NSYS), public :: DCTEST_type = '' !
  character(len=ADM_NSYS), public :: DCTEST_case = '' !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: hs_init
  private :: jbw_init
  private :: tracer_init
  private :: mountwave_init
  private :: gravwave_init
  private :: tomita_init

  private :: tomita_2004
  private :: eta_vert_coord_NW
  private :: steady_state
  private :: geo2prs
  private :: ps_estimation
  private :: perturbation
  private :: conv_vxvyvz
  private :: simpson

  private :: Sp_Unit_East
  private :: Sp_Unit_North

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  ! physical parameters configurations
  real(8), private :: Kap                    ! temporal value
  real(8), private :: d2r                    ! Degree to Radian
  real(8), private :: r2d                    ! Radian to Degree
  real(8), private :: eps = 1.D-14           ! minimum value
  real(8), private :: zero = 0.D0            ! zero

  ! for Held and Suarez
  real(8), private :: deltaT = 60.D0
  real(8), private :: deltaTh = 10.D0
  ! for Jablonowski
  real(8), private :: clat = 40.D0              ! perturbation center: latitude [deg]
  real(8), private :: clon = 20.D0              ! perturbation center: longitude [deg]
  real(8), private :: etaS = 1.D0               ! surface eta level
  real(8), private :: etaT = 0.2d0              ! threashold of vertical profile
  real(8), private :: eta0 = 0.252d0            ! threashold of vertical profile
  real(8), private :: t0 = 288.D0               ! [K]
  real(8), private :: delT = 4.8d+5             ! [K]
  real(8), private :: ganma = 0.005d0           ! [K m^-1]
  real(8), private :: u0 = 35.D0                ! [m s^-1]
  real(8), private :: uP = 1.D0                 ! [m s^-1]
  real(8), private :: p0 = 1.D+5                ! [Pa]
  real(8), private :: ps = 1.D+5                ! [Pa]
  logical, private, parameter :: message = .false.
  integer, private, parameter :: itrmax = 100       ! # of iteration maximum

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine dycore_input( &
       DIAG_var )
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

    character(len=ADM_NSYS) :: init_type   = ''
    character(len=ADM_NSYS) :: test_case   = ''
    real(8) :: eps_geo2prs = 1.D-2
    logical :: nicamcore = .true.

    namelist / DYCORETESTPARAM / &
       init_type,   &
       test_case,   &
       eps_geo2prs, &
       nicamcore

    integer :: ierr
    !---------------------------------------------------------------------------

    Kap = Rd / Cp
    d2r = pi/180.D0
    r2d = 180.D0/pi

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
    write(ADM_LOG_FID,nml=DYCORETESTPARAM)

    DCTEST_type = init_type
    DCTEST_case = test_case

    write(ADM_LOG_FID,*) '*** type: ', trim(init_type)
    select case(init_type)
!    case ('DCMIP2012-11','DCMIP2012-12','DCMIP2012-13' &
!          'DCMIP2012-200','DCMIP2012-21','DCMIP2012-22')

!       write(ADM_LOG_FID,*) '*** test case: ', trim(test_case)
!       call IDEAL_init_DCMIP2012( ADM_gall, ADM_kall, ADM_lall, init_type, DIAG_var(:,:,:,:) )

    case ('Heldsuarez')

       call hs_init ( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

    case ('Jablonowski')

       write(ADM_LOG_FID,*) '*** test case  : ', trim(test_case)
       write(ADM_LOG_FID,*) '*** eps_geo2prs= ', eps_geo2prs
       write(ADM_LOG_FID,*) '*** nicamcore  = ', nicamcore
       call jbw_init( ADM_gall, ADM_kall, ADM_lall, test_case, eps_geo2prs, nicamcore, DIAG_var(:,:,:,:) )

    case ('Traceradvection')

       write(ADM_LOG_FID,*) '*** test case: ', trim(test_case)
       call tracer_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )

    case ('Mountainwave')

       write(ADM_LOG_FID,*) '*** test case: ', trim(test_case)
       call mountwave_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )

    case ('Gravitywave')

       call gravwave_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

    case ('Tomita2004')

       call tomita_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

    case default

       write(ADM_LOG_FID,*) 'xxx Invalid init_type. STOP.'
       call ADM_proc_stop

    end select

    return
  end subroutine dycore_input

  !-----------------------------------------------------------------------------
  subroutine tracer_input( &
       TRC_var )
    use mod_adm, only: &
       ADM_gall,  &
       ADM_kall,  &
       ADM_KNONE, &
       ADM_lall
    use mod_random, only: &
       RANDOM_get
    use mod_cnst, only: &
       CNST_D2R, &
       CNST_PI
    use mod_gmtr, only: &
       GMTR_lon, &
       GMTR_lat
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(out) :: TRC_var(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)

    real(8) :: random(ADM_gall,ADM_kall,ADM_lall)
    integer :: deg

    integer :: g, k, l, nq, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    do nq = 1, TRC_VMAX
       call RANDOM_get( random(:,:,:) )

       if ( nq == 1 ) then ! vapor (dummy for thermodynamics)
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             TRC_var(g,k,l,nq) = 0.D0
          enddo
          enddo
          enddo
       elseif( nq == 2 ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             deg = nint( GMTR_lon(g,l) / CNST_D2R )
             if ( mod(deg,10) == 0 ) then
                TRC_var(g,k,l,nq) = real(ADM_kall-k+1,kind=8)
             else
                TRC_var(g,k,l,nq) = 0.D0
             endif
          enddo
          enddo
          enddo
       elseif( nq == 3 ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             deg = nint( GMTR_lat(g,l) / CNST_D2R )
             if ( mod(deg,10) == 0 ) then
                TRC_var(g,k,l,nq) = real(ADM_kall-k+1,kind=8)
             else
                TRC_var(g,k,l,nq) = 0.D0
             endif
          enddo
          enddo
          enddo
       else
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             TRC_var(g,k,l,nq) = random(g,k,l) * exp(real(nq,kind=8)/10.D0)
          enddo
          enddo
          enddo
       endif
    enddo

    return
  end subroutine tracer_input

  !-----------------------------------------------------------------------------
  subroutine hs_init( &
       ijdim,   &
       kdim,    &
       lall,    &
       DIAG_var )
    use mod_adm, only: &
       ADM_kmin,  &
       ADM_kmax
    use mod_grd, only: &
       GRD_Z,    &
       GRD_vz
    use mod_gmtr, only: &
       GMTR_lat, &
       GMTR_lon
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

    real(8), parameter :: deltaT  = 60.D0
    real(8), parameter :: deltaTh = 10.D0
    real(8), parameter :: eps_hs  = 1.D-7

    real(8) :: f, df
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

       lat = GMTR_lat(n,l)
       lon = GMTR_lon(n,l)

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

          f  = log(pre(k)/pre_sfc) / dz(k) + g / ( Rd * 0.5D0 * (tem(k)+tem_sfc) )
          df = 1.D0 / (pre(k)*dz(k))

          pre(k) = pre(k) - f / df
!          tem(k) = 300.D0 * ( pre(k)/PRE00 )**KAPPA
          tem(k) = ( 315.D0 - deltaT*sin(lat)**2 - deltaTh*log(pre(k)/PRE00)*cos(lat)**2 ) &
                 * ( pre(k)/PRE00 )**KAPPA
          tem(k) = max( 200.D0, tem(k) )

          if( abs(pre_save-pre(k)) <= eps_hs ) exit
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

             f  = log(pre(k)/pre(k-1)) / dz(k) + g / ( Rd * 0.5D0 * (tem(k)+tem(k-1)) )
             df = 1.D0 / (pre(k)*dz(k))

             pre(k) = pre(k) - f / df
!             tem(k) = 300.D0 * ( pre(k)/PRE00 )**KAPPA
             tem(k) = ( 315.D0 - deltaT*sin(lat)**2 - deltaTh*log(pre(k)/PRE00)*cos(lat)**2 ) &
                    * ( pre(k)/PRE00 )**KAPPA
             tem(k) = max( 200.D0, tem(k) )

             if( abs(pre_save-pre(k)) <= eps_hs ) exit

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
       ijdim,        &
       kdim,         &
       lall,         &
       test_case,    &
       eps_geo2prs,  &
       nicamcore,    &
       DIAG_var      )
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_kmin,      &
       ADM_kmax
    use mod_grd, only: &
       GRD_Z,          &
       GRD_ZH, &
       GRD_vz
    use mod_gmtr, only: &
       GMTR_lat, &
       GMTR_lon
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: lall
    character(len=*), intent(in)  :: test_case
    real(8), intent(in) :: eps_geo2prs
    logical,          intent(in)  :: nicamcore
    real(8), intent(out) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    ! work paramters
    real(8) :: lat, lon                 ! latitude, longitude on Icosahedral grid
    real(8) :: eta(kdim,2), geo(kdim)   ! eta & geopotential in ICO-grid field
    real(8) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)   ! zonal/meridional wind components in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)
    real(8) :: u(kdim)
    real(8) :: v(kdim)
    real(8) :: ps

    logical :: signal       ! if true, continue iteration
    logical :: pertb        ! if true, with perturbation
    logical :: psgm         ! if true, PS Gradient Method
    logical :: eta_limit    ! if true, value of eta is limited upto 1.0
    logical :: logout       ! log output switch for Pressure Convert

    integer :: n, k, l, itr
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.D0

    eta_limit = .true.
    psgm = .false.
    logout = .true.

    select case( trim(test_case) )
    case ('1', '4-1')  ! with perturbation
       write(ADM_LOG_FID,*) "Jablonowski Initialize - case 1: with perturbation (no rebalance)"
       pertb = .true.
    case ('2', '4-2')  ! without perturbation
       write(ADM_LOG_FID,*) "Jablonowski Initialize - case 2: without perturbation (no rebalance)"
       pertb = .false.
    case ('3')  ! with perturbation (PS Distribution Method)
       write(ADM_LOG_FID,*) "Jablonowski Initialize - PS Distribution Method: with perturbation"
       write(ADM_LOG_FID,*) "### DO NOT INPUT ANY TOPOGRAPHY ###"
       pertb = .true.
       psgm = .true.
       eta_limit = .false.
    case ('4')  ! without perturbation (PS Distribution Method)
       write(ADM_LOG_FID,*) "Jablonowski Initialize - PS Distribution Method: without perturbation"
       write(ADM_LOG_FID,*) "### DO NOT INPUT ANY TOPOGRAPHY ###"
       pertb = .false.
       psgm = .true.
       eta_limit = .false.
    case default
       write(ADM_LOG_FID,*) "Unknown test_case: '"//trim(test_case)//"' specified."
       write(ADM_LOG_FID,*) "Force changed to case 1 (with perturbation)"
       pertb = .true.
    end select
    write(ADM_LOG_FID,*) " | eps for geo2prs: ", eps_geo2prs
    write(ADM_LOG_FID,*) " | nicamcore switch for geo2prs: ", nicamcore

    do l = 1, lall
    do n = 1, ijdim
       z_local(ADM_kmin-1) = GRD_vz(n,2,l,GRD_ZH)
       do k = ADM_kmin, ADM_kmax+1
          z_local(k) = GRD_vz(n,k,l,GRD_Z)
       enddo

       lat = GMTR_lat(n,l)
       lon = GMTR_lon(n,l)

       signal = .true.
       ! iteration -----
       do itr = 1, itrmax

          if( itr == 1 ) then
             eta(:,:) = 1.D-7 ! This initial value is recommended by Jablonowsky.
          else
             call eta_vert_coord_NW( kdim, itr, z_local, tmp, geo, eta_limit, eta, signal )
          endif

          call steady_state( kdim, lat, eta, wix, wiy, tmp, geo )

          if( .NOT. signal ) exit
       enddo

       if ( itr > itrmax ) then
          write(*          ,*) 'ETA ITERATION ERROR: NOT CONVERGED', n, l
          write(ADM_LOG_FID,*) 'ETA ITERATION ERROR: NOT CONVERGED', n, l
          call ADM_proc_stop
       endif

       if (psgm) then
          call ps_estimation ( kdim, lat, eta(:,1), tmp, geo, wix, ps, nicamcore )
          call geo2prs ( kdim, ps, lat, tmp, geo, wix, prs, eps_geo2prs, nicamcore, logout )
       else
          call geo2prs ( kdim, p0, lat, tmp, geo, wix, prs, eps_geo2prs, nicamcore, logout )
       endif
       logout = .false.

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

    if (pertb) call perturbation( ijdim, kdim, lall, 5, DIAG_var(:,:,:,1:5) )

    write (ADM_LOG_FID,*) " |            Vertical Coordinate used in JBW initialization              |"
    write (ADM_LOG_FID,*) " |------------------------------------------------------------------------|"
    do k=1, kdim
       write (ADM_LOG_FID,'(3X,"(k=",I3,") HGT:",F8.2," [m]",2X,"PRS: ",F9.2," [Pa]",2X,"GH: ",F8.2," [m]",2X,"ETA: ",F9.5)') &
       k, z_local(k), prs(k), geo(k)/g, eta(k,1)
    enddo
    write (ADM_LOG_FID,*) " |------------------------------------------------------------------------|"

    return
  end subroutine jbw_init

  !-----------------------------------------------------------------------------
  subroutine tracer_init( &
       ijdim,      &
       kdim,       &
       lall,       &
       test_case,  &
       DIAG_var    )
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_kmin,      &
       ADM_kmax
    use mod_cnst, only: &
       UNDEF => CNST_UNDEF
    use mod_grd, only: &
       GRD_gz, &
       GRD_Z,  &
       GRD_ZH, &
       GRD_vz
    use mod_gmtr, only: &
       GMTR_lon, &
       GMTR_lat
    use mod_runconf, only: &
       TRC_vmax, &
       NCHEM_STR
    use mod_chemvar, only: &
       chemvar_getid
    use dcmip_initial_conditions_test_1_2_3, only: &
       test1_advection_deformation, &
       test1_advection_hadley,      &
       test1_advection_orography
    implicit none

    integer,                 intent(in)    :: ijdim
    integer,                 intent(in)    :: kdim
    integer,                 intent(in)    :: lall
    character(len=ADM_NSYS), intent(in)    :: test_case
    real(8),                 intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(8) :: lon     ! longitude            [rad]
    real(8) :: lat     ! latitude             [rad]
    real(8) :: z(kdim) ! Height               [m]
    real(8) :: p(kdim) ! pressure             [Pa]
    real(8) :: u(kdim) ! zonal      wind      [m/s]
    real(8) :: v(kdim) ! meridional wind      [m/s]
    real(8) :: w(kdim) ! vertical   wind      [m/s]
    real(8) :: t(kdim) ! temperature          [K]
    real(8) :: phis    ! surface geopotential [m2/s2], not in use
    real(8) :: ps      ! surface pressure     [Pa]   , not in use
    real(8) :: rho     ! density              [kg/m3], not in use
    real(8) :: q       ! specific humidity    [kg/kg], not in use

    real(8) :: q1(kdim) ! passive tracer       [kg/kg]
    real(8) :: q2(kdim) ! passive tracer       [kg/kg]
    real(8) :: q3(kdim) ! passive tracer       [kg/kg]
    real(8) :: q4(kdim) ! passive tracer       [kg/kg]

    real(8) :: vx(kdim)
    real(8) :: vy(kdim)
    real(8) :: vz(kdim)

    logical, parameter :: hybrid_eta = .false. ! dont use hybrid sigma-p (eta) coordinate
    integer, parameter :: zcoords    = 1       ! if zcoords = 1, then we use z and output p
    integer, parameter :: cfv        = 2       ! if cfv = 2 then our velocities follow Gal-Chen coordinates and we need to specify w
    real(8)            :: hyam       = UNDEF   ! dont use hybrid sigma-p (eta) coordinate
    real(8)            :: hybm       = UNDEF   ! dont use hybrid sigma-p (eta) coordinate
    real(8)            :: gc                   ! bar{z} for Gal-Chen coordinate

    integer :: I_pasv1, I_pasv2
    integer :: I_pasv3, I_pasv4
    integer :: n, k, l
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.D0

    I_pasv1 = 6 + NCHEM_STR + chemvar_getid( "passive001" ) - 1
    I_pasv2 = 6 + NCHEM_STR + chemvar_getid( "passive002" ) - 1
    I_pasv3 = 6 + NCHEM_STR + chemvar_getid( "passive003" ) - 1
    I_pasv4 = 6 + NCHEM_STR + chemvar_getid( "passive004" ) - 1

    select case(test_case)
    case ('1', '1-1') ! DCMIP 2012 Test 1-1: 3D Deformational Flow

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.D0

          lat = GMTR_lat(n,l)
          lon = GMTR_lon(n,l)

          do k=1, kdim
             call test1_advection_deformation( lon,     & ! [IN]
                                               lat,     & ! [IN]
                                               p(k),    & ! [INOUT]
                                               z(k),    & ! [IN]
                                               zcoords, & ! [IN]
                                               u(k),    & ! [OUT]
                                               v(k),    & ! [OUT]
                                               w(k),    & ! [OUT]
                                               t(k),    & ! [OUT]
                                               phis,    & ! [OUT]
                                               ps,      & ! [OUT]
                                               rho,     & ! [OUT]
                                               q,       & ! [OUT]
                                               q1(k),   & ! [OUT]
                                               q2(k),   & ! [OUT]
                                               q3(k),   & ! [OUT]
                                               q4(k)    ) ! [OUT]
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k=1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q1(k)
             DIAG_var(n,k,l,I_pasv2) = q2(k)
             DIAG_var(n,k,l,I_pasv3) = q3(k)
             DIAG_var(n,k,l,I_pasv4) = q4(k)
          enddo
       enddo
       enddo

    case ('2', '1-2') ! DCMIP 2012 Test 1-2: Hadley-like Meridional Circulation

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.D0

          lat = GMTR_lat(n,l)
          lon = GMTR_lon(n,l)

          do k=1, kdim
             call test1_advection_hadley( lon,     & ! [IN]
                                          lat,     & ! [IN]
                                          p(k),    & ! [INOUT]
                                          z(k),    & ! [IN]
                                          zcoords, & ! [IN]
                                          u(k),    & ! [OUT]
                                          v(k),    & ! [OUT]
                                          w(k),    & ! [OUT]
                                          t(k),    & ! [OUT]
                                          phis,    & ! [OUT]
                                          ps,      & ! [OUT]
                                          rho,     & ! [OUT]
                                          q,       & ! [OUT]
                                          q1(k)    ) ! [OUT]

             q2(k) = 0.D0
             q3(k) = 0.D0
             q4(k) = 0.D0
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k=1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q1(k)
             DIAG_var(n,k,l,I_pasv2) = q2(k)
             DIAG_var(n,k,l,I_pasv3) = q3(k)
             DIAG_var(n,k,l,I_pasv4) = q4(k)
          enddo
       enddo
       enddo

    case ('3', '1-3') ! DCMIP 2012 Test 1-3: Horizontal advection of thin cloud-like tracers in the presence of orography

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.D0

          lat = GMTR_lat(n,l)
          lon = GMTR_lon(n,l)

          do k=1, kdim
             gc = GRD_gz(k)

             call test1_advection_orography(lon,        & ! [IN]
                                            lat,        & ! [IN]
                                            p(k),       & ! [INOUT]
                                            z(k),       & ! [IN]
                                            zcoords,    & ! [IN]
                                            cfv,        & ! [IN]
                                            hybrid_eta, & ! [IN]
                                            hyam,       & ! [IN]
                                            hybm,       & ! [IN]
                                            gc,         & ! [IN]
                                            u(k),       & ! [OUT]
                                            v(k),       & ! [OUT]
                                            w(k),       & ! [OUT]
                                            t(k),       & ! [OUT]
                                            phis,       & ! [OUT]
                                            ps,         & ! [OUT]
                                            rho,        & ! [OUT]
                                            q,          & ! [OUT]
                                            q1(k),      & ! [OUT]
                                            q2(k),      & ! [OUT]
                                            q3(k),      & ! [OUT]
                                            q4(k)       ) ! [OUT]
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k=1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q1(k)
             DIAG_var(n,k,l,I_pasv2) = q2(k)
             DIAG_var(n,k,l,I_pasv3) = q3(k)
             DIAG_var(n,k,l,I_pasv4) = q4(k)
          enddo
       enddo
       enddo

    case default
       write(*          ,*) "Unknown test_case: ", trim(test_case)," specified. STOP"
       write(ADM_LOG_FID,*) "Unknown test_case: ", trim(test_case)," specified. STOP"
       call ADM_proc_stop
    end select

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
    real(8) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)   ! zonal/meridional wind components in ICO-grid field
    real(8) :: wiz(kdim)                ! vertical wind components in ICO-grid field
    real(8) :: q(kdim),     rho(kdim)   ! tracer and rho in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)

    integer :: I_pasv1
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

    I_pasv1 = 6 + chemvar_getid( "passive001" ) + NCHEM_STR - 1

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
          DIAG_var(n,k,l,I_pasv1) = q(k)
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
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_kmin,      &
       ADM_kmax
    use mod_grd, only: &
       GRD_Z,          &
       GRD_ZH, &
       GRD_vz
    use mod_gmtr, only: &
       GMTR_lon, &
       GMTR_lat
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: lall
    real(8), intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    integer, parameter :: zcoords = 1

    real(8) :: lon     ! longitude            [rad]
    real(8) :: lat     ! latitude             [rad]
    real(8) :: z(kdim) ! Height               [m]
    real(8) :: p(kdim) ! pressure             [Pa]
    real(8) :: u(kdim) ! zonal      wind      [m/s]
    real(8) :: v(kdim) ! meridional wind      [m/s]
    real(8) :: w(kdim) ! vertical   wind      [m/s]
    real(8) :: t(kdim) ! temperature          [K]
    real(8) :: phis    ! surface geopotential [m2/s2], not in use
    real(8) :: ps      ! surface pressure     [Pa]   , not in use
    real(8) :: rho     ! density              [kg/m3], not in use
    real(8) :: q       ! specific humidity    [kg/kg], not in use

    real(8) :: vx(kdim)
    real(8) :: vy(kdim)
    real(8) :: vz(kdim)

    integer :: n, k, l
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.D0

    do l = 1, lall
    do n = 1, ijdim
       z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
       do k = ADM_kmin, ADM_kmax+1
          z(k) = GRD_vz(n,k,l,GRD_Z)
       enddo
       p(:) = 0.D0

       lat = GMTR_lat(n,l)
       lon = GMTR_lon(n,l)

       do k=1, kdim
          call test3_gravity_wave( lon,     & ! [IN]
                                   lat,     & ! [IN]
                                   p(k),    & ! [INOUT]
                                   z(k),    & ! [IN]
                                   zcoords, & ! [IN]
                                   u(k),    & ! [OUT]
                                   v(k),    & ! [OUT]
                                   w(k),    & ! [OUT]
                                   t(k),    & ! [OUT]
                                   phis,    & ! [OUT]
                                   ps,      & ! [OUT]
                                   rho,     & ! [OUT]
                                   q        ) ! [OUT]
       enddo

       call conv_vxvyvz( kdim, lat, lon, u, v, vx, vy, vz )

       do k=1, kdim
          DIAG_var(n,k,l,1) = p (k)
          DIAG_var(n,k,l,2) = t (k)
          DIAG_var(n,k,l,3) = vx(k)
          DIAG_var(n,k,l,4) = vy(k)
          DIAG_var(n,k,l,5) = vz(k)
          DIAG_var(n,k,l,6) = w (k)
       enddo
    enddo
    enddo

    return
  end subroutine gravwave_init

  !-----------------------------------------------------------------------------
  subroutine tomita_init( &
       ijdim,      &
       kdim,       &
       lall,       &
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
    real(8), intent(out) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    ! work paramters
    real(8) :: lat, lon                 ! latitude, longitude on Icosahedral grid
    real(8) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)   ! zonal/meridional wind components in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)

    integer :: n, l, k, K0
    logical :: logout
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE
    logout = .true.

    DIAG_var(:,:,:,:) = 0.D0

    write(ADM_LOG_FID,*) "Qian98 Like Mountain Wave Exp. (Tomita and Satoh 2004)"

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

       call tomita_2004( kdim, lat, z_local, wix, wiy, tmp, prs, logout )
       logout = .false.
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

    write (ADM_LOG_FID,*) " |            Vertical Coordinate used in JBW initialization              |"
    write (ADM_LOG_FID,*) " |------------------------------------------------------------------------|"
    do k=1, kdim
       write (ADM_LOG_FID,'(3X,"(k=",I3,") HGT:",F8.2," [m]",2X,"PRS: ",F9.2," [Pa]")') &
       k, z_local(k), prs(k)
    enddo
    write (ADM_LOG_FID,*) " |------------------------------------------------------------------------|"

    return
  end subroutine tomita_init

  !-----------------------------------------------------------------------------
  ! estimation ps distribution by using topography
  subroutine tomita_2004( &
      kdim,     &  !--- IN : # of z dimension
      lat,      &  !--- IN : latitude
      z_local,  &  !--- IN : z vertical coordinate
      wix,      &  !--- INOUT : zonal wind field
      wiy,      &  !--- INOUT : meridional wind field
      tmp,      &  !--- INOUT : temperature
      prs,      &  !--- INOUT : pressure
      logout    )  !--- IN : log output switch
    use mod_adm, only :  &
       ADM_LOG_FID
    implicit none
    integer, intent(in) :: kdim
    real(8), intent(in) :: lat
    real(8), intent(in) :: z_local(kdim)
    real(8), intent(inout) :: wix(kdim)
    real(8), intent(inout) :: wiy(kdim)
    real(8), intent(inout) :: tmp(kdim)
    real(8), intent(inout) :: prs(kdim)
    logical, intent(in) :: logout

    integer :: i, k
    real(8) :: g1, g2, Gphi, Gzero, Pphi
    real(8), parameter :: N = 0.0187D0        ! Brunt-Vaisala Freq.
    real(8), parameter :: prs0 = 1.D5         ! pressure at the equator [Pa]
    real(8), parameter :: ux0 = 40.D0         ! zonal wind at the equator [ms-1]
    real(8) :: N2                              ! Square of Brunt-Vaisala Freq.
    real(8) :: work
    !-----

    if (logout) then
       write(ADM_LOG_FID, '("| == Tomita 2004 Mountain Wave Exp.  ( Z levels:", I4, ")")') kdim
       write(ADM_LOG_FID, '("| -- Brunt-Vaisala Freq.:", F20.10)') N
       write(ADM_LOG_FID, '("| -- Earth Angular Velocity:", F20.10)') omega
       write(ADM_LOG_FID, '("| -- Earth Radius:", F20.10)') a
       write(ADM_LOG_FID, '("| -- Earth Gravity Accel.:", F20.10)') g
    endif

    N2 = N**2.D0
    work = (N2*a) / (4.D0*g*Kap)

    g1 =   2.D0 * ( 3.D0 + 4.D0*cos(2.D0*zero) + cos(4.D0*zero) ) * (ux0**4.D0)   &
        +  8.D0 * ( 3.D0 + 4.D0*cos(2.D0*zero) + cos(4.D0*zero) ) * (ux0**3.D0)*a*omega   &
        +  8.D0 * ( 3.D0 + 4.D0*cos(2.D0*zero) + cos(4.D0*zero) ) * (ux0**2.D0)*(a**2.D0)*(omega**2.D0)   &
        - 16.D0 * ( 1.D0 + cos(2.D0*zero) ) * (ux0**2.D0)*a*g   &
        - 32.D0 * ( 1.D0 + cos(2.D0*zero) ) * ux0*(a**2.D0)*g*omega   &
        + 16.D0 * (a**2.D0) * (g**2.D0)
    g2 = (ux0**4.D0) + 4.D0*a*omega*(ux0**3.D0) + 4.D0*(a**2.D0)*(omega**2.D0)*(ux0**2.D0)
    Gzero = ( g1 / g2 )**work

    g1 =   2.D0 * ( 3.D0 + 4.D0*cos(2.D0*lat) + cos(4.D0*lat) ) * (ux0**4.D0)   &
        +  8.D0 * ( 3.D0 + 4.D0*cos(2.D0*lat) + cos(4.D0*lat) ) * (ux0**3.D0)*a*omega   &
        +  8.D0 * ( 3.D0 + 4.D0*cos(2.D0*lat) + cos(4.D0*lat) ) * (ux0**2.D0)*(a**2.D0)*(omega**2.D0)   &
        - 16.D0 * ( 1.D0 + cos(2.D0*lat) ) * (ux0**2.D0)*a*g   &
        - 32.D0 * ( 1.D0 + cos(2.D0*lat) ) * ux0*(a**2.D0)*g*omega   &
        + 16.D0 * (a**2.D0) * (g**2.D0)
    g2 = (ux0**4.D0) + 4.D0*a*omega*(ux0**3.D0) + 4.D0*(a**2.D0)*(omega**2.D0)*(ux0**2.D0)
    Gphi = ( g1 / g2 )**work

    Pphi = prs0 * ( Gzero / Gphi )

    do k=1, kdim
       wix(k) = ux0 * cos(lat)
       prs(k) = Pphi * exp( (-1.D0*N2*z_local(k)) / (g*Kap) )
       tmp(k) = (g * Kap * ( g - (wix(k)**2.D0)/a - 2.D0*omega*wix(k)*cos(lat) )) / (N2*Rd)
    enddo

    return
  end subroutine tomita_2004

  !-----------------------------------------------------------------------------
  ! eta vertical coordinate by Newton Method
  subroutine eta_vert_coord_NW( &
      kdim,      &  !--- IN : # of z dimension
      itr,       &  !--- IN : iteration number
      z,         &  !--- IN : z-height vertical coordinate
      tmp,       &  !--- IN : guessed temperature
      geo,       &  !--- IN : guessed geopotential
      eta_limit, &  !--- IN : eta limitation flag
      eta,       &  !--- INOUT : eta level vertical coordinate
      signal     )  !--- INOUT : iteration signal
    !
    implicit none
    integer, intent(in) :: itr
    integer, intent(in) :: kdim
    real(8), intent(in) :: z(kdim)
    real(8), intent(in) :: geo(kdim), tmp(kdim)
    real(8), intent(inout) :: eta(kdim,2)
    logical, intent(in) :: eta_limit
    logical, intent(inout) :: signal
    integer :: k
    real(8) :: diffmax, diff(kdim)
    real(8) :: F(kdim), Feta(kdim)
    !
    do k=1, kdim
       F(k) = -g*z(k) + geo(k)
       Feta(k) = -1.D0 * ( Rd/eta(k,1) ) * tmp(k)
       eta(k,2) = eta(k,1) - ( F(k)/Feta(k) )
       if (eta_limit) then                     ! [add] for PSDM (2013/12/20 R.Yoshida)
          if(eta(k,2) > 1.D0) eta(k,2) = 1.D0   ! not allow over 1.0 for eta
       endif
       if(eta(k,2) < 0.D0) eta(k,2) = 1.D-20
       diff(k) = abs( eta(k,2) - eta(k,1) )
    enddo
    !
    eta(:,1) = eta(:,2)
    diffmax = maxval(diff)
    if (message) write (ADM_LOG_FID, '("| Eta  ",I4,": -- MAX: ",F23.20,3X,"MIN: ",F23.20)') &
                 itr, maxval(eta(:,1)), minval(eta(:,1))
    if (message) write (ADM_LOG_FID, '("| Diff ",I4,": -- MAX: ",F23.20,3X,"MIN: ",F23.20)') &
                 itr, diffmax, minval(diff)
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
  ! calculation of steady state
  subroutine steady_state( &
      kdim, &  !--- IN : # of z dimension
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
    real(8), intent(in) :: lat
    real(8), intent(in) :: eta(kdim,2)
    real(8), intent(inout) :: wix(kdim)
    real(8), intent(inout) :: wiy(kdim)
    real(8), intent(inout) :: tmp(kdim)
    real(8), intent(inout) :: geo(kdim)
    real(8) :: eta_v
    real(8) :: work1, work2
    !
    ! ---------- horizontal mean
    work1 = pi/2.D0
    work2 = Rd*ganma/g
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
          geo(k) = t0*g/ganma * ( 1.D0 - eta(k,1)**work2 ) - Rd * delT *                              &
                  ( ( log(eta(k,1)/etaT) + 137.D0/60.D0 )*etaT**5.D0 - 5.D0*(etaT**4.D0)*eta(k,1)     &
                    + 5.D0*(etaT**3.D0)*(eta(k,1)**2.D0) - (10.D0/3.D0)*(etaT**2.D0)*(eta(k,1)**3.D0) &
                    + (5.D0/4.D0)*etaT*(eta(k,1)**4.D0) - (1.D0/5.D0)*(eta(k,1)**5.D0)                &
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
    work2 = 3.D0/4.D0 * ( pi*u0 / Rd )
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
      kdim,         &  !--- IN : # of z dimension
      ps,           &  !--- IN : surface pressure
      lat,          &  !--- IN : latitude
      tmp,          &  !--- IN : temperature
      geo,          &  !--- IN : geopotential height at full height
      wix,          &  !--- IN : zonal wind
      prs,          &  !--- IN : pressure
      eps_geo2prs,  &  !--- IN : eps
      nicamcore,    &  !--- IN : nicamcore switch
      logout        )  !--- IN : switch of log output
    !
    implicit none
    integer, intent(in) :: kdim
    real(8), intent(in) :: ps
    real(8), intent(in) :: lat
    real(8), intent(in) :: tmp(kdim)
    real(8), intent(in) :: geo(kdim)
    real(8), intent(in) :: wix(kdim)
    real(8), intent(inout) :: prs(kdim)
    real(8), intent(in) :: eps_geo2prs
    logical, intent(in) :: nicamcore
    logical, intent(in) :: logout

    integer :: i, k
    integer, parameter :: limit = 400
    real(8) :: dz, uave, diff
    real(8) :: f_cf(3), rho(3)
    real(8) :: pp(kdim)
    logical :: iteration = .false.
    logical :: do_iter = .true.
    !-----

    pp(1) = ps

    ! first guess (upward: trapezoid)
    do k=2, kdim
       dz = (geo(k) - geo(k-1))/g
       if (nicamcore) then
          uave = (wix(k) + wix(k-1)) * 0.5D0
          f_cf(1) = 2.D0*omega*uave*cos(lat) + (uave**2.D0)/a
       else
          f_cf(1) = 0.D0
       endif

       pp(k) = pp(k-1) * ( 1.D0 + dz*(f_cf(1) - g)/(2.D0*Rd*tmp(k-1)) ) &
                       / ( 1.D0 - dz*(f_cf(1) - g)/(2.D0*Rd*tmp(k)) )
    enddo
    prs(:) = pp(:)

    ! iteration (simpson)
    if (iteration) then
    do i=1, limit
       prs(1) = ps
       do k=3, kdim        ! upward
          pp(k) = simpson( prs(k), prs(k-1), prs(k-2), tmp(k), tmp(k-1), tmp(k-2), &
                  wix(k), wix(k-1), wix(k-2), geo(k), geo(k-2), lat, .false., nicamcore )
       enddo
       prs(:) = pp(:)
       do k=kdim-2, 1, -1  ! downward
          pp(k) = simpson( prs(k+2), prs(k+1), prs(k), tmp(k+2), tmp(k+1), tmp(k), &
                  wix(k+2), wix(k+1), wix(k), geo(k+2), geo(k), lat, .true., nicamcore )
       enddo
       prs(:) = pp(:)
       diff = pp(1) - ps

       if ( abs(diff) < eps_geo2prs ) then
          do_iter = .false.
          exit
       endif
    enddo
    else
       !write(ADM_LOG_FID,*) 'ETA ITERATION SKIPPED'
       do_iter = .false.
    endif

    if (do_iter) then
       write(ADM_LOG_FID,*) 'ETA ITERATION ERROR: NOT CONVERGED at GEO2PRS', diff
       stop
    endif

    ! fininalize
    prs(1) = ps
    pp(1) = ps
    do k=3, kdim        ! upward
       pp(k) = simpson( prs(k), prs(k-1), prs(k-2), tmp(k), tmp(k-1), tmp(k-2), &
               wix(k), wix(k-1), wix(k-2), geo(k), geo(k-2), lat, .false., nicamcore )
    enddo
    prs(:) = pp(:)

    if (logout) then
       if (iteration) then
          write(ADM_LOG_FID, *) " | diff (guess - ps) : ", diff, "[Pa]  --  itr times: ", (i-1)
       else
          write(ADM_LOG_FID, *) " | no iteration in geo2prs"
       endif
    endif
    if (message) then
       write (ADM_LOG_FID,'(A)') "| ----- Pressure (Final Guess) -----"
       do k=1, kdim
          write(ADM_LOG_FID, '("| K(",I3,") -- ",F20.13)') k, prs(k)
       enddo
    endif

    return
  end subroutine geo2prs

  !-----------------------------------------------------------------------------
  ! estimation ps distribution by using topography
  subroutine ps_estimation( &
      kdim,   &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      eta,  &  !--- IN : eta coordinate
      tmp,  &  !--- IN : temperature
      geo,  &  !--- IN : geopotential height at full height
      wix,  &  !--- IN : zonal wind speed
      ps,   &  !--- OUT : surface pressure
      nicamcore )  !--- IN : nicamcore switch
    use mod_adm, only :  &
       ADM_LOG_FID
    implicit none
    integer, intent(in) :: kdim
    real(8), intent(in) :: lat
    real(8), intent(in) :: eta(kdim)
    real(8), intent(in) :: tmp(kdim)
    real(8), intent(in) :: geo(kdim)
    real(8), intent(in) :: wix(kdim)
    real(8), intent(out) :: ps
    logical, intent(in) :: nicamcore

    integer :: k
    real(8), parameter :: lat0 = 0.691590985442682
    real(8) :: cs32ev, f1, f2
    real(8) :: eta_v, tmp0, tmp1
    real(8) :: ux1, ux2, hgt0, hgt1
    real(8) :: dz, uave
    real(8) :: f_cf(3), rho(3)
    real(8), parameter :: eta1 = 1.0D0
    !-----

    eta_v = (eta1 - eta0)*(pi*0.5D0)

    ! temperature at bottom of eta-grid
    tmp0 = t0                                                   &
           + (3.D0/4.D0 * (pi*u0/Rd))*eta1 * sin(eta_v) * (cos(eta_v))**0.5d0                &
           * ( ( -2.D0 * (sin(lat0))**6.D0 * (cos(lat0)**2.D0 + 1.D0/3.D0) + 10.D0/63.D0 )   &
                * 2.D0*u0*(cos(eta_v))**1.5d0                   &
                + ( 8.D0/5.D0 * (cos(lat0))**3.D0 * ((sin(lat0))**2.D0 + 2.D0/3.D0) - pi/4.D0 ) &
                * a*omega  )
    tmp1 = tmp(1)

    ! wind speed at bottom of eta-grid
    ux1 = (u0 * cos(eta_v)**1.5d0) * (sin(2.D0*lat0))**2.D0
    ux2 = wix(1)

    ! topography calculation (imported from mod_grd.f90)
    cs32ev = ( cos( (1.D0-0.252D0) * pi * 0.5D0 ) )**1.5D0
    f1 = 10.D0/63.D0 - 2.D0 * sin(lat)**6 * ( cos(lat)**2 + 1.D0/3.D0 )
    f2 = 1.6D0 * cos(lat)**3 * ( sin(lat)**2 + 2.D0/3.D0 ) - 0.25D0 * pi
    hgt1 = -1.D0 * u0 * cs32ev * ( f1*u0*cs32ev + f2*a*omega ) / g
    hgt0 = 0.D0

    ! ps estimation
    dz = hgt1 - hgt0
    if (nicamcore) then
       uave = (ux1 + ux2) * 0.5D0
       f_cf(1) = 2.D0*omega*uave*cos(lat) + (uave**2.D0)/a
    else
       f_cf(1) = 0.D0
    endif
    ps = p0 * ( 1.D0 + dz*(f_cf(1) - g)/(2.D0*Rd*tmp0) ) &
            / ( 1.D0 - dz*(f_cf(1) - g)/(2.D0*Rd*tmp1) )

    return
  end subroutine ps_estimation

  !-----------------------------------------------------------------------------
  ! setting perturbation
  subroutine perturbation( &
      ijdim,     &  !--- IN : # of ij dimension
      kdim,      &  !--- IN : # of z dimension
      lall,      &  !--- IN : # of region
      vmax,      &  !--- IN : # of variables
      DIAG_var   )  !--- INOUT : variables container
    use mod_grd, only: &
       GRD_x,      &
       GRD_x_pl,   &
       GRD_XDIR,   &
       GRD_YDIR,   &
       GRD_ZDIR
    use mod_misc, only: &
       MISC_get_latlon
    use mod_adm, only: &
       K0 => ADM_KNONE
    implicit none
    integer, intent(in) :: ijdim, kdim, lall, vmax
    real(8), intent(inout) :: DIAG_var(ijdim,kdim,lall,vmax)

    integer, parameter :: ID_vx  = 3
    integer, parameter :: ID_vy  = 4
    integer, parameter :: ID_vz  = 5
    integer :: n, k, l
    real(8) :: lat, lon
    real(8) :: r, rr, rbyrr, cla, clo
    real(8) :: ptb_wix(kdim), ptb_wiy(kdim)
    real(8) :: ptb_vx(kdim), ptb_vy(kdim), ptb_vz(kdim)

    cla = clat * d2r
    clo = clon * d2r

    do l = 1, lall
    do n = 1, ijdim
       call MISC_get_latlon( lat, lon,              &
                             GRD_x(n,K0,l,GRD_XDIR), &
                             GRD_x(n,K0,l,GRD_YDIR), &
                             GRD_x(n,K0,l,GRD_ZDIR)  )
       r = a * acos( sin(cla)*sin(lat) + cos(cla)*cos(lat)*cos(lon-clo) )
       rr = a / 10.D0
       rbyrr = r/rr
       do k = 1, kdim
          ptb_wix(k) = uP * exp( -1.D0*rbyrr**2.D0 )
          ptb_wiy(k) = 0.0d0
       enddo

       call conv_vxvyvz( kdim, lat, lon, ptb_wix, ptb_wiy, ptb_vx, ptb_vy, ptb_vz )
       do k = 1, kdim
          DIAG_var(n,k,l,ID_vx) = DIAG_var(n,k,l,ID_vx) + ptb_vx(k)
          DIAG_var(n,k,l,ID_vy) = DIAG_var(n,k,l,ID_vy) + ptb_vy(k)
          DIAG_var(n,k,l,ID_vz) = DIAG_var(n,k,l,ID_vz) + ptb_vz(k)
       enddo
    enddo
    enddo
    !
    return
  end subroutine perturbation

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
    real(8), intent(in)    :: lat
    real(8), intent(in)    :: lon
    real(8), intent(in)    :: wix(kdim)
    real(8), intent(in)    :: wiy(kdim)
    real(8), intent(inout) :: vx1d(kdim)
    real(8), intent(inout) :: vy1d(kdim)
    real(8), intent(inout) :: vz1d(kdim)
    !
    integer :: k
    real(8) :: unit_east(3), unit_north(3)
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
  function simpson( &
      pin1,        &  !--- IN : pressure (top)
      pin2,        &  !--- IN : pressure (middle)
      pin3,        &  !--- IN : pressure (bottom)
      t1,          &  !--- IN : temperature (top)
      t2,          &  !--- IN : temperature (middle)
      t3,          &  !--- IN : temperature (bottom)
      u1,          &  !--- IN : zonal wind (top)
      u2,          &  !--- IN : zonal wind (middle)
      u3,          &  !--- IN : zonal wind (bottom)
      geo1,        &  !--- IN : geopotential (top)
      geo3,        &  !--- IN : geopotential (bottom)
      lat,         &  !--- IN : latitude
      downward,    &  !--- IN : downward switch
      nicamcore )  &  !--- IN : nicamcore switch
      result (pout)
    !
    implicit none
    real(8), intent(in) :: pin1, pin2, pin3
    real(8), intent(in) :: t1, t2, t3
    real(8), intent(in) :: u1, u2, u3
    real(8), intent(in) :: geo1, geo3, lat
    logical, intent(in) :: downward
    logical, intent(in) :: nicamcore
    !
    real(8) :: dz, pout
    real(8) :: f_cf(3), rho(3)
    !---------------------------------------------------------------------------

    dz = (geo1-geo3) / g * 0.5D0

    if (nicamcore) then
       f_cf(1) = 2.D0*omega*u1*cos(lat) + (u1**2.D0)/a
       f_cf(2) = 2.D0*omega*u2*cos(lat) + (u2**2.D0)/a
       f_cf(3) = 2.D0*omega*u3*cos(lat) + (u3**2.D0)/a
    else
       f_cf(:) = 0.D0
    endif
    rho(1) = pin1 / ( Rd*t1 )
    rho(2) = pin2 / ( Rd*t2 )
    rho(3) = pin3 / ( Rd*t3 )

    if (downward) then
       pout = pin1 - ( (1.D0/3.D0) * rho(1) * ( f_cf(1) - g ) &
                       + (4.D0/3.D0) * rho(2) * ( f_cf(2) - g ) &
                     + (1.D0/3.D0) * rho(3) * ( f_cf(3) - g ) ) * dz
    else
       pout = pin3 + ( (1.D0/3.D0) * rho(1) * ( f_cf(1) - g ) &
                       + (4.D0/3.D0) * rho(2) * ( f_cf(2) - g ) &
                     + (1.D0/3.D0) * rho(3) * ( f_cf(3) - g ) ) * dz
    endif

    return
  end function simpson

  !-----------------------------------------------------------------------------
  function Sp_Unit_East( lon ) result( unit_east )
    implicit none

    real(8), intent(in) :: lon ! [rad]
    real(8)             :: unit_east(3)
    !---------------------------------------------------------------------------

    unit_east(1) = -sin(lon) ! x-direction
    unit_east(2) =  cos(lon) ! y-direction
    unit_east(3) = 0.D0      ! z-direction

    return
  end function Sp_Unit_East

  !-----------------------------------------------------------------------------
  function Sp_Unit_North( lon, lat ) result( unit_north )
    implicit none

    real(8), intent(in) :: lon, lat ! [rad]
    real(8)             :: unit_north(3)
    !---------------------------------------------------------------------------

    unit_north(1) = -sin(lat) * cos(lon) ! x-direction
    unit_north(2) = -sin(lat) * sin(lon) ! y-direction
    unit_north(3) =  cos(lat)            ! z-direction

    return
  end function Sp_Unit_North

end module mod_ideal_init
