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
  !      0.01      13-06-12   Test cases in DCMIP2012 were imported
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
  use mod_dcmip, only: &
     test1_advection_deformation,  &
     test1_advection_hadley,       &
     test2_steady_state_mountain,  &
     test2_schaer_mountain,        &
     test3_gravity_wave
  use mod_cnst, only: &
     Rd    => CNST_RAIR,    &  ! ideal gas constant of dry air
     Cp    => CNST_CP,      &  ! heat capacity of const. pressure
     g     => CNST_EGRAV,   &  ! gravity accelaration [ms^-2]
     a     => CNST_ERADIUS, &  ! mean radis of the earth [m]
     omega => CNST_EOHM,    &  ! rotation of the earth [s^-1]
     pi    => CNST_PI
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! physical parameters configurations
  real(8), private, save :: Kap                    ! temporal value
  real(8), private, save :: d2r                    ! Degree to Radian
  real(8), private, save :: r2d                    ! Radian to Degree
  real(8), private, save :: eps = 1.D-14           ! minimum value
  real(8), private, save :: zero = 0.D0            ! zero

  ! test configurations
  integer, private, parameter :: PRCS_D = 8
  ! for Held and Suarez
  real(8), private, save :: deltaT = 60.D0
  real(8), private, save :: deltaTh = 10.D0
  ! for Jablonowski
  real(8), private, save :: clat = 40.D0              ! perturbation center: latitude [deg]
  real(8), private, save :: clon = 20.D0              ! perturbation center: longitude [deg]
  real(8), private, save :: etaS = 1.D0               ! surface eta level
  real(8), private, save :: etaT = 0.2d0              ! threashold of vertical profile
  real(8), private, save :: eta0 = 0.252d0            ! threashold of vertical profile
  real(8), private, save :: t0 = 288.D0               ! [K]
  real(8), private, save :: delT = 4.8d+5             ! [K]
  real(8), private, save :: ganma = 0.005d0           ! [K m^-1]
  real(8), private, save :: u0 = 35.D0                ! [m s^-1]
  real(8), private, save :: uP = 1.D0                 ! [m s^-1]
  real(8), private, save :: p0 = 1.D+5                ! [Pa]
  real(8), private, save :: ps = 1.D+5                ! [Pa]
  logical, private, parameter :: message = .false.
  integer, private, parameter :: itrmax = 100       ! # of iteration maximum

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: dycore_input
  !
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
  private :: sphere_xyz_to_lon
  private :: sphere_xyz_to_lat
  private :: eta_vert_coord_NW
  private :: steady_state
  private :: geo2prs
  private :: geost_rebalance
  private :: perturbation
  private :: conv_vxvyvz
  private :: Sp_Unit_East
  private :: Sp_Unit_North
  !
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
    character(len=ADM_NSYS) :: test_case

    namelist / DYCORETESTPARAM / &
       init_type, &
       test_case

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
    write(ADM_LOG_FID,DYCORETESTPARAM)

    write(ADM_LOG_FID,*) '*** type: ', trim(init_type)

    if (      init_type == "Jablonowski"     &
         .OR. init_type == "Traceradvection" &
         .OR. init_type == "Mountainwave"    ) then

       write(ADM_LOG_FID,*) '*** test case: ', trim(test_case)

    endif

    if ( init_type == "Heldsuarez" ) then

       call hs_init ( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

    elseif( init_type == "Jablonowski" ) then

       call jbw_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )

    elseif( init_type == "Traceradvection" ) then

       call tracer_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )

    elseif( init_type == "Mountainwave" ) then

       call mountwave_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )

    elseif( init_type == "Gravitywave" ) then

       call gravwave_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

    elseif( init_type == "Tomita2004" ) then

       call tomita_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

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
       ADM_KNONE, &
       ADM_kmin,  &
       ADM_kmax
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV, &
       RAIR  => CNST_RAIR,  &
       KAPPA => CNST_KAPPA, &
       PRE00 => CNST_PRE00
    use mod_grd, only: &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_Z,    &
       GRD_x,    &
       GRD_vz
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
    real(8) :: eps_hs = 1.0d-7

    integer :: n, k, l, itr
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.D0

    do l = 1, lall
    do n = 1, ijdim

       dz(ADM_kmin) = GRD_vz(n,ADM_kmin,l,GRD_Z)
       do k = ADM_kmin+1, ADM_kmax+1
          dz(k) = GRD_vz(n,k,l,GRD_Z) - GRD_vz(n,k-1,l,GRD_Z)
       enddo

       call MISC_get_latlon( lat, lon,                      &
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

             f  = log(pre(k)/pre(k-1)) / dz(k) + GRAV / ( RAIR * 0.5D0 * (tem(k)+tem(k-1)) )
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
    real(PRCS_D) :: wix(kdim),   wiy(kdim)   ! zonal/meridional wind components in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)
    real(8) :: u(kdim)
    real(8) :: v(kdim)

    logical :: signal       ! if true, continue iteration
    logical :: pertb        ! if true, with perturbation
    logical :: psgm         ! if true, PS Gradient Method
    logical :: eta_limit    ! if true, value of eta is limited upto 1.0
    logical :: logout       ! log output switch for Pressure Convert
    integer :: n, l, k, itr, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

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
             call eta_vert_coord_NW( kdim, itr, z_local, tmp, geo, eta_limit, eta, signal )
          endif

          call steady_state( kdim, lat, eta, wix, wiy, tmp, geo )

          if( .NOT. signal ) exit
       enddo

       if ( itr > itrmax ) then
          write(ADM_LOG_FID,*) 'ETA ITERATION ERROR: NOT CONVERGED', n, l
          stop
       endif

       if (prsph) then
          call prs_phi_profile( kdim, nlat, eta, latphi, prsphi )
          prsph = .false.
       endif

       call geo2prs ( kdim, lat, tmp, geo, wix, prs, logout )
       logout = .false.

       if(psgm) call ps_estimation ( kdim, lat, eta(:,1), tmp, geo, wix, prs )

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
  !
  !-----------------------------------------------------------------------------
  subroutine tracer_init( &
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
    real(8) :: lat, lon                  ! latitude, longitude on Icosahedral grid
    real(8) :: prs(kdim),   tmp(kdim)    ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)    ! zonal/meridional wind components in ICO-grid field
    real(8) :: wiz(kdim)                 ! vertical wind components in ICO-grid field
    real(8) :: q(kdim),     rho(kdim)    ! Qvapor and rho in ICO-grid field
    real(8) :: q1(kdim), q2(kdim)        ! passive tracer in ICO-grid field
    real(8) :: q3(kdim), q4(kdim)        ! passive tracer in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)

    integer :: I_passive1, I_passive2
    integer :: I_passive3, I_passive4
    integer :: n, l, k, K0

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

    I_passive1 = 6 + NCHEM_STR + chemvar_getid( "passive1" ) - 1
    I_passive2 = 6 + NCHEM_STR + chemvar_getid( "passive2" ) - 1
    I_passive3 = 6 + NCHEM_STR + chemvar_getid( "passive3" ) - 1
    I_passive4 = 6 + NCHEM_STR + chemvar_getid( "passive4" ) - 1

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
       ! DCMIP: TEST CASE 11 - Pure Advection - 3D deformational flow
       case ('1', '1-1')
          do k=1, kdim
             call test1_advection_deformation (lon, lat, prs(k), z_local(k), zcoords,    &
                                                wix(k), wiy(k), wiz(k), tmp(k), phis, ps, &
                                                rho(k), q(k), q1(k), q2(k), q3(k), q4(k)  )
          enddo
       ! DCMIP: TEST CASE 12 - Pure Advection - 3D HADLEY-like flow
       case ('2', '1-2')
          do k=1, kdim
             call test1_advection_hadley (lon, lat, prs(k), z_local(k), zcoords,    &
                                           wix(k), wiy(k), wiz(k), tmp(k), phis, ps, &
                                           rho(k), q(k), q1(k)  )
             q2(k) = 0.0d0
             q3(k) = 0.0d0
             q4(k) = 0.0d0
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
          !DIAG_var(n,k,l,XXXXXXXXX) = q(k)  ! this has not prepared yet (20130612)
          DIAG_var(n,k,l,I_passive1) = q1(k)
          DIAG_var(n,k,l,I_passive2) = q2(k)
          DIAG_var(n,k,l,I_passive3) = q3(k)
          DIAG_var(n,k,l,I_passive4) = q4(k)
       enddo
    enddo
    enddo

    return
  end subroutine tracer_init
  !-----------------------------------------------------------------------------
  !
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
    real(8) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)   ! zonal/meridional wind components in ICO-grid field
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
  !
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
    real(8) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(8) :: wix(kdim),   wiy(kdim)   ! zonal/meridional wind components in ICO-grid field
    real(8) :: wiz(kdim)                ! vertical wind components in ICO-grid field
    real(8) :: q(kdim),     rho(kdim)   ! tracer and rho in ICO-grid field

    real(8) :: z_local (kdim)
    real(8) :: vx_local(kdim)
    real(8) :: vy_local(kdim)
    real(8) :: vz_local(kdim)

    integer :: I_passive1
    integer :: n, l, k, K0

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
  !
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
    real(PRCS_D) :: lat, lon                 ! latitude, longitude on Icosahedral grid
    real(PRCS_D) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(PRCS_D) :: wix(kdim),   wiy(kdim)   ! zonal/meridional wind components in ICO-grid field

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
  !
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
    real(PRCS_D), intent(in) :: z(kdim)
    real(PRCS_D), intent(in) :: geo(kdim), tmp(kdim)
    real(PRCS_D), intent(inout) :: eta(kdim,2)
    logical, intent(in) :: eta_limit
    logical, intent(inout) :: signal
    integer :: k
    real(PRCS_D) :: diffmax, diff(kdim)
    real(PRCS_D) :: F(kdim), Feta(kdim)
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
  !
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
  !
  !-----------------------------------------------------------------------------
  ! convert geopotential height to pressure
  subroutine geo2prs( &
      kdim,   &  !--- IN : # of z dimension
      lat,    &  !--- IN : latitude
      tmp,    &  !--- IN : temperature
      geo,    &  !--- IN : geopotential height at full height
      wix,    &  !--- IN : zonal wind
      prs,    &  !--- INOUT : pressure
      logout  )  !--- IN : switch of log output
    !
    implicit none
    integer :: k
    integer, intent(in) :: kdim
    real(PRCS_D), intent(in) :: lat
    real(PRCS_D), intent(in) :: tmp(kdim)
    real(PRCS_D), intent(in) :: geo(kdim)
    real(PRCS_D), intent(in) :: wix(kdim)
    real(PRCS_D), intent(inout) :: prs(kdim)
    logical, intent(in) :: logout

    real(PRCS_D) :: dz, uave
    real(PRCS_D) :: f_cf(3), rho(3)
    logical :: nicamcore = .true.
    !-----

    prs(1) = p0

    ! first guess upper bottom level
    do k=2, kdim
       dz = (geo(k) - geo(k-1))/g
       if (nicamcore) then
          uave = (wix(k) + wix(k-1)) * 0.5D0
          f_cf(1) = 2.D0*omega*uave*cos(lat) + (uave**2.D0)/a
       else
          f_cf(1) = 0.D0
       endif

       prs(k) = prs(k-1) * ( 1.D0 + dz*(f_cf(1) - g)/(2.D0*Rd*tmp(k-1)) ) &
                         / ( 1.D0 - dz*(f_cf(1) - g)/(2.D0*Rd*tmp(k)) )
    enddo

    ! final guess
    do k=kdim-2, 1
       dz = (geo(k+2) - geo(k))/g
       if (nicamcore) then
          uave = (wix(k) + wix(k-1)) * 0.5D0
          f_cf(1) = 2.D0*omega*wix(k+2)*cos(lat) + (wix(k+2)**2.D0)/a
          f_cf(2) = 2.D0*omega*wix(k+1)*cos(lat) + (wix(k+1)**2.D0)/a
          f_cf(3) = 2.D0*omega*wix( k )*cos(lat) + (wix( k )**2.D0)/a
       else
          f_cf(:) = 0.D0
       endif
       rho(1) = prs(k+2) / ( Rd*tmp(k+2) )
       rho(2) = prs(k+1) / ( Rd*tmp(k+1) )
       rho(3) = prs( k ) / ( Rd*tmp( k) )

       prs(k) = prs(k+1) - ( &
                               (1.D0/3.D0) * rho(1) * ( f_cf(1) - g ) &
                             + (4.D0/3.D0) * rho(2) * ( f_cf(2) - g ) &
                             + (1.D0/3.D0) * rho(3) * ( f_cf(3) - g ) &
                           ) * dz
    enddo

    if (logout) write(ADM_LOG_FID, *) " | diff(guess - p0): ", (prs(1) - p0)
    if (message) then
       write (ADM_LOG_FID,'(A)') "| ----- Pressure (Final Guess) -----"
       do k=1, kdim
          write(ADM_LOG_FID, '("| K(",I3,") -- ",F20.13)') k, prs(k)
       enddo
    endif

    return
  end subroutine geo2prs
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  ! estimation ps distribution by using topography
  subroutine ps_estimation( &
      kdim,   &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      eta,  &  !--- IN : eta coordinate
      tmp,  &  !--- IN : temperature
      geo,  &  !--- IN : geopotential height at full height
      wix,  &  !--- IN : zonal wind speed
      prs   )  !--- INOUT : pressure
    use mod_adm, only :  &
       ADM_LOG_FID
    implicit none
    integer, intent(in) :: kdim
    real(PRCS_D), intent(in) :: lat
    real(PRCS_D), intent(in) :: eta(kdim)
    real(PRCS_D), intent(in) :: tmp(kdim)
    real(PRCS_D), intent(in) :: geo(kdim)
    real(PRCS_D), intent(in) :: wix(kdim)
    real(PRCS_D), intent(inout) :: prs(kdim)

    integer :: k
    real(PRCS_D), parameter :: lat0 = 0.691590985442682
    real(PRCS_D) :: cs32ev, f1, f2
    real(PRCS_D) :: eta_v, tmp0, tmp1, ux1, ux2, hgt0, hgt1
    real(PRCS_D) :: dz, uave
    real(PRCS_D) :: f_cf(3), rho(3)
    logical :: nicamcore = .true.
    !-----

    ! temperature at bottom of eta-grid
    tmp0 = t0
    tmp1 = tmp(1)

    ! wind speed at bottom of eta-grid
    eta_v = (eta(1) - eta0)*(pi*0.5D0)
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
    prs(1) = p0 * ( 1.D0 + dz*(f_cf(1) - g)/(2.D0*Rd*tmp0) ) &
                / ( 1.D0 - dz*(f_cf(1) - g)/(2.D0*Rd*tmp1) )

    ! first guess upper bottom level
    do k=2, kdim
       dz = (geo(k) - geo(k-1))/g
       if (nicamcore) then
          uave = (wix(k) + wix(k-1)) * 0.5D0
          f_cf(1) = 2.D0*omega*uave*cos(lat) + (uave**2.D0)/a
       else
          f_cf(1) = 0.D0
       endif

       prs(k) = prs(k-1) * ( 1.D0 + dz*(f_cf(1) - g)/(2.D0*Rd*tmp(k-1)) ) &
                         / ( 1.D0 - dz*(f_cf(1) - g)/(2.D0*Rd*tmp(k)) )
    enddo

    ! final guess
    do k=kdim-2, 1
       dz = (geo(k+2) - geo(k))/g
       if (nicamcore) then
          uave = (wix(k) + wix(k-1)) * 0.5D0
          f_cf(1) = 2.D0*omega*wix(k+2)*cos(lat) + (wix(k+2)**2.D0)/a
          f_cf(2) = 2.D0*omega*wix(k+1)*cos(lat) + (wix(k+1)**2.D0)/a
          f_cf(3) = 2.D0*omega*wix( k )*cos(lat) + (wix( k )**2.D0)/a
       else
          f_cf(:) = 0.D0
       endif
       rho(1) = prs(k+2) / ( Rd*tmp(k+2) )
       rho(2) = prs(k+1) / ( Rd*tmp(k+1) )
       rho(3) = prs( k ) / ( Rd*tmp( k) )

       prs(k) = prs(k+1) - ( &
                               (1.D0/3.D0) * rho(1) * ( f_cf(1) - g ) &
                             + (4.D0/3.D0) * rho(2) * ( f_cf(2) - g ) &
                             + (1.D0/3.D0) * rho(3) * ( f_cf(3) - g ) &
                           ) * dz
    enddo

    if (message) then
       write (*,'(A)') "| ----- Pressure (Final Guess) -----"
       do k=1, kdim
          write(ADM_LOG_FID, '("| K(",I3,") -- ",F20.13)') k, prs(k)
       enddo
    endif

    return
  end subroutine ps_estimation
  !-----------------------------------------------------------------------------
  !
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
    real(PRCS_D), intent(in) :: lat
    real(PRCS_D), intent(in) :: z_local(kdim)
    real(PRCS_D), intent(inout) :: wix(kdim)
    real(PRCS_D), intent(inout) :: wiy(kdim)
    real(PRCS_D), intent(inout) :: tmp(kdim)
    real(PRCS_D), intent(inout) :: prs(kdim)
    logical, intent(in) :: logout

    integer :: i, k
    real(PRCS_D) :: g1, g2, Gphi, Gzero, Pphi
    real(PRCS_D), parameter :: N = 0.0187D0        ! Brunt-Vaisala Freq.
    real(PRCS_D), parameter :: prs0 = 1.D5         ! pressure at the equator [Pa]
    real(PRCS_D), parameter :: ux0 = 40.D0         ! zonal wind at the equator [ms-1]
    real(PRCS_D) :: N2                              ! Square of Brunt-Vaisala Freq.
    real(PRCS_D) :: work
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
  !
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

    integer, parameter :: ID_prs = 1
    integer, parameter :: ID_tmp = 2
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
  !
  !  "subroutine surface_height" was moved to share/mod_grd.f90
  !-----------------------------------------------------------------------------
  !
end module mod_dycoretest
!-------------------------------------------------------------------------------
