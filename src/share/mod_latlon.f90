!-------------------------------------------------------------------------------
!>
!! Geodesic (Lat-Lon) grid module
!!
!! @par Description
!!         This module contains the tools to convert between icosaheral grid
!!         and lat-lon grid
!!
!! @author S.Iga
!!
!! @par History
!! @li      2004-02-17 (S.Iga)    Imported from igdc-4.39
!! @li      2004-03-05 (S.Iga)    'mod_latlon2.f90' is merged into this module.
!! @li      2004-05-31 (H.Tomita) Delete debug write statements
!! @li      2005-11-10 (M.Satoh)  bug fix: output_lldata_type_in
!! @li      2005-12-17 (M.Satoh)  add namelist options for lat/lon max/min_deg
!! @li      2006-02-10 (S.Iga)    bug fix: for the case LL grid is near to
!! @li                            ICO grid (in the past, for gl11, weight at
!! @li                            ix=8197,iy=4176 was NaN)
!! @li      2007-07-12 (T.Mitsui) bug fix: "fid" had been undefined in mkllmap.
!! @li      2009-07-17 (Y.Yamada) bug fix: negative area had existed in mkllmap.
!! @li      2011-01-11 (S.Iga)    handling "lon>180"
!! @li      2011-11-09  H.Yashiro [mod] Avoid arc-cos, precise calculation
!!
!<
module mod_latlon
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use mod_precision
  use mod_stdio
  use mod_debug
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LATLON_ico_setup
  public :: LATLON_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: GMTR_P_nmax_var = 2
  integer, public, parameter :: GMTR_P_LAT = 1
  integer, public, parameter :: GMTR_P_LON = 2

  real(RP), public, allocatable :: GMTR_P_ll   (:,:,:,:)
  real(RP), public, allocatable :: GMTR_P_ll_pl(:,:,:,:)

  character(len=H_SHORT),  public :: polygon_type = 'ON_SPHERE' ! triangle is fit to the sphere
  !                                                 'ON_PLANE'  ! triangle is treated as 2D

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: setup_latlon
  private :: set_equidist_grid
  private :: set_gaussian_grid
  private :: mkrelmap_ico2ll

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: latlon_type = 'EQUIDIST' ! grid type ( equidist or gaussian )
  integer,                private :: imax        = 360        ! number of longitude
  integer,                private :: jmax        = 180        ! number of latitude
  real(RP),               private :: lonmin      = -999.0_RP  ! minimun longitude of region window in deg
  real(RP),               private :: lonmax      = -999.0_RP  ! maximun longitude of region window in deg
  real(RP),               private :: latmin      = -999.0_RP  ! minimun latitude of region window in deg
  real(RP),               private :: latmax      = -999.0_RP  ! maximun latitude of region window in deg
  logical,                private :: lon_offset  = .true.     ! logitude offset
  real(RP),               private :: polar_limit = -999.0_RP  ! search all longitude if abs(lat) > polar_limit

  character(len=H_LONG),  private :: SAMPLE_OUT_BASENAME = ''
  character(len=H_SHORT), private :: SAMPLE_io_mode      = 'ADVANCED'

  real(RP), private, allocatable :: lat(:)
  real(RP), private, allocatable :: lon(:)

  integer,  private              :: nmax_llgrid
  integer,  private, allocatable :: nmax_llgrid_rgn(:)

  integer,  private, allocatable :: lon_index(:)
  integer,  private, allocatable :: lat_index(:)
  integer,  private, allocatable :: l_index  (:)
  integer,  private, allocatable :: t_index  (:)
  integer,  private, allocatable :: n1_index (:)
  integer,  private, allocatable :: n2_index (:)
  integer,  private, allocatable :: n3_index (:)
  real(RP), private, allocatable :: w1       (:)
  real(RP), private, allocatable :: w2       (:)
  real(RP), private, allocatable :: w3       (:)

  real(4),  private, allocatable :: checkmap   (:,:)
  real(4),  private, allocatable :: checkmapsum(:,:)

  logical,  private :: debug = .false.
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> setup lat/lon value of the ico-grid (without mod_gmtr)
  subroutine LATLON_ico_setup
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_KNONE,       &
       ADM_GSLF_PL,     &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo
    use mod_vector, only: &
       VECTR_xyz2latlon
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_x,    &
       GRD_x_pl
    implicit none

    integer :: ij, n, k, l
    !---------------------------------------------------------------------------

    k = ADM_KNONE

    !--- setup point data
    allocate( GMTR_P_ll   (ADM_gall,   ADM_KNONE,ADM_lall,   GMTR_P_nmax_var) )
    allocate( GMTR_P_ll_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,GMTR_P_nmax_var) )
    GMTR_P_ll   (:,:,:,:) = 0.0_RP
    GMTR_P_ll_pl(:,:,:,:) = 0.0_RP

    do l = 1, ADM_lall
       do n = 1, ADM_IooJoo_nmax
          ij = ADM_IooJoo(n,ADM_GIoJo)
          call VECTR_xyz2latlon( GMTR_P_ll(ij,k,l,GMTR_P_LAT), &
                                 GMTR_P_ll(ij,k,l,GMTR_P_LON), &
                                 GRD_x    (ij,k,l,GRD_XDIR),   &
                                 GRD_x    (ij,k,l,GRD_YDIR),   &
                                 GRD_x    (ij,k,l,GRD_ZDIR)    )
       enddo ! ij loop
    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_GSLF_PL
       do l = 1,ADM_lall_pl
          call VECTR_xyz2latlon( GMTR_P_ll_pl(n,k,l,GMTR_P_LAT), &
                                 GMTR_P_ll_pl(n,k,l,GMTR_P_LON), &
                                 GRD_x_pl    (n,k,l,GRD_XDIR),   &
                                 GRD_x_pl    (n,k,l,GRD_YDIR),   &
                                 GRD_x_pl    (n,k,l,GRD_ZDIR)    )
       enddo ! l loop
    endif

    !--- communication of point data
    call COMM_data_transfer( GMTR_P_ll, GMTR_P_ll_pl )

    return
  end subroutine LATLON_ico_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LATLON_setup( output_dirname )
    use mod_adm, only: &
       ADM_proc_stop,      &
       ADM_COMM_WORLD,     &
       ADM_prc_all,        &
       ADM_prc_tab,        &
       ADM_prc_run_master, &
       ADM_prc_me,         &
       ADM_lall
    use mod_cnst, only: &
       CNST_PI
    implicit none

    character(len=*), intent(in) :: output_dirname

    real(RP) :: latmax_deg      =   90.0_RP
    real(RP) :: latmin_deg      =  -90.0_RP
    real(RP) :: lonmax_deg      =  180.0_RP
    real(RP) :: lonmin_deg      = -180.0_RP
    real(RP) :: polar_limit_deg =   89.0_RP

    namelist / LATLONPARAM / &
         latlon_type,         &
         imax,                &
         jmax,                &
         lonmin_deg,          &
         lonmax_deg,          &
         latmin_deg,          &
         latmax_deg,          &
         lon_offset,          &
         polar_limit_deg,     &
         SAMPLE_OUT_BASENAME, &
         SAMPLE_io_mode,      &
         debug

    character(len=H_LONG) :: fname

    real(RP) :: d2r

    integer :: globalsum
    integer :: sendbuf(1)
    integer :: recvbuf(ADM_prc_all)

    integer :: fid, ierr
    integer :: nstart, nend
    integer :: n, l, rgnid, i, j
    !---------------------------------------------------------------------------

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[latlon]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=LATLONPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** LATLONPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist LATLONPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist LATLONPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(IO_FID_LOG,nml=LATLONPARAM)

    d2r    = CNST_PI / 180.0_RP
    latmax = latmax_deg * d2r
    latmin = latmin_deg * d2r
    lonmax = lonmax_deg * d2r
    lonmin = lonmin_deg * d2r

    polar_limit = abs(polar_limit_deg) * d2r

    !--- setup latitude-longitude grid
    allocate( lat(jmax) )
    allocate( lon(imax) )

    call setup_latlon

    if( ADM_prc_me == ADM_prc_run_master ) then
       fid = IO_get_available_fid()
       open( unit   = fid,                                 &
             file   = trim(output_dirname)//'/llmap.info', &
             form   = 'unformatted',                       &
             status = 'unknown'                            )
          write(fid) imax
          write(fid) lon(:)
          write(fid) jmax
          write(fid) lat(:)
       close(fid)
    endif

    allocate( checkmap   (imax,jmax) )
    allocate( checkmapsum(imax,jmax) )
    checkmap   (:,:) = 0.0
    checkmapsum(:,:) = 0.0

    write(IO_FID_LOG,*) '====== Lat-Lon grid info. ======'
    if ( latlon_type == 'EQUIDIST' ) then
       write(IO_FID_LOG,*) '--- Latitude  type   : Equal distance'
    elseif( latlon_type == 'GAUSSIAN' ) then
       write(IO_FID_LOG,*) '--- Latitude  type   : Gaussian'
    endif
    if ( lon_offset ) then
       write(IO_FID_LOG,*) '--- Longitude offset : yes'
    else
       write(IO_FID_LOG,*) '--- Longitude offset : no'
    endif
    write(IO_FID_LOG,*)    '--- # of Latitude    :', jmax
    write(IO_FID_LOG,*)    '--- # of Longitude   :', imax
    write(IO_FID_LOG,*)    '--- Latitude  range  :', lat(1)/d2r,' - ', lat(jmax)/d2r
    write(IO_FID_LOG,*)    '--- Longitude range  :', lon(1)/d2r,' - ', lon(imax)/d2r

    allocate( nmax_llgrid_rgn(ADM_lall) )

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Start counting grid.'

    ! count lat-lon number
    call mkrelmap_ico2ll( 'GET_NUM' )

    if ( debug ) then
       call MPI_Allreduce( checkmap(1,1),      &
                           checkmapsum(1,1),   &
                           imax*jmax,          &
                           MPI_REAL,           &
                           MPI_SUM,            &
                           ADM_COMM_WORLD, &
                           ierr                )
    endif

    write(IO_FID_LOG,*) '# of managing llgrid'
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       write(IO_FID_LOG,*) 'region=', rgnid, ', llgrid=', nmax_llgrid_rgn(l)

       if ( debug ) then
          do j = 1, jmax
          do i = 1, imax
             if    ( checkmapsum(i,j) >  1.0 ) then
                write(IO_FID_LOG,*) 'dupicate! (i,j)=', i, j, checkmapsum(i,j)
             elseif( checkmapsum(i,j) == 0.0 ) then
                write(IO_FID_LOG,*) 'missed!   (i,j)=', i, j, checkmapsum(i,j)
             endif
          enddo
          enddo

          call IO_make_idstr(fname,trim(output_dirname)//'/checkmap','grd',rgnid)

          fid = IO_get_available_fid()
          open( unit   = fid,           &
                file   = trim(fname),   &
                form   = 'unformatted', &
                access = 'direct',      &
                recl   = imax*jmax*4,   &
                status = 'unknown'      )
             write(fid,rec=1) checkmap(:,:)
          close(fid)
       endif
    enddo

    ! check total lat-lon number
    sendbuf(1) = nmax_llgrid

    call MPI_Allgather( sendbuf,            &
                        1,                  &
                        MPI_INTEGER,        &
                        recvbuf,            &
                        1,                  &
                        MPI_INTEGER,        &
                        ADM_COMM_WORLD, &
                        ierr                )

    globalsum = sum( recvbuf(:) )

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) 'imax x jmax                    = ', imax*jmax
    write(IO_FID_LOG,*) 'global total of counted llgrid = ', globalsum
    if ( globalsum /= imax*jmax ) then
       write(*,          *) 'counted llgrid does not match!'
       write(IO_FID_LOG,*) 'counted llgrid does not match!'
!       call ADM_proc_stop
    endif

    allocate( lon_index(nmax_llgrid) )
    allocate( lat_index(nmax_llgrid) )
    allocate( l_index  (nmax_llgrid) )
    allocate( t_index  (nmax_llgrid) )
    allocate( n1_index (nmax_llgrid) )
    allocate( n2_index (nmax_llgrid) )
    allocate( n3_index (nmax_llgrid) )
    allocate( w1       (nmax_llgrid) )
    allocate( w2       (nmax_llgrid) )
    allocate( w3       (nmax_llgrid) )

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Start calc relation map.'

    ! calc relation map
    call mkrelmap_ico2ll( 'SET_INDEX' )

    write(IO_FID_LOG,*) '*** OK.'

    ! output relation map
    if ( debug ) then
       do l = 1, ADM_lall
          nend   = sum(nmax_llgrid_rgn(1:l))
          nstart = nend - nmax_llgrid_rgn(l) + 1

          do n = nstart, nend
             if ( abs(w1(n)+w2(n)+w3(n)-1.0_RP) > 1.E-15_RP ) then
                write(IO_FID_LOG,'(A,2I6,E30.20)') '(lat,lon,area)=', &
                lat_index(n),lon_index(n),w1(n)+w2(n)+w3(n)
             endif
          enddo
       enddo
    endif

    ! output relation map
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       call IO_make_idstr(fname,trim(output_dirname)//'/llmap','rgn',rgnid)

       fid = IO_get_available_fid()
       open( unit   = fid,           &
             file   = trim(fname),   &
             form   = 'unformatted', &
             status = 'unknown'      )

          nend   = sum(nmax_llgrid_rgn(1:l))
          nstart = nend - nmax_llgrid_rgn(l) + 1

          write(fid) nmax_llgrid_rgn(l)
          if ( nmax_llgrid_rgn(l) /= 0 ) then
             write(fid) lon_index(nstart:nend)
             write(fid) lat_index(nstart:nend)
             write(fid) n1_index (nstart:nend)
             write(fid) n2_index (nstart:nend)
             write(fid) n3_index (nstart:nend)
             write(fid) w1       (nstart:nend)
             write(fid) w2       (nstart:nend)
             write(fid) w3       (nstart:nend)
          endif
       close(fid)
    enddo

    if ( SAMPLE_OUT_BASENAME /= '' ) then
       call LL_outputsample
    endif

    ! ASCII output for debug
!    do l = 1, ADM_lall
!       rgnid = ADM_prc_tab(l,ADM_prc_me)
!
!       call IO_make_idstr(fname,trim(output_dirname)//'/llmap','rgntxt',rgnid)
!
!       fid = IO_get_available_fid()
!       open( unit   = fid,           &
!             file   = trim(fname),   &
!             form   = 'formatted', &
!             status = 'unknown'      )
!
!          nend   = sum(nmax_llgrid_rgn(1:l))
!          nstart = nend - nmax_llgrid_rgn(l) + 1
!
!          write(fid,*) nmax_llgrid_rgn(l)
!          if ( nmax_llgrid_rgn(l) /= 0 ) then
!             write(fid,*) lon_index(nstart:nend)
!             write(fid,*) lat_index(nstart:nend)
!             write(fid,*) n1_index (nstart:nend)
!             write(fid,*) n2_index (nstart:nend)
!             write(fid,*) n3_index (nstart:nend)
!             write(fid,*) w1       (nstart:nend)
!             write(fid,*) w2       (nstart:nend)
!             write(fid,*) w3       (nstart:nend)
!          endif
!       close(fid)
!    enddo

    return
  end subroutine LATLON_setup

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine mkrelmap_ico2ll
  !>
  subroutine mkrelmap_ico2ll( what_is_done )
    use mod_adm, only: &
       ADM_proc_stop,     &
       ADM_prc_tab,       &
       ADM_prc_me,        &
       ADM_rgnid_npl_mng, &
       ADM_rgnid_spl_mng, &
       ADM_TI,            &
       ADM_TJ,            &
       ADM_lall,          &
       ADM_gall_1d,       &
       ADM_gmax,          &
       ADM_gmin,          &
       ADM_KNONE,         &
       ADM_IooJoo_nmax,   &
       ADM_IooJoo,        &
       ADM_GIoJo,         &
       ADM_GIpJo,         &
       ADM_GIpJp,         &
       ADM_GIoJp
    use mod_cnst, only: &
       CNST_PI
    use mod_vector, only: &
       VECTR_triangle, &
       VECTR_cross,    &
       VECTR_dot,      &
       VECTR_abs
    use mod_grd, only: &
       GRD_rscale, &
       GRD_x
    implicit none

    character(len=*), intent(in) :: what_is_done

    real(RP), parameter :: rscale = 1.0_RP

    real(RP), parameter :: o(3) = 0.0_RP
    real(RP) :: r0(3), r1(3), r2(3), r3(3)
    real(RP) :: nvec(3), v12xv10(3), v23xv20(3), v31xv30(3)
    real(RP) :: judge12, judge23, judge31
    real(RP) :: rf, rn, ip, len
    real(RP) :: v01(3)

    real(RP) :: coslat(jmax), sinlat(jmax)
    real(RP) :: coslon(imax), sinlon(imax)
    real(RP) :: lat1, lat2, lat3
    real(RP) :: lon1, lon2, lon3
    real(RP) :: latmin_l,latmax_l
    real(RP) :: lonmin_l,lonmax_l
    logical :: near_pole

    real(RP) :: area_total, area1, area2, area3

    real(RP) :: eps_judge  = 1.E-18_RP ! marginal value for inner products
    real(RP) :: eps_latlon = 1.E-15_RP ! marginal square near grid points (in radian)
    real(RP) :: eps_vertex = 1.E-15_RP ! marginal value for vartex
    real(RP) :: eps_area   = 0.0_RP    ! marginal value for triangle area

    integer :: rgnid
    integer :: n, k, l, t, i, j
    !---------------------------------------------------------------------------

    k = ADM_KNONE

    do i=1,imax
       coslon(i) = cos(lon(i))
       sinlon(i) = sin(lon(i))
    enddo
    do j=1,jmax
       coslat(j) = cos(lat(j))
       sinlat(j) = sin(lat(j))
    enddo

    nmax_llgrid        = 0
    nmax_llgrid_rgn(:) = 0

    do l = 1, ADM_lall
    do n = 1, ADM_IooJoo_nmax
    do t = ADM_TI, ADM_TJ

       if ( t == ADM_TI ) then
          r1(:) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),k,l,:) / GRD_rscale
          r2(:) = GRD_x(ADM_IooJoo(n,ADM_GIpJo),k,l,:) / GRD_rscale
          r3(:) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),k,l,:) / GRD_rscale

          lat1 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIoJo),k,l,GMTR_P_LAT)
          lon1 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIoJo),k,l,GMTR_P_LON)
          lat2 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIpJo),k,l,GMTR_P_LAT)
          lon2 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIpJo),k,l,GMTR_P_LON)
          lat3 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIpJp),k,l,GMTR_P_LAT)
          lon3 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIpJp),k,l,GMTR_P_LON)
       else !--- ADM_TJ
          r1(:) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),k,l,:) / GRD_rscale
          r2(:) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),k,l,:) / GRD_rscale
          r3(:) = GRD_x(ADM_IooJoo(n,ADM_GIoJp),k,l,:) / GRD_rscale

          lat1 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIoJo),k,l,GMTR_P_LAT)
          lon1 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIoJo),k,l,GMTR_P_LON)
          lat2 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIpJp),k,l,GMTR_P_LAT)
          lon2 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIpJp),k,l,GMTR_P_LON)
          lat3 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIoJp),k,l,GMTR_P_LAT)
          lon3 = GMTR_P_ll(ADM_IooJoo(n,ADM_GIoJp),k,l,GMTR_P_LON)
       endif

       latmax_l = max(lat1,lat2,lat3) + eps_latlon
       latmin_l = min(lat1,lat2,lat3) - eps_latlon

       if( latmin_l >  polar_limit ) latmax_l =  CNST_PI
       if( latmax_l < -polar_limit ) latmin_l = -CNST_PI

       lonmax_l = max(lon1,lon2,lon3)
       lonmin_l = min(lon1,lon2,lon3)
       if ( lonmax_l-lonmin_l > CNST_PI ) then
          if( lon1 < 0 ) lon1 = lon1 + 2.0_RP * CNST_PI
          if( lon2 < 0 ) lon2 = lon2 + 2.0_RP * CNST_PI
          if( lon3 < 0 ) lon3 = lon3 + 2.0_RP * CNST_PI

          lonmax_l = max(lon1,lon2,lon3)
          lonmin_l = min(lon1,lon2,lon3)
       endif
       lonmax_l = lonmax_l + eps_latlon
       lonmin_l = lonmin_l - eps_latlon

       do j = 1, jmax

          if( lat(j) > latmax_l ) cycle
          if( lat(j) < latmin_l ) cycle

          near_pole = .false.
          if( lat(j) >  polar_limit ) near_pole = .true.
          if( lat(j) < -polar_limit ) near_pole = .true.

          do i = 1, imax

             if ( .NOT. near_pole ) then
                if ( .NOT. (      (       ( lon(i)                  <= lonmax_l ) &
                                    .AND. ( lon(i)                  >= lonmin_l ) ) &
                             .OR. (       ( lon(i) - 2.0_RP*CNST_PI <= lonmax_l ) &
                                    .AND. ( lon(i) - 2.0_RP*CNST_PI >= lonmin_l ) ) &
                             .OR. (       ( lon(i) + 2.0_RP*CNST_PI <= lonmax_l ) &
                                    .AND. ( lon(i) + 2.0_RP*CNST_PI >= lonmin_l ) ) ) ) then
                   cycle
                endif
             endif

             !--- target latlon point on the sphere
             r0(1) = coslat(j) * coslon(i)
             r0(2) = coslat(j) * sinlon(i)
             r0(3) = sinlat(j)

             !--- remove the case inner product is negative
             call VECTR_dot( ip, o(:), r1(:), o(:), r0(:) )
             if( ip < 0.0_RP ) cycle
             v01(:) = r1(:) - r0(:)

             !--- normal vector
             call VECTR_cross( nvec(:), r1(:), r2(:), r2(:), r3(:) )
             call VECTR_abs( len, nvec(:) )

             nvec(:) = nvec(:) / len

             !------ distance from origin to a plane with r0.
             call VECTR_dot( rf, o(:), nvec(:), o(:), r0(:) )
             !------ distance from origin to a plane with r1 or (r2,r3).
             call VECTR_dot( rn, o(:), nvec(:), o(:), r1(:) )

             !------ mapping r0
             r0(1) = r0(1) * (rn/rf)
             r0(2) = r0(2) * (rn/rf)
             r0(3) = r0(3) * (rn/rf)

             !--- calculate vectors from triangler points
             call VECTR_cross( v12xv10(:), r1(:), r2(:), r0(:), r1(:) )
             call VECTR_cross( v23xv20(:), r2(:), r3(:), r0(:), r2(:) )
             call VECTR_cross( v31xv30(:), r3(:), r1(:), r0(:), r3(:) )

             call VECTR_dot( judge12, o(:), nvec(:), o(:), v12xv10(:) )
             call VECTR_dot( judge23, o(:), nvec(:), o(:), v23xv20(:) )
             call VECTR_dot( judge31, o(:), nvec(:), o(:), v31xv30(:) )

             if (       judge12 < eps_judge &
                  .AND. judge23 < eps_judge &
                  .AND. judge31 < eps_judge ) then ! in the triangle

                select case( trim(what_is_done) )
                case( 'GET_NUM' )

                   nmax_llgrid        = nmax_llgrid        + 1
                   nmax_llgrid_rgn(l) = nmax_llgrid_rgn(l) + 1
                   checkmap(i,j) = checkmap(i,j) + 1.0

                case('SET_INDEX')

                   nmax_llgrid        = nmax_llgrid        + 1
                   nmax_llgrid_rgn(l) = nmax_llgrid_rgn(l) + 1

                   lon_index(nmax_llgrid) = i
                   lat_index(nmax_llgrid) = j
                   l_index  (nmax_llgrid) = l
                   t_index  (nmax_llgrid) = t
                   if ( t == ADM_TI ) then
                      n1_index(nmax_llgrid) = ADM_IooJoo(n,ADM_GIoJo)
                      n2_index(nmax_llgrid) = ADM_IooJoo(n,ADM_GIpJo)
                      n3_index(nmax_llgrid) = ADM_IooJoo(n,ADM_GIpJp)
                   else !--- ADM_TJ
                      n1_index(nmax_llgrid) = ADM_IooJoo(n,ADM_GIoJo)
                      n2_index(nmax_llgrid) = ADM_IooJoo(n,ADM_GIpJp)
                      n3_index(nmax_llgrid) = ADM_IooJoo(n,ADM_GIoJp)
                   endif

                   area1 = VECTR_triangle( r0(:), r2(:), r3(:), &
                                           polygon_type, rscale )
                   area2 = VECTR_triangle( r0(:), r3(:), r1(:), &
                                           polygon_type, rscale )
                   area3 = VECTR_triangle( r0(:), r1(:), r2(:), &
                                           polygon_type, rscale )

                   if (      area1 * 0.0_RP /= 0.0_RP &
                        .OR. area2 * 0.0_RP /= 0.0_RP &
                        .OR. area3 * 0.0_RP /= 0.0_RP ) then ! Nan?
                      write(*,          *) 'Nan! (i,j,n,t,l)=', i,j,n,t,l
                      write(*,          *) '(area1,area2,area3)=', area1,area2,area3
                      write(IO_FID_LOG,*) 'Nan! (i,j,n,t,l)=', i,j,n,t,l
                      write(IO_FID_LOG,*) '(area1,area2,area3)=', area1,area2,area3
                      call ADM_proc_stop
                   endif

                   area_total = area1 + area2 + area3

                   w1(nmax_llgrid) = area1 / area_total
                   w2(nmax_llgrid) = area2 / area_total
                   w3(nmax_llgrid) = area3 / area_total
                endselect

                cycle

             elseif(       t == ADM_TI              &
                     .AND. abs(v01(1)) < eps_vertex &
                     .AND. abs(v01(2)) < eps_vertex &
                     .AND. abs(v01(3)) < eps_vertex ) then ! on the triangle vertex

                select case( trim(what_is_done) )
                case( 'GET_NUM' )

                   nmax_llgrid        = nmax_llgrid        + 1
                   nmax_llgrid_rgn(l) = nmax_llgrid_rgn(l) + 1
                   checkmap(i,j) = checkmap(i,j) + 1.0

                case('SET_INDEX')

                   nmax_llgrid        = nmax_llgrid        + 1
                   nmax_llgrid_rgn(l) = nmax_llgrid_rgn(l) + 1

                   lon_index(nmax_llgrid) = i
                   lat_index(nmax_llgrid) = j
                   l_index  (nmax_llgrid) = l
                   t_index  (nmax_llgrid) = t
                   n1_index (nmax_llgrid) = ADM_IooJoo(n,ADM_GIoJo)
                   n2_index (nmax_llgrid) = ADM_IooJoo(n,ADM_GIpJo)
                   n3_index (nmax_llgrid) = ADM_IooJoo(n,ADM_GIpJp)
                   w1       (nmax_llgrid) = 1.0_RP
                   w2       (nmax_llgrid) = 0.0_RP
                   w3       (nmax_llgrid) = 0.0_RP
                endselect

             endif

          enddo ! i LOOP
       enddo ! j LOOP


    enddo ! TI,TJ
    enddo ! n LOOP
    enddo ! l LOOP

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       if ( rgnid == ADM_rgnid_npl_mng ) then
          n = ADM_gall_1d * ADM_gmax + ADM_gmin
       elseif ( rgnid == ADM_rgnid_spl_mng ) then
          n = ADM_gall_1d * (ADM_gmin-1) + ADM_gmax+1
       else
          cycle
       endif

       r1(:) = GRD_x(n,k,l,:) / GRD_rscale

       do j = 1, jmax
       do i = 1, imax
          !--- target latlon point on the sphere
          r0(1) = coslat(j) * coslon(i)
          r0(2) = coslat(j) * sinlon(i)
          r0(3) = sinlat(j)

          !--- remove the case inner product is negative
          call VECTR_dot( ip, o(:), r1(:), o(:), r0(:) )
          if( ip < 0.0_RP ) cycle
          v01(:) = r1(:) - r0(:)

          if (      abs(v01(1)) < eps_vertex &
              .AND. abs(v01(2)) < eps_vertex &
              .AND. abs(v01(3)) < eps_vertex ) then ! on the pole

             select case( trim(what_is_done) )
             case( 'GET_NUM' )

                nmax_llgrid        = nmax_llgrid        + 1
                nmax_llgrid_rgn(l) = nmax_llgrid_rgn(l) + 1
                checkmap(i,j) = checkmap(i,j) + 1.0

             case('SET_INDEX')

                nmax_llgrid        = nmax_llgrid        + 1
                nmax_llgrid_rgn(l) = nmax_llgrid_rgn(l) + 1

                lon_index(nmax_llgrid) = i
                lat_index(nmax_llgrid) = j
                l_index  (nmax_llgrid) = l
                t_index  (nmax_llgrid) = 0
                n1_index (nmax_llgrid) = n
                n2_index (nmax_llgrid) = n
                n3_index (nmax_llgrid) = n
                w1       (nmax_llgrid) = 1.0_RP
                w2       (nmax_llgrid) = 0.0_RP
                w3       (nmax_llgrid) = 0.0_RP
             endselect

          endif

       enddo ! i LOOP
       enddo ! j LOOP
    enddo ! l LOOP


    return
  end subroutine mkrelmap_ico2ll

  !-----------------------------------------------------------------------------
  !> Output sample output
  subroutine LL_outputsample
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_prc_tab,   &
       ADM_prc_me,    &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_gall_1d,   &
       ADM_gmax,      &
       ADM_gmin,      &
       ADM_KNONE
    use mod_comm, only: &
       COMM_data_transfer
    use mod_io_param, only: &
       IO_REAL8
    use mod_fio, only: &
       FIO_output
    implicit none

    real(RP) :: SAMPLE   ( ADM_gall,   ADM_KNONE,ADM_lall,   4)
    real(RP) :: SAMPLE_pl( ADM_gall_pl,ADM_KNONE,ADM_lall_pl,4)

    character(len=H_LONG) :: fname

    integer :: fid
    integer :: rgnid, prc
    integer :: i, j, ij, k, l
    !---------------------------------------------------------------------------

    k = ADM_KNONE

    SAMPLE   (:,:,:,:) = -999.0_RP
    SAMPLE_pl(:,:,:,:) = -999.0_RP

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)
       prc   = ADM_prc_me

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij = ADM_gall_1d * (j-1) + i

          SAMPLE(ij,k,l,1) = real(prc,  kind=RP)
          SAMPLE(ij,k,l,2) = real(rgnid,kind=RP)
          SAMPLE(ij,k,l,3) = real(i,    kind=RP)
          SAMPLE(ij,k,l,4) = real(j,    kind=RP)
       enddo
       enddo
    enddo

    do l = 1, ADM_lall_pl
       rgnid = l
       prc   = ADM_prc_me

       do ij = 1, ADM_gall_pl
          SAMPLE_pl(ij,k,l,1) = real(prc,   kind=RP)
          SAMPLE_pl(ij,k,l,2) = real(-rgnid,kind=RP)
          SAMPLE_pl(ij,k,l,3) = real( ij,   kind=RP)
          SAMPLE_pl(ij,k,l,4) = real(-ij,   kind=RP)
       enddo
    enddo

    call COMM_data_transfer( SAMPLE, SAMPLE_pl )

    if ( SAMPLE_io_mode == 'ADVANCED' ) then

       call FIO_output( SAMPLE(:,:,:,1), SAMPLE_OUT_BASENAME, "", "", &
                       "sample1", "sample data(prc)", "", "NIL",      &
                       IO_REAL8, "ZSSFC1", k, k, 1, 0.0_DP, 0.0_DP    )
       call FIO_output( SAMPLE(:,:,:,2), SAMPLE_OUT_BASENAME, "", "", &
                       "sample2", "sample data(rgn)", "", "NIL",      &
                       IO_REAL8, "ZSSFC1", k, k, 1, 0.0_DP, 0.0_DP    )
       call FIO_output( SAMPLE(:,:,:,3), SAMPLE_OUT_BASENAME, "", "", &
                       "sample3", "sample data(i)", "", "NIL",        &
                       IO_REAL8, "ZSSFC1", k, k, 1, 0.0_DP, 0.0_DP    )
       call FIO_output( SAMPLE(:,:,:,4), SAMPLE_OUT_BASENAME, "", "", &
                       "sample4", "sample data(j)", "", "NIL",        &
                       IO_REAL8, "ZSSFC1", k, k, 1, 0.0_DP, 0.0_DP    )

    elseif( sample_io_mode == 'LEGACY' ) then

       do l = 1, ADM_lall
          rgnid = ADM_prc_tab(l,ADM_prc_me)
          call IO_make_idstr(fname,trim(sample_io_mode),'rgn',rgnid)

          fid = IO_get_available_fid()
          open( unit = fid, &
               file=trim(fname),   &
               form='unformatted', &
               access='direct',    &
               recl=ADM_gall*8     )

             write(fid,rec=1) SAMPLE(:,k,l,1)
             write(fid,rec=2) SAMPLE(:,k,l,2)
             write(fid,rec=3) SAMPLE(:,k,l,3)
             write(fid,rec=4) SAMPLE(:,k,l,4)
          close(fid)
       enddo
    else
       write(IO_FID_LOG,*) 'Invalid io_mode!'
       call ADM_proc_stop
    endif

    return
  end subroutine LL_outputsample

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine setup_latlon
  !>
  subroutine setup_latlon
    use mod_cnst, only: &
         CNST_PI
    implicit none
    !---------------------------------------------------------------------------

    if ( latlon_type == 'EQUIDIST' ) then

       call set_equidist_grid

    elseif( latlon_type == 'GAUSSIAN' ) then

       !--- if GAUSSIAN grid is selected, latmax and latmin are overwritten.
       latmax =  CNST_PI*0.5_RP
       latmin = -CNST_PI*0.5_RP
       call set_gaussian_grid

    endif

    return
  end subroutine setup_latlon

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine set_equidist_grid
  !>
  subroutine set_equidist_grid
    implicit none

    real(RP) :: dlat, dlon
    integer :: i, j
    !---------------------------------------------------------------------------

    dlat = ( latmax - latmin ) / real(jmax,kind=RP)

    do j = 1, jmax
       lat(j) = latmin + dlat * ( real(j,kind=RP) - 0.5_RP )
    enddo

    dlon = ( lonmax - lonmin ) / real(imax,kind=RP)

    if ( lon_offset ) then
       do i = 1, imax
          lon(i) = lonmin + dlon * ( real(i,kind=RP) - 0.5_RP )
       enddo
    else
       do i = 1, imax
          lon(i) = lonmin + dlon * ( real(i,kind=RP) - 1.0_RP )
       enddo
    endif

    return
  end subroutine set_equidist_grid

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine set_gaussian_grid
  !>
  subroutine set_gaussian_grid
    use mod_cnst, only: &
         CNST_PI
    implicit none

    integer, parameter :: nb = 256

    real(RP) :: e(nb)
    real(RP) :: eps

    real(RP) :: mu0, dmu, dP0
    real(RP) :: mu(jmax)
    real(RP) :: P0(0:jmax)

    real(RP) :: dlon

    integer :: n, i, j
    !---------------------------------------------------------------------------

    !--- calculate machine eps.
    e(:) = 0.0_RP
    eps  = 1.0_RP

    do i = 1, nb
       eps  = eps * 0.5_RP
       e(i) = eps+1
    enddo

    i = 1
    eps = 1.0_RP

    do i = 1, nb
       if( e(i) > 1.0_RP ) exit
       eps = eps * 0.5_RP
    enddo
    eps = eps * 4.0_RP

    !---  calculate gausian_grid by Newton-Rapson method
    do j = 1, jmax

       mu0=sin(CNST_PI*real(jmax+1-2*j,kind=RP)/real(2*jmax+1,kind=RP))

       loopeps:do

          P0(0)=1.0_RP
          P0(1)=mu0
          do n=1,jmax-1
             P0(n+1)=(real(2*n+1,kind=RP)*mu0*P0(n)-n*P0(n-1)) / real(n+1,kind=RP)
          enddo
          dP0=jmax*(P0(jmax-1)-mu0*P0(jmax))/(1-mu0*mu0)
          dmu=P0(jmax)/dP0
          mu0=mu0-dmu
          if (abs(dmu)<eps) then
             mu(j)=mu0
             P0(0)=1.0_RP
             P0(1)=mu(j)
             do n=1,jmax-1
                P0(n+1)=(real(2*n+1,kind=RP)*mu(j)*P0(n)-n*P0(n-1)) / real(n+1,kind=RP)
             enddo
             exit loopeps
          endif

       enddo loopeps

    enddo

    do j = 1, jmax
       lat(j) = -asin(mu(j))
    enddo

    dlon = ( lonmax - lonmin) / real(imax,kind=RP)

    if ( lon_offset ) then
       do i = 1, imax
          lon(i) = lonmin + dlon * ( real(i,kind=RP) - 0.5_RP )
       enddo
    else
       do i = 1, imax
          lon(i) = lonmin + dlon * ( real(i,kind=RP) - 1.0_RP )
       enddo
    endif

    return
  end subroutine set_gaussian_grid

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine intrpl_2
  !>
!  subroutine intrpl_2( var_ll, var, var_pl, kmin, kmax )
!    !
!    use mod_adm, only :              &
!         ADM_prc_me,               &
!         ADM_prc_pl,               &
!         ADM_GSLF_PL,             &
!         ADM_gall_pl,              &
!         ADM_lall_pl,              &
!         ADM_IooJoo_nmax,          &
!         ADM_IooJoo,               &
!         ADM_GIoJo,                &
!         ADM_GIpJo,                &
!         ADM_GIpJp,                &
!         ADM_GIoJp,                &
!         ADM_GImJo,                &
!         ADM_GIoJm,                &
!         ADM_GImJm,                &
!         ADM_KNONE,                &
!         ADM_kall,                 &
!         ADM_gall,                 &
!         ADM_lall
!    use mod_grd, only :               &
!         GRD_x,GRD_x_pl,           &
!         GRD_XDIR,                 &
!         GRD_YDIR,                 &
!         GRD_ZDIR
!    use mod_cnst, only :               &
!         CNST_UNDEF
!    use mod_oprt, only :           &
!         OPRT_gradient
!    use mod_comm, only :               &
!         COMM_data_transfer
!    !
!    implicit none
!    !
!    integer, intent(in) :: kmin,kmax
!    real(RP), intent(out) :: var_ll(nmax_llgrid,kmin:kmax)
!    !    real(RP), intent(in)  :: var(ADM_gall,kmin:kmax,ADM_lall)
!    !    real(RP), intent(in)  :: var_pl(ADM_gall_pl,kmin:kmax,ADM_lall_pl)
!    !    real(RP), intent(in)  :: var(:,:,:)
!    !    real(RP), intent(in)  :: var_pl(:,:,:)
!    real(RP), intent(in) :: var(ADM_gall,ADM_kall,ADM_lall)
!    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!    !
!    integer, parameter :: ix = 1
!    integer, parameter :: iy = 2
!    integer, parameter :: iz = 3
!    !
!    real(RP) ::  grd(ADM_gall,ADM_kall,ADM_lall,ix:iz)
!    real(RP) ::  grd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,ix:iz)
!    logical ::  def_grd(ADM_gall,kmin:kmax,ADM_lall)
!    !
!    integer :: n0,n1,n2,n3,n4,n5,n6
!    real(RP) :: vec12(ix:iz),vec23(ix:iz),vec31(ix:iz)
!    real(RP) :: v1,v2,v3,v12,v23,v31
!    !
!    integer :: k,m,l,n
!    !
!    grd(:,:,:,:)=CNST_UNDEF
!    grd_pl(:,:,:,:)=CNST_UNDEF
!    call OPRT_gradient(&
!         grd(:,:,:,ix), grd_pl(:,:,:,ix),&
!         grd(:,:,:,iy), grd_pl(:,:,:,iy),&
!         grd(:,:,:,iz), grd_pl(:,:,:,iz),&
!         var,  var_pl)
!    !
!    def_grd(:,:,:) = .false.
!    !
!    ! cheange by kgoto
!    ! check gradient operator by checking variable value
!    do l=1,ADM_lall
!       do k=kmin,kmax
!          do n=1,ADM_IooJoo_nmax
!             n0=ADM_IooJoo(n,ADM_GIoJo)
!             n1=ADM_IooJoo(n,ADM_GIpJo)
!             n2=ADM_IooJoo(n,ADM_GIpJp)
!             n3=ADM_IooJoo(n,ADM_GIoJp)
!             n4=ADM_IooJoo(n,ADM_GImJo)
!             n5=ADM_IooJoo(n,ADM_GImJm)
!             n6=ADM_IooJoo(n,ADM_GIoJm)
!             if ( (var(n0,k,l)/=CNST_UNDEF) .and. &
!                  (var(n1,k,l)/=CNST_UNDEF) .and. &
!                  (var(n2,k,l)/=CNST_UNDEF) .and. &
!                  (var(n3,k,l)/=CNST_UNDEF) .and. &
!                  (var(n4,k,l)/=CNST_UNDEF) .and. &
!                  (var(n5,k,l)/=CNST_UNDEF) .and. &
!                  (var(n6,k,l)/=CNST_UNDEF) ) then
!             else
!                grd(n0,k,l,ix:iz)=CNST_UNDEF
!             endif
!          enddo
!       enddo
!    enddo
!    if (adm_prc_me==adm_prc_pl) then
!       do l=1,ADM_lall_pl
!          do k=kmin,kmax
!             if ( (var_pl(1,k,l)/=CNST_UNDEF) .and. &
!                  (var_pl(2,k,l)/=CNST_UNDEF) .and. &
!                  (var_pl(3,k,l)/=CNST_UNDEF) .and. &
!                  (var_pl(4,k,l)/=CNST_UNDEF) .and. &
!                  (var_pl(5,k,l)/=CNST_UNDEF) .and. &
!                  (var_pl(6,k,l)/=CNST_UNDEF) ) then
!             else
!                grd_pl(adm_gslf_pl,k,l,ix:iz)=CNST_UNDEF
!             endif
!          enddo
!       enddo
!    endif
!    !
!    call COMM_data_transfer(grd,grd_pl)
!    !
!    do k=kmin,kmax
!       do m=1,nmax_llgrid
!          n1 = n1_index(m)
!          n2 = n2_index(m)
!          n3 = n3_index(m)
!          l = l_index(m)
!          if ( (grd(n1,k,l,ix)/=CNST_UNDEF) .and. &
!               (grd(n1,k,l,iy)/=CNST_UNDEF) .and. &
!               (grd(n1,k,l,iz)/=CNST_UNDEF) ) then
!             def_grd(n1,k,l)=.true.
!          endif
!          if ( (grd(n2,k,l,ix)/=CNST_UNDEF) .and. &
!               (grd(n2,k,l,iy)/=CNST_UNDEF) .and. &
!               (grd(n2,k,l,iz)/=CNST_UNDEF) ) then
!             def_grd(n2,k,l)=.true.
!          endif
!          if ( (grd(n3,k,l,ix)/=CNST_UNDEF) .and. &
!               (grd(n3,k,l,iy)/=CNST_UNDEF) .and. &
!               (grd(n3,k,l,iz)/=CNST_UNDEF) ) then
!             def_grd(n3,k,l)=.true.
!          endif
!       enddo
!    enddo
!
!    do k=kmin, kmax
!       do m=1,nmax_llgrid
!          n1 = n1_index(m)
!          n2 = n2_index(m)
!          n3 = n3_index(m)
!          l = l_index(m)
!          !
!          !--- vector n1 to n2
!          vec12(ix) &
!               = GRD_x(n2,ADM_KNONE,l,GRD_XDIR)&
!               - GRD_x(n1,ADM_KNONE,l,GRD_XDIR)
!          vec12(iy) &
!               = GRD_x(n2,ADM_KNONE,l,GRD_YDIR)&
!               - GRD_x(n1,ADM_KNONE,l,GRD_YDIR)
!          vec12(iz) &
!               = GRD_x(n2,ADM_KNONE,l,GRD_ZDIR)&
!               - GRD_x(n1,ADM_KNONE,l,GRD_ZDIR)
!          !
!          !--- vector n2 to n3
!          vec23(ix) &
!               = GRD_x(n3,ADM_KNONE,l,GRD_XDIR)&
!               - GRD_x(n2,ADM_KNONE,l,GRD_XDIR)
!          vec23(iy) &
!               = GRD_x(n3,ADM_KNONE,l,GRD_YDIR)&
!               - GRD_x(n2,ADM_KNONE,l,GRD_YDIR)
!          vec23(iz) &
!               = GRD_x(n3,ADM_KNONE,l,GRD_ZDIR)&
!               - GRD_x(n2,ADM_KNONE,l,GRD_ZDIR)
!          !
!          !--- vector n3 to n1
!          vec31(ix) &
!               = GRD_x(n1,ADM_KNONE,l,GRD_XDIR)&
!               - GRD_x(n3,ADM_KNONE,l,GRD_XDIR)
!          vec31(iy) &
!               = GRD_x(n1,ADM_KNONE,l,GRD_YDIR)&
!               - GRD_x(n3,ADM_KNONE,l,GRD_YDIR)
!          vec31(iz) &
!               = GRD_x(n1,ADM_KNONE,l,GRD_ZDIR)&
!               - GRD_x(n3,ADM_KNONE,l,GRD_ZDIR)
!          !
!          !--- calculate values at the nodes
!          v1=var(n1,k,l)
!          v2=var(n2,k,l)
!          v3=var(n3,k,l)
!          !
!          !--- midpoint value between n1 and n2
!          !--- based on Hermitian interpolation.
!          v12 = (v1+v2)*0.5_RP                                   &
!               + 0.125_RP                                        &
!               * (( grd(n1,k,l,ix) - grd(n2,k,l,ix) )*vec12(ix) &
!                 +( grd(n1,k,l,iy) - grd(n2,k,l,iy) )*vec12(iy) &
!                 +( grd(n1,k,l,iz) - grd(n2,k,l,iz) )*vec12(iz))
!          !
!          !--- midpoint value between n2 and n3
!          !--- based on Hermitian interpolation.
!          v23 = (v2+v3)*0.5_RP                                   &
!               + 0.125_RP                                        &
!               * (( grd(n2,k,l,ix) - grd(n3,k,l,ix) )*vec23(ix) &
!                 +( grd(n2,k,l,iy) - grd(n3,k,l,iy) )*vec23(iy) &
!                 +( grd(n2,k,l,iz) - grd(n3,k,l,iz) )*vec23(iz))
!          !
!          !--- midpoint value between n3 and n1
!          !--- based on Hermitian interpolation.
!          v31 = (v3+v1)*0.5_RP                                   &
!               + 0.125_RP                                        &
!               * (( grd(n3,k,l,ix) - grd(n1,k,l,ix) )*vec31(ix) &
!                 +( grd(n3,k,l,iy) - grd(n1,k,l,iy) )*vec31(iy) &
!                 +( grd(n3,k,l,iz) - grd(n1,k,l,iz) )*vec31(iz))
!          if(  def_grd(n1,k,l).and. &
!               def_grd(n2,k,l).and. &
!               def_grd(n3,k,l) ) then
!             var_ll(m,k)                          &
!                  = w1(m)*(2.0_RP*w1(m)-1.0_RP)*v1  &
!                  + w2(m)*(2.0_RP*w2(m)-1.0_RP)*v2  &
!                  + w3(m)*(2.0_RP*w3(m)-1.0_RP)*v3  &
!                  + 4.0_RP*w1(m)*w2(m)        *v12 &
!                  + 4.0_RP*w2(m)*w3(m)        *v23 &
!                  + 4.0_RP*w3(m)*w1(m)        *v31
!          else
!             var_ll(m,k) = CNST_UNDEF
!          endif
!       enddo
!    enddo
!    !
!    return
!    !
!  end subroutine intrpl_2

end module mod_latlon
!-------------------------------------------------------------------------------
