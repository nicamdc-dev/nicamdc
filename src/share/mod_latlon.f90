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
  use mod_adm, only: &
     ADM_NSYS,    &
     ADM_MAXFNAME
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LATLON_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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
  character(len=ADM_NSYS), private, save :: latlon_type = 'EQUIDIST' ! grid type ( equidist or gaussian )
  integer,                 private, save :: imax       = 360         ! number of longitude
  integer,                 private, save :: jmax       = 180         ! number of latitude
  real(8),                 private, save :: lonmin     = -999.D0     ! minimun longitude of region window in deg
  real(8),                 private, save :: lonmax     = -999.D0     ! maximun longitude of region window in deg
  real(8),                 private, save :: latmin     = -999.D0     ! minimun latitude of region window in deg
  real(8),                 private, save :: latmax     = -999.D0     ! maximun latitude of region window in deg
  logical,                 private, save :: lon_offset = .true.      ! logitude offset

  character(len=ADM_NSYS),     private, save :: output_lldata_type = 'mkllmap'

  real(8), private, allocatable, save :: lat(:)
  real(8), private, allocatable, save :: lon(:)

  integer,              private, save :: nmax_llgrid
  integer, private, allocatable, save :: nmax_llgrid_rgn(:)

  integer, private, allocatable, save :: lon_index(:)
  integer, private, allocatable, save :: lat_index(:)
  integer, private, allocatable, save :: l_index  (:)
  integer, private, allocatable, save :: t_index  (:)
  integer, private, allocatable, save :: n1_index (:)
  integer, private, allocatable, save :: n2_index (:)
  integer, private, allocatable, save :: n3_index (:)
  real(8), private, allocatable, save :: w1       (:)
  real(8), private, allocatable, save :: w2       (:)
  real(8), private, allocatable, save :: w3       (:)

  logical, private, save :: debug = .false.
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine LATLON_setup
  !>
  subroutine LATLON_setup( output_dirname, output_lldata_type_in )
    use mod_misc, only: &
       MISC_get_available_fid, &
       MISC_make_idstr
    use mod_adm, only: &
       ADM_prc_all,    &
       ADM_COMM_RUN_WORLD, &
       ADM_LOG_FID,        &
       ADM_CTL_FID,        &
       ADM_proc_stop,      &
       ADM_prc_tab,        &
       ADM_prc_run_master, &
       ADM_prc_me,         &
       ADM_lall
    use mod_cnst, only: &
       CNST_PI
    implicit none

    character(len=*), intent(in) :: output_dirname
    character(len=*), intent(in) :: output_lldata_type_in

    real(8) :: latmax_deg =   90.D0
    real(8) :: latmin_deg =  -90.D0
    real(8) :: lonmax_deg =  180.D0
    real(8) :: lonmin_deg = -180.D0

    namelist / LATLONPARAM / &
         latlon_type, &
         imax,        &
         jmax,        &
         lonmin_deg,  &
         lonmax_deg,  &
         latmin_deg,  &
         latmax_deg,  &
         lon_offset,  &
         debug

    character(len=ADM_MAXFNAME) :: fname

    real(8) :: d2r

    integer :: globalsum
    integer :: gathersum(ADM_prc_all)

    integer :: fid, ierr
    integer :: nstart, nend
    integer :: n, l, rgnid
    !---------------------------------------------------------------------------

    output_lldata_type = output_lldata_type_in

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[latlon]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=LATLONPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** LATLONPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist LATLONPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist LATLONPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,LATLONPARAM)

    d2r    = CNST_PI / 180.D0
    latmax = latmax_deg * d2r
    latmin = latmin_deg * d2r
    lonmax = lonmax_deg * d2r
    lonmin = lonmin_deg * d2r

    !--- setup latitude-longitude grid
    allocate( lat(jmax) )
    allocate( lon(imax) )

    call setup_latlon

    if( ADM_prc_me == ADM_prc_run_master ) then
       fid = MISC_get_available_fid()
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

    write(ADM_LOG_FID,*) '====== Lat-Lon grid info. ======'
    if ( latlon_type == 'EQUIDIST' ) then
       write(ADM_LOG_FID,*) '--- Latitude  type   : Equal distance'
    elseif( latlon_type == 'GAUSSIAN' ) then
       write(ADM_LOG_FID,*) '--- Latitude  type   : Gaussian'
    endif
    if ( lon_offset ) then
       write(ADM_LOG_FID,*) '--- Longitude offset : yes'
    else
       write(ADM_LOG_FID,*) '--- Longitude offset : no'
    endif
    write(ADM_LOG_FID,*)    '--- # of Latitude    :', jmax
    write(ADM_LOG_FID,*)    '--- # of Longitude   :', imax
    write(ADM_LOG_FID,*)    '--- Latitude  range  :', lat(1)/d2r,' - ', lat(jmax)/d2r
    write(ADM_LOG_FID,*)    '--- Longitude range  :', lon(1)/d2r,' - ', lon(imax)/d2r

    allocate( nmax_llgrid_rgn(ADM_lall) )

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '*** Start counting grid.'

    ! count lat-lon number
    call mkrelmap_ico2ll( 'GET_NUM' )

    write(ADM_LOG_FID,*) '# of managing llgrid'
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       write(ADM_LOG_FID,*) 'region=', rgnid, ', llgrid=', nmax_llgrid_rgn(l)
    enddo

    ! check total lat-lon number
    gathersum(ADM_prc_me) = nmax_llgrid

    call MPI_Barrier(ADM_COMM_RUN_WORLD,ierr)

    call MPI_Allgather( MPI_IN_PLACE,       &
                        1,                  &
                        MPI_INTEGER,        &
                        gathersum(1),       &
                        1,                  &
                        MPI_INTEGER,        &
                        ADM_COMM_RUN_WORLD, &
                        ierr                )

    globalsum = sum( gathersum(:) )

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) 'imax x jmax                    = ', imax*jmax
    write(ADM_LOG_FID,*) 'global total of counted llgird = ', globalsum
    if ( globalsum /= imax*jmax ) then
       write(*,          *) 'counted llgrid does not match!'
       write(ADM_LOG_FID,*) 'counted llgrid does not match!'
       call ADM_proc_stop
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

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '*** Start calc relation map.'

    ! calc relation map
    call mkrelmap_ico2ll( 'SET_INDEX' )

    write(ADM_LOG_FID,*) '*** OK.'

    ! output relation map
    if ( debug ) then
       do l = 1, ADM_lall
          nend   = sum(nmax_llgrid_rgn(1:l))
          nstart = nend - nmax_llgrid_rgn(l) + 1

          do n = nstart, nend
             if ( abs(w1(n)+w2(n)+w3(n)-1.D0) > 1.D-15 ) then
                write(ADM_LOG_FID,'(A,2I6,E30.20)') '(lat,lon,area)=', &
                lat_index(n),lon_index(n),w1(n)+w2(n)+w3(n)
             endif
          enddo
       enddo
    endif

    ! output relation map
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       call MISC_make_idstr(fname,trim(output_dirname)//'/llmap','rgn',rgnid)

       fid = MISC_get_available_fid()
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
    ! output relation map
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       call MISC_make_idstr(fname,trim(output_dirname)//'/llmap','rgntxt',rgnid)

       fid = MISC_get_available_fid()
       open( unit   = fid,           &
             file   = trim(fname),   &
             form   = 'formatted', &
             status = 'unknown'      )

          nend   = sum(nmax_llgrid_rgn(1:l))
          nstart = nend - nmax_llgrid_rgn(l) + 1

          write(fid,*) nmax_llgrid_rgn(l)
          if ( nmax_llgrid_rgn(l) /= 0 ) then
             write(fid,*) lon_index(nstart:nend)
             write(fid,*) lat_index(nstart:nend)
             write(fid,*) n1_index (nstart:nend)
             write(fid,*) n2_index (nstart:nend)
             write(fid,*) n3_index (nstart:nend)
             write(fid,*) w1       (nstart:nend)
             write(fid,*) w2       (nstart:nend)
             write(fid,*) w3       (nstart:nend)
          endif
       close(fid)
    enddo
 
    return
  end subroutine LATLON_setup

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine mkrelmap_ico2ll
  !>
  subroutine mkrelmap_ico2ll( what_is_done )
    use mod_misc, only: &
       MISC_triangle_area_q, &
       MISC_3dvec_triangle,  &
       MISC_3dvec_cross,     &
       MISC_3dvec_dot,       &
       MISC_3dvec_abs
    use mod_adm, only: &
       ADM_proc_stop,   &
       ADM_LOG_FID,     &
       ADM_TI,          &
       ADM_TJ,          &
       ADM_lall,        &
       ADM_KNONE,       &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GIpJo,       &
       ADM_GIpJp,       &
       ADM_GIoJp
    use mod_cnst, only: &
       CNST_PI
    use mod_grd, only: &
       GRD_rscale,       &
       GRD_x
    use mod_gmtr, only: &
       P_LAT => GMTR_P_LAT, &
       P_LON => GMTR_P_LON, &
       GMTR_polygon_type,   &
       GMTR_P_var
    implicit none

    character(len=*), intent(in) :: what_is_done

    real(8), parameter :: rscale = 1.D0

    real(8), parameter :: o(3) = 0.D0
    real(8) :: r0(3), r1(3), r2(3), r3(3)
    real(8) :: nvec(3), v12xv10(3), v23xv20(3), v31xv30(3)
    real(8) :: judge12, judge23, judge31
    real(8) :: rf, rn, ip, len
    real(8) :: v01(3)

    real(8) :: coslat(jmax), sinlat(jmax)
    real(8) :: coslon(imax), sinlon(imax)
    real(8) :: lat1, lat2, lat3, lat4
    real(8) :: lon1,lon2,lon3,lon4
    real(8) :: latmin_l,latmax_l
    real(8) :: lonmin_l,lonmax_l

    real(8) :: area_total, area1, area2, area3

    real(8) :: eps_judge  = 1.D-10 ! marginal value for inner products
    real(8) :: eps_latlon = 1.D-10 ! marginal square near grid points (in radian)
    real(8) :: epsi_area  = 0.D0   ! marginal value for triangle area

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

       lat1 = GMTR_P_var(ADM_IooJoo(n,ADM_GIoJo),k,l,P_LAT)
       lat2 = GMTR_P_var(ADM_IooJoo(n,ADM_GIpJo),k,l,P_LAT)
       lat3 = GMTR_P_var(ADM_IooJoo(n,ADM_GIpJp),k,l,P_LAT)
       lat4 = GMTR_P_var(ADM_IooJoo(n,ADM_GIoJp),k,l,P_LAT)

       lon1 = GMTR_P_var(ADM_IooJoo(n,ADM_GIoJo),k,l,P_LON)
       lon2 = GMTR_P_var(ADM_IooJoo(n,ADM_GIpJo),k,l,P_LON)
       lon3 = GMTR_P_var(ADM_IooJoo(n,ADM_GIpJp),k,l,P_LON)
       lon4 = GMTR_P_var(ADM_IooJoo(n,ADM_GIoJp),k,l,P_LON)

       latmin_l = min(lat1,lat2,lat3,lat4)
       latmax_l = max(lat1,lat2,lat3,lat4)

       lonmin_l = min(lon1,lon2,lon3,lon4) - eps_latlon
       lonmax_l = max(lon1,lon2,lon3,lon4) + eps_latlon
       if ( lonmax_l-lonmin_l > CNST_PI ) then
          if(lon1 < 0 ) lon1 = lon1 + 2.D0 * CNST_PI
          if(lon2 < 0 ) lon2 = lon2 + 2.D0 * CNST_PI
          if(lon3 < 0 ) lon3 = lon3 + 2.D0 * CNST_PI
          if(lon4 < 0 ) lon4 = lon4 + 2.D0 * CNST_PI

          lonmin_l = min(lon1,lon2,lon3,lon4)
          lonmax_l = max(lon1,lon2,lon3,lon4)
       endif
       lonmin_l = lonmin_l - eps_latlon
       lonmax_l = lonmax_l + eps_latlon

       do t = ADM_TI, ADM_TJ

          if ( t == ADM_TI ) then
             r1(:) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),k,l,:) / GRD_rscale
             r2(:) = GRD_x(ADM_IooJoo(n,ADM_GIpJo),k,l,:) / GRD_rscale
             r3(:) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),k,l,:) / GRD_rscale
          else !--- ADM_TJ
             r1(:) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),k,l,:) / GRD_rscale
             r2(:) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),k,l,:) / GRD_rscale
             r3(:) = GRD_x(ADM_IooJoo(n,ADM_GIoJp),k,l,:) / GRD_rscale
          endif

          do i = 1, imax

             if ( .NOT. (      (       ( lonmin_l <= lon(i)                ) &
                                 .AND. ( lonmax_l >= lon(i)                ) ) &
                          .OR. (       ( lonmin_l <= lon(i) - 2.D0*CNST_PI ) &
                                 .AND. ( lonmax_l >= lon(i) - 2.D0*CNST_PI ) ) &
                          .OR. (       ( lonmin_l <= lon(i) + 2.D0*CNST_PI ) &
                                 .AND. ( lonmax_l >= lon(i) + 2.D0*CNST_PI ) ) ) ) then
                cycle
             endif

          do j = 1, jmax

             if( lat(j) < latmin_l ) cycle
             if( lat(j) > latmax_l ) cycle

             !--- target latlon point on the sphere
             r0(1) = coslat(j) * coslon(i)
             r0(2) = coslat(j) * sinlon(i)
             r0(3) = sinlat(j)

             !--- remove the case inner product is negative
             call MISC_3dvec_dot( ip, o(:), r1(:), o(:), r0(:) )
             if( ip < 0.D0 ) cycle
             v01(:) = r1(:) - r0(:)

             !--- normal vector
             call MISC_3dvec_cross( nvec(:), r1(:), r2(:), r2(:), r3(:) )
             call MISC_3dvec_abs( len, nvec(:) )

             nvec(:) = nvec(:) / len

             !------ distance from origin to a plane with r0.
             call MISC_3dvec_dot( rf, o(:), nvec(:), o(:), r0(:) )
             !------ distance from origin to a plane with r1 or (r2,r3).
             call MISC_3dvec_dot( rn, o(:), nvec(:), o(:), r1(:) )

             !------ mapping r0 
             r0(1) = r0(1) * (rn/rf)
             r0(2) = r0(2) * (rn/rf)
             r0(3) = r0(3) * (rn/rf)

             !--- calculate vectors from triangler points 
             call MISC_3dvec_cross( v12xv10(:), r1(:), r2(:), r0(:), r1(:) )
             call MISC_3dvec_cross( v23xv20(:), r2(:), r3(:), r0(:), r2(:) )
             call MISC_3dvec_cross( v31xv30(:), r3(:), r1(:), r0(:), r3(:) )

             call MISC_3dvec_dot( judge12, o(:), nvec(:), o(:), v12xv10(:) )
             call MISC_3dvec_dot( judge23, o(:), nvec(:), o(:), v23xv20(:) )
             call MISC_3dvec_dot( judge31, o(:), nvec(:), o(:), v31xv30(:) )

             if (       judge12 <= -eps_judge &
                  .AND. judge23 <= -eps_judge &
                  .AND. judge31 <= -eps_judge ) then ! in the triangle

                select case( trim(what_is_done) )
                case( 'GET_NUM' )

                   nmax_llgrid        = nmax_llgrid        + 1
                   nmax_llgrid_rgn(l) = nmax_llgrid_rgn(l) + 1

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
                      n2_index(nmax_llgrid) = ADM_IooJoo(n,ADM_GIpJo)
                      n3_index(nmax_llgrid) = ADM_IooJoo(n,ADM_GIpJp)
                   endif

                   if ( output_lldata_type == 'mkllmap_q' ) then ! quad precision
                      area1 = MISC_triangle_area_q( r0(:), r2(:), r3(:),       &
                                                    GMTR_polygon_type, rscale, &
                                                    critical=epsi_area         )
                      area2 = MISC_triangle_area_q( r0(:), r3(:), r1(:),       &
                                                    GMTR_polygon_type, rscale, &
                                                    critical=epsi_area         )
                      area3 = MISC_triangle_area_q( r0(:), r1(:), r2(:),       &
                                                    GMTR_polygon_type, rscale, &
                                                    critical=epsi_area         )
                   else ! double precision
                      area1 = MISC_3Dvec_triangle( r0(:), r2(:), r3(:),       &
                                                   GMTR_polygon_type, rscale  )
                      area2 = MISC_3Dvec_triangle( r0(:), r3(:), r1(:),       &
                                                   GMTR_polygon_type, rscale  )
                      area3 = MISC_3Dvec_triangle( r0(:), r1(:), r2(:),       &
                                                   GMTR_polygon_type, rscale  )
                   endif

                   if (      area1 * 0.D0 /= 0.D0 &
                        .OR. area2 * 0.D0 /= 0.D0 &
                        .OR. area3 * 0.D0 /= 0.D0 ) then ! Nan?
                      write(*,          *) 'Nan! (i,j,n,t,l)=', i,j,n,t,l
                      write(*,          *) '(area1,area2,area3)=', area1,area2,area3
                      write(ADM_LOG_FID,*) 'Nan! (i,j,n,t,l)=', i,j,n,t,l
                      write(ADM_LOG_FID,*) '(area1,area2,area3)=', area1,area2,area3
                      call ADM_proc_stop
                   endif

                   area_total = area1 + area2 + area3

                   w1(nmax_llgrid) = area1 / area_total
                   w2(nmax_llgrid) = area2 / area_total
                   w3(nmax_llgrid) = area3 / area_total
                endselect

                cycle

             elseif(       t == ADM_TI               &
                     .AND. abs(v01(1)) <= eps_latlon &
                     .AND. abs(v01(2)) <= eps_latlon &
                     .AND. abs(v01(3)) <= eps_latlon ) then ! on the triangle vertex

                select case( trim(what_is_done) )
                case( 'GET_NUM' )

                   nmax_llgrid        = nmax_llgrid        + 1
                   nmax_llgrid_rgn(l) = nmax_llgrid_rgn(l) + 1

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
                   w1       (nmax_llgrid) = 1.D0
                   w2       (nmax_llgrid) = 0.D0
                   w3       (nmax_llgrid) = 0.D0
                endselect

             endif

          enddo ! i LOOP
          enddo ! j LOOP

       enddo ! TI,TJ

    enddo ! n LOOP
    enddo ! l LOOP

    return
  end subroutine mkrelmap_ico2ll

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
       latmax =  CNST_PI*0.5D0
       latmin = -CNST_PI*0.5D0
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

    real(8) :: dlat, dlon
    integer :: i, j
    !---------------------------------------------------------------------------

    dlat = ( latmax - latmin ) / real(jmax,kind=8)

    do j = 1, jmax
       lat(j) = latmin + dlat * ( real(j,kind=8) - 0.5D0 )
    enddo

    dlon = ( lonmax - lonmin ) / real(imax,kind=8)

    if ( lon_offset ) then
       do i = 1, imax
          lon(i) = lonmin + dlon * ( real(i,kind=8) - 0.5D0 )
       enddo
    else
       do i = 1, imax
          lon(i) = lonmin + dlon * ( real(i,kind=8) - 1.D0 )
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

    real(8) :: e(nb)
    real(8) :: eps

    real(8) :: mu0, dmu, dP0
    real(8) :: mu(jmax)
    real(8) :: P0(0:jmax)

    real(8) :: dlon

    integer :: n, i, j
    !---------------------------------------------------------------------------

    !--- calculate machine eps.
    e(:) = 0.D0
    eps  = 1.D0

    do i = 1, nb
       eps  = eps * 0.5D0
       e(i) = eps+1
    enddo

    i = 1
    eps = 1.D0

    do i = 1, nb
       if( e(i) > 1.D0 ) exit
       eps = eps * 0.5D0
    enddo
    eps = eps * 4.D0

    !---  calculate gausian_grid by Newton-Rapson method
    do j = 1, jmax

       mu0=sin(CNST_PI*dble(jmax+1-2*j)/dble(2*jmax+1))

       loopeps:do

          P0(0)=1.0d0
          P0(1)=mu0
          do n=1,jmax-1
             P0(n+1)=(dble(2*n+1)*mu0*P0(n)-n*P0(n-1))/dble(n+1)
          enddo
          dP0=jmax*(P0(jmax-1)-mu0*P0(jmax))/(1-mu0*mu0)
          dmu=P0(jmax)/dP0
          mu0=mu0-dmu
          if (abs(dmu)<eps) then
             mu(j)=mu0
             P0(0)=1.0d0
             P0(1)=mu(j)
             do n=1,jmax-1
                P0(n+1)=(dble(2*n+1)*mu(j)*P0(n)-n*P0(n-1))/dble(n+1)
             enddo
             exit loopeps
          endif

       enddo loopeps

    enddo

    do j = 1, jmax
       lat(j) = -asin(mu(j))
    enddo

    dlon = ( lonmax - lonmin) / real(imax,kind=8)

    if ( lon_offset ) then
       do i = 1, imax
          lon(i) = lonmin + dlon * ( real(i,kind=8) - 0.5D0 )
       enddo
    else
       do i = 1, imax
          lon(i) = lonmin + dlon * ( real(i,kind=8) - 1.D0 )
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
!         ADM_GALL_PL,              &
!         ADM_LALL_PL,              &
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
!    real(8), intent(out) :: var_ll(nmax_llgrid,kmin:kmax)
!    !    real(8), intent(in)  :: var(ADM_gall,kmin:kmax,ADM_lall)
!    !    real(8), intent(in)  :: var_pl(ADM_gall_pl,kmin:kmax,ADM_lall_pl)
!    !    real(8), intent(in)  :: var(:,:,:)
!    !    real(8), intent(in)  :: var_pl(:,:,:)
!    real(8), intent(in) :: var(ADM_gall,ADM_kall,ADM_lall)
!    real(8), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!    !
!    integer, parameter :: ix = 1
!    integer, parameter :: iy = 2
!    integer, parameter :: iz = 3
!    !
!    real(8) ::  grd(ADM_gall,ADM_kall,ADM_lall,ix:iz)
!    real(8) ::  grd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,ix:iz)
!    logical ::  def_grd(ADM_gall,kmin:kmax,ADM_lall)
!    !
!    integer :: n0,n1,n2,n3,n4,n5,n6
!    real(8) :: vec12(ix:iz),vec23(ix:iz),vec31(ix:iz)
!    real(8) :: v1,v2,v3,v12,v23,v31
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
!          v12 = (v1+v2)*0.5D0                                   &
!               + 0.125D0                                        &
!               * (( grd(n1,k,l,ix) - grd(n2,k,l,ix) )*vec12(ix) &
!                 +( grd(n1,k,l,iy) - grd(n2,k,l,iy) )*vec12(iy) &
!                 +( grd(n1,k,l,iz) - grd(n2,k,l,iz) )*vec12(iz))
!          !
!          !--- midpoint value between n2 and n3
!          !--- based on Hermitian interpolation.
!          v23 = (v2+v3)*0.5D0                                   &
!               + 0.125D0                                        &
!               * (( grd(n2,k,l,ix) - grd(n3,k,l,ix) )*vec23(ix) &
!                 +( grd(n2,k,l,iy) - grd(n3,k,l,iy) )*vec23(iy) &
!                 +( grd(n2,k,l,iz) - grd(n3,k,l,iz) )*vec23(iz))
!          !
!          !--- midpoint value between n3 and n1
!          !--- based on Hermitian interpolation.
!          v31 = (v3+v1)*0.5D0                                   &
!               + 0.125D0                                        &
!               * (( grd(n3,k,l,ix) - grd(n1,k,l,ix) )*vec31(ix) &
!                 +( grd(n3,k,l,iy) - grd(n1,k,l,iy) )*vec31(iy) &
!                 +( grd(n3,k,l,iz) - grd(n1,k,l,iz) )*vec31(iz))
!          if(  def_grd(n1,k,l).and. &
!               def_grd(n2,k,l).and. &
!               def_grd(n3,k,l) ) then
!             var_ll(m,k)                          &
!                  = w1(m)*(2.0D0*w1(m)-1.0D0)*v1  &
!                  + w2(m)*(2.0D0*w2(m)-1.0D0)*v2  &
!                  + w3(m)*(2.0D0*w3(m)-1.0D0)*v3  &
!                  + 4.0D0*w1(m)*w2(m)        *v12 &
!                  + 4.0D0*w2(m)*w3(m)        *v23 &
!                  + 4.0D0*w3(m)*w1(m)        *v31
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
