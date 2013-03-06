!-------------------------------------------------------------------------------
!>
!! Geometrics module
!!
!! @par Description
!!          In this module, the geometrics of the icosahedral grid such as
!!          area are calculated.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2010-06-08 (S.Iga)     a new grid is implemented
!! @li      2011-07-22 (T.Ohno)    a new grid is implemented
!! @li      2011-08-18 (T.Ohno)    bugfix
!!
!<
module mod_gmtr
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GMTR_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: GMTR_T_nmax_var = 7

  integer, public, parameter :: GMTR_T_AREA  = 1
  integer, public, parameter :: GMTR_T_RAREA = 2
  integer, public, parameter :: GMTR_T_W1    = 3
  integer, public, parameter :: GMTR_T_W2    = 4
  integer, public, parameter :: GMTR_T_W3    = 5
  integer, public, parameter :: GMTR_T_LAT   = 6
  integer, public, parameter :: GMTR_T_LON   = 7

  integer, public, parameter :: GMTR_P_nmax_var = 10

  integer, public, parameter :: GMTR_P_AREA  = 1
  integer, public, parameter :: GMTR_P_RAREA = 2
  integer, public, parameter :: GMTR_P_IX    = 3
  integer, public, parameter :: GMTR_P_IY    = 4
  integer, public, parameter :: GMTR_P_IZ    = 5
  integer, public, parameter :: GMTR_P_JX    = 6
  integer, public, parameter :: GMTR_P_JY    = 7
  integer, public, parameter :: GMTR_P_JZ    = 8
  integer, public, parameter :: GMTR_P_LAT   = 9
  integer, public, parameter :: GMTR_P_LON   = 10

  integer, public, parameter :: GMTR_A_nmax_var    = 12
  integer, public, parameter :: GMTR_A_nmax_var_pl = 18

  integer, public, parameter :: GMTR_A_HNX  = 1
  integer, public, parameter :: GMTR_A_HNY  = 2
  integer, public, parameter :: GMTR_A_HNZ  = 3
  integer, public, parameter :: GMTR_A_HTX  = 4
  integer, public, parameter :: GMTR_A_HTY  = 5
  integer, public, parameter :: GMTR_A_HTZ  = 6
  integer, public, parameter :: GMTR_A_TNX  = 7
  integer, public, parameter :: GMTR_A_TNY  = 8
  integer, public, parameter :: GMTR_A_TNZ  = 9
  integer, public, parameter :: GMTR_A_TTX  = 10
  integer, public, parameter :: GMTR_A_TTY  = 11
  integer, public, parameter :: GMTR_A_TTZ  = 12

  integer, public, parameter :: GMTR_A_TN2X = 13
  integer, public, parameter :: GMTR_A_TN2Y = 14
  integer, public, parameter :: GMTR_A_TN2Z = 15
  integer, public, parameter :: GMTR_A_TT2X = 16
  integer, public, parameter :: GMTR_A_TT2Y = 17
  integer, public, parameter :: GMTR_A_TT2Z = 18

  real(8), public, allocatable, save :: GMTR_P_var   (:,:,:,:)
  real(8), public, allocatable, save :: GMTR_P_var_pl(:,:,:,:)

  real(8), public, allocatable, save :: GMTR_T_var   (:,:,:,:,:)
  real(8), public, allocatable, save :: GMTR_T_var_pl(:,:,:,:)

  real(8), public, allocatable, save :: GMTR_A_var   (:,:,:,:,:)
  real(8), public, allocatable, save :: GMTR_A_var_pl(:,:,:,:)

  real(8), public, allocatable, save :: GMTR_area   (:,:)
  real(8), public, allocatable, save :: GMTR_area_pl(:,:)
  real(8), public, allocatable, save :: GMTR_lat    (:,:)
  real(8), public, allocatable, save :: GMTR_lat_pl (:,:)
  real(8), public, allocatable, save :: GMTR_lon    (:,:)
  real(8), public, allocatable, save :: GMTR_lon_pl (:,:)

  character(len=ADM_NSYS),  public, save :: GMTR_polygon_type = 'ON_SPHERE'
  !                                       ! 'ON_SPHERE' : triangle is fit to the sphere
  !                                       ! 'ON_PLANE'  : triangle is treated as 2D

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: xcalc_gmtr_t
  private :: xcalc_gmtr_p
  private :: xcalc_gmtr_a

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine GMTR_setup
  !>
  subroutine GMTR_setup
    use mod_adm, only: &
       ADM_CTL_FID,      &
       ADM_LOG_FID,      &
       ADM_proc_stop,    &
       ADM_TI,           &
       ADM_TJ,           &
       ADM_AI,           &
       ADM_AJ,           &
       ADM_gmin,         &
       ADM_gmax,         &
       ADM_gall,         &
       ADM_gall_1d,      &
       ADM_gall_pl,      &
       ADM_lall,         &
       ADM_lall_pl,      &
       ADM_KNONE
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    character(len=ADM_NSYS) :: polygon_type !--- polygon type

    namelist / GMTRPARAM / &
        polygon_type

    integer :: ierr
    integer :: K0
    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    polygon_type = GMTR_polygon_type

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[gmtr]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=GMTRPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** GMTRPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist GMTRPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist GMTRPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,GMTRPARAM)

    GMTR_polygon_type = polygon_type

    K0 = ADM_KNONE

    ! --- setup triangle data
    allocate( GMTR_T_var   (ADM_gall,   K0,ADM_lall,   ADM_TI:ADM_TJ,GMTR_T_nmax_var) )
    allocate( GMTR_T_var_pl(ADM_gall_pl,K0,ADM_lall_pl,              GMTR_T_nmax_var) )
    GMTR_T_var   (:,:,:,:,:) = 0.D0
    GMTR_T_var_pl(:,:,:,:)   = 0.D0

    call xcalc_gmtr_t

    !--- setup point data
    allocate( GMTR_P_var   (ADM_gall,   K0,ADM_lall,   GMTR_P_nmax_var) )
    allocate( GMTR_P_var_pl(ADM_gall_pl,K0,ADM_lall_pl,GMTR_P_nmax_var) )
    GMTR_P_var   (:,:,:,:) = 0.D0
    GMTR_P_var_pl(:,:,:,:) = 0.D0

    call xcalc_gmtr_p

    !--- communication of point data
    call COMM_data_transfer( GMTR_P_var, GMTR_P_var_pl )
    ! fill unused grid (dummy)
    GMTR_P_var(suf(ADM_gall_1d,1),:,:,:) = GMTR_P_var(suf(ADM_gmax+1,ADM_gmin),:,:,:)
    GMTR_P_var(suf(1,ADM_gall_1d),:,:,:) = GMTR_P_var(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    !--- for simple use
    allocate( GMTR_area   (ADM_gall,   ADM_lall   ) )
    allocate( GMTR_area_pl(ADM_gall_pl,ADM_lall_pl) )
    allocate( GMTR_lat    (ADM_gall,   ADM_lall   ) )
    allocate( GMTR_lat_pl (ADM_gall_pl,ADM_lall_pl) )
    allocate( GMTR_lon    (ADM_gall,   ADM_lall   ) )
    allocate( GMTR_lon_pl (ADM_gall_pl,ADM_lall_pl) )

    GMTR_area   (:,:) = GMTR_P_var   (:,K0,:,GMTR_P_AREA)
    GMTR_area_pl(:,:) = GMTR_P_var_pl(:,K0,:,GMTR_P_AREA)
    GMTR_lat    (:,:) = GMTR_P_var   (:,K0,:,GMTR_P_LAT )
    GMTR_lat_pl (:,:) = GMTR_P_var_pl(:,K0,:,GMTR_P_LAT )
    GMTR_lon    (:,:) = GMTR_P_var   (:,K0,:,GMTR_P_LON )
    GMTR_lon_pl (:,:) = GMTR_P_var_pl(:,K0,:,GMTR_P_LON )

    !--- setup arc data
    allocate( GMTR_A_var   (ADM_gall,   K0,ADM_lall,   ADM_AI:ADM_AJ,GMTR_A_nmax_var   ) )
    allocate( GMTR_A_var_pl(ADM_gall_pl,K0,ADM_lall_pl,              GMTR_A_nmax_var_pl) )
    GMTR_A_var   (:,:,:,:,:) = 0.D0
    GMTR_A_var_pl(:,:,:,:)   = 0.D0

    call xcalc_gmtr_a

    return
  end subroutine GMTR_setup

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine xcalc_gmtr_t
  !>
  subroutine xcalc_gmtr_t
    use mod_misc, only: &
       MISC_triangle_area, &
       MISC_get_latlon
    use mod_adm, only: &
       ADM_prc_me,      &
       ADM_prc_pl,      &
       ADM_prc_tab,     &
       ADM_rgn_vnum,    &
       ADM_W,           &
       ADM_TI,          &
       ADM_TJ,          &
       ADM_gmin,        &
       ADM_gmax,        &
       ADM_gall,        &
       ADM_gall_1d,     &
       ADM_gall_pl,     &
       ADM_GSLF_PL,     &
       ADM_GMIN_PL,     &
       ADM_GMAX_PL,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_KNONE,       &
       ADM_ImoJmo_nmax, &
       ADM_ImoJmo,      &
       ADM_GIoJo,       &
       ADM_GIpJo,       &
       ADM_GIpJp,       &
       ADM_GIoJp
    use mod_grd, only: &
       GRD_XDIR,      &
       GRD_YDIR,      &
       GRD_ZDIR,      &
       GRD_x,         &
       GRD_x_pl,      &
       GRD_xt,        &
       GRD_xt_pl,     &
       GRD_grid_type, & ! [add] T.Ohno 110722
       GRD_rscale
    implicit none

    real(8) :: v   (GRD_XDIR:GRD_ZDIR,0:3,ADM_gall,   ADM_TI:ADM_TJ)
    real(8) :: v_pl(GRD_XDIR:GRD_ZDIR,0:3,ADM_gall_pl)

    real(8) :: area, area1, area2, area3
    integer :: l, d, t, n
    integer :: rgnid, ij, K0

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    do l = 1,ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          v(GRD_XDIR,0,ij,ADM_TI) = GRD_xt(ij,K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,1,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,3,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_XDIR)

          v(GRD_XDIR,0,ij,ADM_TJ) = GRD_xt(ij,K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,1,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_XDIR)
          v(GRD_XDIR,3,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJp),K0,l,GRD_XDIR)

          v(GRD_YDIR,0,ij,ADM_TI) = GRD_xt(ij,K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,1,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,3,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_YDIR)

          v(GRD_YDIR,0,ij,ADM_TJ) = GRD_xt(ij,K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,1,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_YDIR)
          v(GRD_YDIR,3,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJp),K0,l,GRD_YDIR)

          v(GRD_ZDIR,0,ij,ADM_TI) = GRD_xt(ij,K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,1,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,3,ij,ADM_TI) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_ZDIR)

          v(GRD_ZDIR,0,ij,ADM_TJ) = GRD_xt(ij,K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,1,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,3,ij,ADM_TJ) = GRD_x(ADM_ImoJmo(n,ADM_GIoJp),K0,l,GRD_ZDIR)
       enddo

       !--- treat unused point
       v(:,:,suf(ADM_gmax,ADM_gmin-1),ADM_TI) = v(:,:,suf(ADM_gmax,ADM_gmin-1),ADM_TJ)
       v(:,:,suf(ADM_gmin-1,ADM_gmax),ADM_TJ) = v(:,:,suf(ADM_gmin-1,ADM_gmax),ADM_TI)

       !--- exception for the west
       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          v(:,:,suf(ADM_gmin-1,ADM_gmin-1),ADM_TI) = v(:,:,suf(ADM_gmin,ADM_gmin-1),ADM_TJ)
       endif

       do t = ADM_TI,ADM_TJ
       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722
             area1 = get_triangle_are_on_plane( v(:,0,ij,t), v(:,2,ij,t), v(:,3,ij,t) )
             area2 = get_triangle_are_on_plane( v(:,0,ij,t), v(:,3,ij,t), v(:,1,ij,t) )
             area3 = get_triangle_are_on_plane( v(:,0,ij,t), v(:,1,ij,t), v(:,2,ij,t) )
          else
             area1 = MISC_triangle_area( v(:,0,ij,t), v(:,2,ij,t), v(:,3,ij,t), GMTR_polygon_type, GRD_rscale )
             area2 = MISC_triangle_area( v(:,0,ij,t), v(:,3,ij,t), v(:,1,ij,t), GMTR_polygon_type, GRD_rscale )
             area3 = MISC_triangle_area( v(:,0,ij,t), v(:,1,ij,t), v(:,2,ij,t), GMTR_polygon_type, GRD_rscale )
          endif

          area = area1 + area2 + area3

          GMTR_T_var(ij,K0,l,t,GMTR_T_AREA)  = area
          GMTR_T_var(ij,K0,l,t,GMTR_T_RAREA) = 1.D0 / area

          GMTR_T_var(ij,K0,l,t,GMTR_T_W1)    = area1 / area
          GMTR_T_var(ij,K0,l,t,GMTR_T_W2)    = area2 / area
          GMTR_T_var(ij,K0,l,t,GMTR_T_W3)    = area3 / area

          call MISC_get_latlon( GMTR_T_var(ij,K0,l,t,GMTR_T_LAT), &
                                GMTR_T_var(ij,K0,l,t,GMTR_T_LON), &
                                GRD_xt    (ij,K0,l,t,GRD_XDIR),   &
                                GRD_xt    (ij,K0,l,t,GRD_YDIR),   &
                                GRD_xt    (ij,K0,l,t,GRD_ZDIR)    )
       enddo
       enddo

    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1,ADM_lall_pl

          do n = ADM_GMIN_PL, ADM_GMAX_PL-1
             do d = GRD_XDIR, GRD_ZDIR
                v_pl(d,0,n) = GRD_xt_pl(n,K0,l,d)
                v_pl(d,1,n) = GRD_x_pl(ADM_GSLF_PL,K0,l,d)
                v_pl(d,2,n) = GRD_x_pl(n  ,        K0,l,d)
                v_pl(d,3,n) = GRD_x_pl(n+1,        K0,l,d)
             enddo
          enddo

          n = ADM_GMAX_PL
          do d = GRD_XDIR, GRD_ZDIR
             v_pl(d,0,n) = GRD_xt_pl(n,K0,l,d)
             v_pl(d,1,n) = GRD_x_pl(ADM_GSLF_PL,K0,l,d)
             v_pl(d,2,n) = GRD_x_pl(n,          K0,l,d)
             v_pl(d,3,n) = GRD_x_pl(ADM_GMIN_PL,K0,l,d)
          enddo

          do n = ADM_GMIN_PL, ADM_GMAX_PL
             area1 = MISC_triangle_area( v_pl(:,0,n), v_pl(:,2,n), v_pl(:,3,n), GMTR_polygon_type, GRD_rscale )
             area2 = MISC_triangle_area( v_pl(:,0,n), v_pl(:,3,n), v_pl(:,1,n), GMTR_polygon_type, GRD_rscale )
             area3 = MISC_triangle_area( v_pl(:,0,n), v_pl(:,1,n), v_pl(:,2,n), GMTR_polygon_type, GRD_rscale )

             area = area1 + area2 + area3

             GMTR_T_var_pl(n,K0,l,GMTR_T_AREA)  = area
             GMTR_T_var_pl(n,K0,l,GMTR_T_RAREA) = 1.D0 / area

             GMTR_T_var_pl(n,K0,l,GMTR_T_W1)    = area1 / area
             GMTR_T_var_pl(n,K0,l,GMTR_T_W2)    = area2 / area
             GMTR_T_var_pl(n,K0,l,GMTR_T_W3)    = area3 / area


             call MISC_get_latlon( GMTR_T_var_pl(n,K0,l,GMTR_T_LAT), &
                                   GMTR_T_var_pl(n,K0,l,GMTR_T_LON), &
                                   GRD_xt_pl    (n,K0,l,GRD_XDIR),   &
                                   GRD_xt_pl    (n,K0,l,GRD_YDIR),   &
                                   GRD_xt_pl    (n,K0,l,GRD_ZDIR)    )
          enddo

       enddo
    endif

    return
  end subroutine xcalc_gmtr_t

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine xcalc_gmtr_p
  !>
  subroutine xcalc_gmtr_p
    use mod_misc, only: &
       MISC_triangle_area, &
       MISC_get_latlon
    use mod_adm, only: &
       ADM_prc_me,      &
       ADM_prc_pl,      &
       ADM_prc_tab,     &
       ADM_rgn_vnum,    &
       ADM_W,           &
       ADM_TI,          &
       ADM_TJ,          &
       ADM_gmin,        &
       ADM_gall,        &
       ADM_gall_1d,     &
       ADM_GSLF_PL,     &
       ADM_GMIN_PL,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_KNONE,       &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GIoJm,       &
       ADM_GImJo,       &
       ADM_GImJm,       &
       ADM_XTMS_K         ! S.Iga100608
    use mod_grd, only: &
       GRD_XDIR,      &
       GRD_YDIR,      &
       GRD_ZDIR,      &
       GRD_x,         &
       GRD_x_pl,      &
       GRD_xt,        &
       GRD_xt_pl,     &
       GRD_grid_type, & ! [add] T.Ohno 110722
       GRD_rscale
    implicit none

    real(8) :: v   (GRD_XDIR:GRD_ZDIR,0:7,ADM_gall)
    real(8) :: v_pl(GRD_XDIR:GRD_ZDIR,0:ADM_XTMS_K+1)

    real(8) :: area
    real(8) :: cos_lam, sin_lam

    integer :: l, n, m
    integer :: rgnid, ij, K0

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do n = 1, ADM_IooJoo_nmax
          ij = ADM_IooJoo(n,ADM_GIoJo)

          v(GRD_XDIR,0,ij) = GRD_x(ij,K0,l,GRD_XDIR)
          v(GRD_XDIR,1,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,3,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,4,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJo),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,5,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,6,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,7,ij) = v(GRD_XDIR,1,ij)

          v(GRD_YDIR,0,ij) = GRD_x(ij,K0,l,GRD_YDIR)
          v(GRD_YDIR,1,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,3,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,4,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJo),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,5,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,6,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,7,ij) = v(GRD_YDIR,1,ij)

          v(GRD_ZDIR,0,ij) = GRD_x(ij,K0,l,GRD_ZDIR)
          v(GRD_ZDIR,1,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,3,ij) = GRD_xt(ADM_IooJoo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,4,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJo),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,5,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,6,ij) = GRD_xt(ADM_IooJoo(n,ADM_GImJm),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,7,ij) = v(GRD_ZDIR,1,ij)
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          v(:,6,suf(ADM_gmin,ADM_gmin)) = v(:,1,suf(ADM_gmin,ADM_gmin))
          v(:,7,suf(ADM_gmin,ADM_gmin)) = v(:,1,suf(ADM_gmin,ADM_gmin))
       endif

       do n = 1, ADM_IooJoo_nmax
          ij = ADM_IooJoo(n,ADM_GIoJo)

          area = 0.D0
          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722
             do m = 1, 6
                area = area + get_triangle_are_on_plane( v(:,0,ij), v(:,m,ij), v(:,m+1,ij) )
             enddo
          else
             do m = 1, 6
                area = area + MISC_triangle_area( v(:,0,ij), v(:,m,ij), v(:,m+1,ij), &
                                                  GMTR_polygon_type, GRD_rscale      )
             enddo
          endif

          GMTR_P_var(ij,K0,l,GMTR_P_AREA)  = area
          GMTR_P_var(ij,K0,l,GMTR_P_RAREA) = 1.D0 / GMTR_P_var(ij,K0,l,GMTR_P_AREA)

          call MISC_get_latlon( GMTR_P_var(ij,K0,l,GMTR_P_LAT), &
                                GMTR_P_var(ij,K0,l,GMTR_P_LON), &
                                GRD_x     (ij,K0,l,GRD_XDIR),   &
                                GRD_x     (ij,K0,l,GRD_YDIR),   &
                                GRD_x     (ij,K0,l,GRD_ZDIR)    )

          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722

             GMTR_P_var(ij,K0,l,GMTR_P_IX) = 1.D0
             GMTR_P_var(ij,K0,l,GMTR_P_IY) = 0.D0
             GMTR_P_var(ij,K0,l,GMTR_P_IZ) = 0.D0
             GMTR_P_var(ij,K0,l,GMTR_P_JX) = 0.D0
             GMTR_P_var(ij,K0,l,GMTR_P_JY) = 1.D0
             GMTR_P_var(ij,K0,l,GMTR_P_JZ) = 0.D0

          else

             sin_lam = sin( GMTR_P_var(ij,K0,l,GMTR_P_LON) )
             cos_lam = cos( GMTR_P_var(ij,K0,l,GMTR_P_LON) )

             GMTR_P_var(ij,K0,l,GMTR_P_IX) = -sin_lam
             GMTR_P_var(ij,K0,l,GMTR_P_IY) =  cos_lam
             GMTR_P_var(ij,K0,l,GMTR_P_IZ) = 0.0D0
             GMTR_P_var(ij,K0,l,GMTR_P_JX) = -( GRD_x(ij,K0,l,GRD_ZDIR) * cos_lam ) / GRD_rscale
             GMTR_P_var(ij,K0,l,GMTR_P_JY) = -( GRD_x(ij,K0,l,GRD_ZDIR) * sin_lam ) / GRD_rscale
             GMTR_P_var(ij,K0,l,GMTR_P_JZ) =  ( GRD_x(ij,K0,l,GRD_XDIR) * cos_lam &
                                              + GRD_x(ij,K0,l,GRD_YDIR) * sin_lam ) / GRD_rscale
          endif

       enddo ! ij loop
    enddo ! l loop

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_GSLF_PL

       do l = 1,ADM_lall_pl

          v_pl(GRD_XDIR,0) = GRD_x_pl(n,K0,l,GRD_XDIR)
          v_pl(GRD_YDIR,0) = GRD_x_pl(n,K0,l,GRD_YDIR)
          v_pl(GRD_ZDIR,0) = GRD_x_pl(n,K0,l,GRD_ZDIR)
          do m = 1, ADM_XTMS_K ! (ICO=5)
             v_pl(GRD_XDIR,m) = GRD_xt_pl(m+ADM_GMIN_PL-1,K0,l,GRD_XDIR)
             v_pl(GRD_YDIR,m) = GRD_xt_pl(m+ADM_GMIN_PL-1,K0,l,GRD_YDIR)
             v_pl(GRD_ZDIR,m) = GRD_xt_pl(m+ADM_GMIN_PL-1,K0,l,GRD_ZDIR)
          enddo
          v_pl(:,ADM_XTMS_K+1) = v_pl(:,1)

          area = 0.D0
          do m = 1, ADM_XTMS_K ! (ICO=5)
             area = area + MISC_triangle_area( v_pl(:,0), v_pl(:,m), v_pl(:,m+1), &
                                               GMTR_polygon_type, GRD_rscale      )
          enddo

          GMTR_P_var_pl(n,K0,l,GMTR_P_AREA)  = area
          GMTR_P_var_pl(n,K0,l,GMTR_P_RAREA) = 1.D0 / GMTR_P_var_pl(n,K0,l,GMTR_P_AREA)

          call MISC_get_latlon( GMTR_P_var_pl(n,K0,l,GMTR_P_LAT), &
                                GMTR_P_var_pl(n,K0,l,GMTR_P_LON), &
                                GRD_x_pl     (n,K0,l,GRD_XDIR),   &
                                GRD_x_pl     (n,K0,l,GRD_YDIR),   &
                                GRD_x_pl     (n,K0,l,GRD_ZDIR)    )

          sin_lam = sin( GMTR_P_var_pl(n,K0,l,GMTR_P_LON) )
          cos_lam = cos( GMTR_P_var_pl(n,K0,l,GMTR_P_LON) )

          GMTR_P_var_pl(n,K0,l,GMTR_P_IX) = -sin_lam
          GMTR_P_var_pl(n,K0,l,GMTR_P_IY) =  cos_lam
          GMTR_P_var_pl(n,K0,l,GMTR_P_IZ) = 0.0D0
          GMTR_P_var_pl(n,K0,l,GMTR_P_JX) = -( GRD_x_pl(n,K0,l,GRD_ZDIR) * cos_lam ) / GRD_rscale
          GMTR_P_var_pl(n,K0,l,GMTR_P_JY) = -( GRD_x_pl(n,K0,l,GRD_ZDIR) * sin_lam ) / GRD_rscale
          GMTR_P_var_pl(n,K0,l,GMTR_P_JZ) =  ( GRD_x_pl(n,K0,l,GRD_XDIR) * cos_lam &
                                             + GRD_x_pl(n,K0,l,GRD_YDIR) * sin_lam ) / GRD_rscale
       enddo ! l loop
    endif

    return
  end subroutine xcalc_gmtr_p

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine xcalc_gmtr_a
  !>
  subroutine xcalc_gmtr_a
    use mod_misc, only: &
       MISC_mk_gmtrvec
    use mod_adm, only: &
       ADM_prc_me,      &
       ADM_prc_pl,      &
       ADM_prc_tab,     &
       ADM_rgn_vnum,    &
       ADM_W,           &
       ADM_TI,          &
       ADM_TJ,          &
       ADM_AI,          &
       ADM_AIJ,         &
       ADM_AJ,          &
       ADM_gmin,        &
       ADM_gmax,        &
       ADM_gall,        &
       ADM_gall_1d,     &
       ADM_gall_pl,     &
       ADM_GSLF_PL,     &
       ADM_GMIN_PL,     &
       ADM_GMAX_PL,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_KNONE,       &
       ADM_ImpJmo_nmax, &
       ADM_ImoJoo_nmax, &
       ADM_IooJmo_nmax, &
       ADM_ImoJmp_nmax, &
       ADM_ImoJmo_nmax, &
       ADM_ImpJmo,      &
       ADM_ImoJoo,      &
       ADM_IooJmo,      &
       ADM_ImoJmo,      &
       ADM_ImoJmp,      &
       ADM_GIoJo,       &
       ADM_GIpJo,       &
       ADM_GIpJp,       &
       ADM_GIoJp,       &
       ADM_GIoJm,       &
       ADM_GImJo
    use mod_grd, only: &
       GRD_XDIR,      &
       GRD_YDIR,      &
       GRD_ZDIR,      &
       GRD_x,         &
       GRD_x_pl,      &
       GRD_xt,        &
       GRD_xt_pl,     &
       GRD_grid_type, & ! [add] T.Ohno 110722
       GRD_rscale
    implicit none

    real(8) :: v   (GRD_XDIR:GRD_ZDIR,2,ADM_gall   )
    real(8) :: v_pl(GRD_XDIR:GRD_ZDIR,2,ADM_gall_pl)

    real(8) :: tvec(GRD_XDIR:GRD_ZDIR)
    real(8) :: nvec(GRD_XDIR:GRD_ZDIR)

    integer :: ij, K0, l, d
    integer :: rgnid
    integer :: n

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    !--- Triangle
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       !--- AI
       do n = 1, ADM_ImoJmp_nmax
          ij = ADM_ImoJmp(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIpJo),K0,l,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIpJo),K0,l,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_x(ADM_ImoJmp(n,ADM_GIpJo),K0,l,GRD_ZDIR)
       enddo

       do d = GRD_XDIR, GRD_ZDIR
          !--- execetion for the south.
          v(d,1,suf(ADM_gmax,ADM_gmin-1)) = GRD_x(suf(ADM_gmax,ADM_gmin-1),K0,l,d)
          v(d,2,suf(ADM_gmax,ADM_gmin-1)) = GRD_x(suf(ADM_gmax,ADM_gmin),  K0,l,d)

          !--- execetion for the south.
          v(d,1,suf(ADM_gmin-1,ADM_gmax+1)) = GRD_x(suf(ADM_gmin,ADM_gmax+1),K0,l,d)
          v(d,2,suf(ADM_gmin-1,ADM_gmax+1)) = GRD_x(suf(ADM_gmin,ADM_gmax),  K0,l,d)

          !--- exception for the west
          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
             v(d,1,suf(ADM_gmin-1,ADM_gmin-1)) = GRD_x(suf(ADM_gmin,ADM_gmin-1),K0,l,d)
             v(d,2,suf(ADM_gmin-1,ADM_gmin-1)) = GRD_x(suf(ADM_gmin+1,ADM_gmin),K0,l,d)
          endif
       enddo

       do n = 1, ADM_ImoJmp_nmax
          ij = ADM_ImoJmp(n,ADM_GIoJo)

          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call MISC_mk_gmtrvec( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_TNZ) = nvec(3)
       enddo

       !--- AIJ
       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_x(ADM_ImoJmo(n,ADM_GIpJp),K0,l,GRD_ZDIR)
       enddo

       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call MISC_mk_gmtrvec( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_TNZ) = nvec(3)
       enddo

       !--- AJ
       do n = 1, ADM_ImpJmo_nmax
          ij = ADM_ImpJmo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJo),K0,l,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJp),K0,l,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJo),K0,l,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJp),K0,l,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJo),K0,l,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_x(ADM_ImpJmo(n,ADM_GIoJp),K0,l,GRD_ZDIR)
       enddo

       do d = GRD_XDIR, GRD_ZDIR
          !--- execetion for the south.
          v(d,1,suf(ADM_gmax+1,ADM_gmin-1)) = GRD_x(suf(ADM_gmax+1,ADM_gmin),K0,l,d)
          v(d,2,suf(ADM_gmax+1,ADM_gmin-1)) = GRD_x(suf(ADM_gmax,ADM_gmin),  K0,l,d)

          !--- execetion for the north.
          v(d,1,suf(ADM_gmin-1,ADM_gmax)) = GRD_x(suf(ADM_gmin-1,ADM_gmax),K0,l,d)
          v(d,2,suf(ADM_gmin-1,ADM_gmax)) = GRD_x(suf(ADM_gmin,ADM_gmax),  K0,l,d)
       enddo

       do n = 1, ADM_ImpJmo_nmax
          ij = ADM_ImpJmo(n,ADM_GIoJo)

          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call MISC_mk_gmtrvec( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_TNZ) = nvec(3)
       enddo
    enddo

    !
    ! --- Hexagon
    !
    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       !--- AI
       do n = 1, ADM_ImoJoo_nmax
          ij = ADM_ImoJoo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_xt(ADM_ImoJoo(n,ADM_GIoJm),K0,l,ADM_TJ,GRD_ZDIR)
       enddo

       do n = 1, ADM_ImoJoo_nmax
          ij = ADM_ImoJoo(n,ADM_GIoJo)

          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call MISC_mk_gmtrvec( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AI,GMTR_A_HNZ) = nvec(3)
       enddo

       !--- AIJ
       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_xt(ADM_ImoJmo(n,ADM_GIoJo),K0,l,ADM_TI,GRD_ZDIR)
       enddo

       do d = GRD_XDIR, GRD_ZDIR
          !--- execetion for the south.
          v(d,1,suf(ADM_gmax,ADM_gmin-1)) = GRD_xt(suf(ADM_gmax,ADM_gmin-1),K0,l,ADM_TJ,d)
          v(d,2,suf(ADM_gmax,ADM_gmin-1)) = GRD_xt(suf(ADM_gmax,ADM_gmin),  K0,l,ADM_TI,d)

          !--- execetion for the north.
          v(d,1,suf(ADM_gmin-1,ADM_gmax)) = GRD_xt(suf(ADM_gmin,ADM_gmax),  K0,l,ADM_TJ,d)
          v(d,2,suf(ADM_gmin-1,ADM_gmax)) = GRD_xt(suf(ADM_gmin-1,ADM_gmax),K0,l,ADM_TI,d)
       enddo

       do n = 1, ADM_ImoJmo_nmax
          ij = ADM_ImoJmo(n,ADM_GIoJo)

          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call MISC_mk_gmtrvec( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AIJ,GMTR_A_HNZ) = nvec(3)
       enddo

       !--- AJ
       do n = 1, ADM_IooJmo_nmax
          ij = ADM_IooJmo(n,ADM_GIoJo)

          v(GRD_XDIR,1,ij) = GRD_xt(ADM_IooJmo(n,ADM_GImJo),K0,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,2,ij) = GRD_xt(ADM_IooJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_XDIR)

          v(GRD_YDIR,1,ij) = GRD_xt(ADM_IooJmo(n,ADM_GImJo),K0,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,2,ij) = GRD_xt(ADM_IooJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_YDIR)

          v(GRD_ZDIR,1,ij) = GRD_xt(ADM_IooJmo(n,ADM_GImJo),K0,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,2,ij) = GRD_xt(ADM_IooJmo(n,ADM_GIoJo),K0,l,ADM_TJ,GRD_ZDIR)
       enddo

       !--- exception for the west
       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          do d = GRD_XDIR, GRD_ZDIR
             v(d,1,suf(ADM_gmin,ADM_gmin-1)) = GRD_xt(suf(ADM_gmin,ADM_gmin),  K0,l,ADM_TI,d)
             v(d,2,suf(ADM_gmin,ADM_gmin-1)) = GRD_xt(suf(ADM_gmin,ADM_gmin-1),K0,l,ADM_TJ,d)
          enddo
       endif

       do n = 1, ADM_IooJmo_nmax
          ij = ADM_IooJmo(n,ADM_GIoJo)

          if ( trim(GRD_grid_type) == 'ON_PLANE' ) then ! [add] T.Ohno 20110722
             call mk_gmtrvec_on_plane( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:) )
          else
             call MISC_mk_gmtrvec( v(:,1,ij), v(:,2,ij), tvec(:), nvec(:), GMTR_polygon_type, GRD_rscale )
          endif

          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HTX) = tvec(1)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HTY) = tvec(2)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HTZ) = tvec(3)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HNX) = nvec(1)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HNY) = nvec(2)
          GMTR_A_var(ij,K0,l,ADM_AJ,GMTR_A_HNZ) = nvec(3)
       enddo

    enddo ! l loop

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          do ij = ADM_GMIN_PL, ADM_GMAX_PL
             do d = GRD_XDIR, GRD_ZDIR
                v_pl(d,1,ij) = GRD_x_pl(ADM_GSLF_PL,K0,l,d)
                v_pl(d,2,ij) = GRD_x_pl(ij         ,K0,l,d)
             enddo

             call MISC_mk_gmtrvec( v_pl(:,1,ij), v_pl(:,2,ij), tvec(:), nvec(:), &
                                   GMTR_polygon_type, GRD_rscale                 )

             GMTR_A_var_pl(ij,K0,l,GMTR_A_TTX:GMTR_A_TTZ) = tvec(1:3)
             GMTR_A_var_pl(ij,K0,l,GMTR_A_TNX:GMTR_A_TNZ) = nvec(1:3)
          enddo

          do ij = ADM_GMIN_PL, ADM_GMAX_PL-1
             do d = GRD_XDIR, GRD_ZDIR
                v_pl(d,1,ij) = GRD_x_pl(ij,  K0,l,d)
                v_pl(d,2,ij) = GRD_x_pl(ij+1,K0,l,d)
             enddo
          enddo
          do d = GRD_XDIR, GRD_ZDIR
             v_pl(d,1,ADM_GMAX_PL) = GRD_x_pl(ADM_GMAX_PL,K0,l,d)
             v_pl(d,2,ADM_GMAX_PL) = GRD_x_pl(ADM_GMIN_PL,K0,l,d)
          enddo

          do ij = ADM_GMIN_PL, ADM_GMAX_PL
             call MISC_mk_gmtrvec( v_pl(:,1,ij), v_pl(:,2,ij), tvec(:), nvec(:), &
                                   GMTR_polygon_type, GRD_rscale                 )

             GMTR_A_var_pl(ij,K0,l,GMTR_A_TT2X:GMTR_A_TT2Z) = tvec(1:3)
             GMTR_A_var_pl(ij,K0,l,GMTR_A_TN2X:GMTR_A_TN2Z) = nvec(1:3)
          enddo

          do d = GRD_XDIR, GRD_ZDIR
             v_pl(d,1,ADM_GMIN_PL) = GRD_xt_pl(ADM_GMAX_PL,K0,l,d)
             v_pl(d,2,ADM_GMIN_PL) = GRD_xt_pl(ADM_GMIN_PL,K0,l,d)
          enddo
          do ij = ADM_GMIN_PL+1, ADM_GMAX_PL
             do d = GRD_XDIR, GRD_ZDIR
                v_pl(d,1,ij) = GRD_xt_pl(ij-1,K0,l,d)
                v_pl(d,2,ij) = GRD_xt_pl(ij,  K0,l,d)
             enddo
          enddo

          do ij = ADM_GMIN_PL, ADM_GMAX_PL
             call MISC_mk_gmtrvec( v_pl(:,1,ij), v_pl(:,2,ij), tvec(:), nvec(:), &
                                   GMTR_polygon_type, GRD_rscale                 )

             GMTR_A_var_pl(ij,K0,l,GMTR_A_HTX:GMTR_A_HTZ) = tvec(1:3)
             GMTR_A_var_pl(ij,K0,l,GMTR_A_HNX:GMTR_A_HNZ) = nvec(1:3)
          enddo

       enddo
    endif

    return
  end subroutine xcalc_gmtr_a

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine mk_gmtrvec_on_plane
  !>
  subroutine mk_gmtrvec_on_plane( vFrom, vTo, vT, vN )
    implicit none

    real(8), intent(in)  :: vFrom(3), vTo(3)
    real(8), intent(out) :: vT(3),    vN(3)
    !---------------------------------------------------------------------------

    vT(:) = vTo(:) - vFrom(:)

    vN(1) = -vT(2)
    vN(2) =  vT(1)
    vN(3) =   0.D0

    return
  end subroutine mk_gmtrvec_on_plane

  !-----------------------------------------------------------------------------
  !>
  !> Description of the function %NAME
  !> @return area
  !>
  function get_triangle_are_on_plane( a, b, c ) result(area)
    implicit none

    real(8), intent(in) :: a(3), b(3), c(3)
    real(8)             :: area
    !
    real(8) :: a2b(3), a2c(3)
    real(8) :: len_a2b, len_a2c
    real(8) :: prd
    !---------------------------------------------------------------------------

    a2b(:)  = b(:) - a(:)
    a2c(:)  = c(:) - a(:)
    len_a2b = a2b(1)*a2b(1) + a2b(3)*a2b(3) + a2b(3)*a2b(3) ! |a->b|**2
    len_a2c = a2c(1)*a2c(1) + a2c(3)*a2c(3) + a2c(3)*a2c(3) ! |a->c|**2
    prd     = a2b(1)*a2c(1) + a2b(2)*a2c(2) + a2b(3)*a2c(3) ! (a->b)*(a->c)

    area    = 0.5D0 * sqrt( len_a2b * len_a2c - prd*prd )

  end function get_triangle_are_on_plane

end module mod_gmtr
!-------------------------------------------------------------------------------
