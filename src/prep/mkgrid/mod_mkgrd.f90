!-------------------------------------------------------------------------------
!>
!! Module mkgrd
!!
!! @par Description
!!          Making grid systems based on the icosahedral grid configuration
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2013-05-1  (H.Yashiro) NICAM-DC
!!
!<
module mod_mkgrd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS,    &
     ADM_MAXFNAME
  use mod_grd, only: &
     GRD_XDIR, &
     GRD_YDIR, &
     GRD_ZDIR, &
     GRD_x,    &
     GRD_x_pl, &
     GRD_xt,   &
     GRD_xt_pl
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKGRD_setup
  public :: MKGRD_standard
  public :: MKGRD_spring
  public :: MKGRD_prerotate
  public :: MKGRD_stretch
  public :: MKGRD_shrink
  public :: MKGRD_rotate
  public :: MKGRD_gravcenter
  public :: MKGRD_diagnosis

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=ADM_MAXFNAME), public :: MKGRD_IN_BASENAME  = ''
  character(len=ADM_MAXFNAME), public :: MKGRD_OUT_BASENAME = ''
  character(len=ADM_NSYS),     public :: MKGRD_IN_io_mode   = 'ADVANCED'
  character(len=ADM_NSYS),     public :: MKGRD_OUT_io_mode  = 'ADVANCED'

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: dir_vindex = 3

  integer, private, parameter :: I_Xaxis = 1
  integer, private, parameter :: I_Yaxis = 2
  integer, private, parameter :: I_Zaxis = 3

  logical, private :: MKGRD_DOSPRING    = .true.
  logical, private :: MKGRD_DOPREROTATE = .false.
  logical, private :: MKGRD_DOSTRETCH   = .false.
  logical, private :: MKGRD_DOSHRINK    = .false.
  logical, private :: MKGRD_DOROTATE    = .false.

  real(8), private :: MKGRD_spring_beta      = 1.15D0 ! parameter beta for spring dynamics 
  real(8), private :: MKGRD_prerotation_tilt =   0.D0 ! [deg]
  real(8), private :: MKGRD_stretch_alpha    = 1.00D0 ! parameter alpha for stretch
  integer, private :: MKGRD_shrink_level     =      0 ! shrink level (only for 1-diamond experiment)
  real(8), private :: MKGRD_rotation_lon     =   0.D0 ! [deg]
  real(8), private :: MKGRD_rotation_lat     =  90.D0 ! [deg]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKGRD_setup
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_KNONE,     &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_TI,        &
       ADM_TJ
    implicit none

    namelist / PARAM_MKGRD / &
      MKGRD_DOSPRING,         &
      MKGRD_DOPREROTATE,      &
      MKGRD_DOSTRETCH,        &
      MKGRD_DOSHRINK,         &
      MKGRD_DOROTATE,         &
      MKGRD_IN_BASENAME,      &
      MKGRD_IN_io_mode,       &
      MKGRD_OUT_BASENAME,     &
      MKGRD_OUT_io_mode,      &
      MKGRD_spring_beta,      &
      MKGRD_prerotation_tilt, &
      MKGRD_stretch_alpha,    &
      MKGRD_shrink_level,     &
      MKGRD_rotation_lon,     &
      MKGRD_rotation_lat

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Program[mkgrd]/Category[prep]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=PARAM_MKGRD,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** PARAM_MKGRD is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist PARAM_MKGRD. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist PARAM_MKGRD. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=PARAM_MKGRD)

    allocate( GRD_x    (ADM_gall,   ADM_KNONE,ADM_lall,   dir_vindex) )
    allocate( GRD_x_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,dir_vindex) )

    allocate( GRD_xt   (ADM_gall,   ADM_KNONE,ADM_lall,   ADM_TI:ADM_TJ,dir_vindex) )
    allocate( GRD_xt_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              dir_vindex) )

    return
  end subroutine MKGRD_setup

  !-----------------------------------------------------------------------------
  !> Make standard grid system
  subroutine MKGRD_standard
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_rlevel,  &
       ADM_glevel,  &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_gmax,    &
       ADM_gmin,    &
       ADM_gslf_pl, &
       ADM_NPL,     &
       ADM_SPL
    use mod_cnst, only: &
       PI => CNST_PI
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(8), allocatable :: r0(:,:,:)
    real(8), allocatable :: r1(:,:,:)
    real(8), allocatable :: g0(:,:,:)
    real(8), allocatable :: g1(:,:,:)

    real(8) :: alpha2, phi

    integer :: rgnid, dmd
    real(8) :: rdmd

    integer :: rgn_all_1d, rgn_all
    integer :: rgnid_dmd, ir, jr

    integer :: nmax, nmax_prev, rl, gl
    integer :: i, j, ij, k, l
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** Make standard grid system'

    k = ADM_KNONE

    alpha2 = 2.D0 * PI / 5.D0
    phi    = asin( cos(alpha2) / (1.D0-cos(alpha2) ) )

    rgn_all_1d = 2**ADM_rlevel
    rgn_all    = rgn_all_1d * rgn_all_1d

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       nmax = 2
       allocate( r0(nmax,nmax,3) )
       allocate( r1(nmax,nmax,3) )

       dmd = (rgnid-1) / rgn_all + 1

       if ( dmd <= 5 ) then ! northern hemisphere
          rdmd = real(dmd-1,kind=8)

          r0(1,1,GRD_XDIR) = cos( phi) * cos(alpha2*rdmd)
          r0(1,1,GRD_YDIR) = cos( phi) * sin(alpha2*rdmd)
          r0(1,1,GRD_ZDIR) = sin( phi)

          r0(2,1,GRD_XDIR) = cos(-phi) * cos(alpha2*(rdmd+0.5D0))
          r0(2,1,GRD_YDIR) = cos(-phi) * sin(alpha2*(rdmd+0.5D0))
          r0(2,1,GRD_ZDIR) = sin(-phi)

          r0(1,2,GRD_XDIR) =  0.D0
          r0(1,2,GRD_YDIR) =  0.D0
          r0(1,2,GRD_ZDIR) =  1.D0

          r0(2,2,GRD_XDIR) = cos( phi) * cos(alpha2*(rdmd+1.0D0))
          r0(2,2,GRD_YDIR) = cos( phi) * sin(alpha2*(rdmd+1.0D0))
          r0(2,2,GRD_ZDIR) = sin( phi)
       else ! southern hemisphere
          rdmd = real(dmd-6,kind=8)

          r0(1,1,GRD_XDIR) = cos(-phi) * cos(-alpha2*(rdmd+0.5D0))
          r0(1,1,GRD_YDIR) = cos(-phi) * sin(-alpha2*(rdmd+0.5D0))
          r0(1,1,GRD_ZDIR) = sin(-phi)

          r0(2,1,GRD_XDIR) =  0.D0
          r0(2,1,GRD_YDIR) =  0.D0
          r0(2,1,GRD_ZDIR) = -1.D0

          r0(1,2,GRD_XDIR) = cos( phi) * cos(-alpha2*rdmd)
          r0(1,2,GRD_YDIR) = cos( phi) * sin(-alpha2*rdmd)
          r0(1,2,GRD_ZDIR) = sin( phi)

          r0(2,2,GRD_XDIR) = cos(-phi) * cos(-alpha2*(rdmd-0.5D0))
          r0(2,2,GRD_YDIR) = cos(-phi) * sin(-alpha2*(rdmd-0.5D0))
          r0(2,2,GRD_ZDIR) = sin(-phi)
       endif

       do rl = 1, ADM_rlevel
          nmax_prev = nmax
          nmax = 2 * (nmax-1) + 1

          deallocate( r1 )
          allocate( r1(nmax,nmax,3) )

          call decomposition( nmax_prev, & ! [IN]
                              r0(:,:,:), & ! [IN]
                              nmax,      & ! [IN]
                              r1(:,:,:)  ) ! [OUT]

          deallocate( r0 )
          allocate( r0(nmax,nmax,3) )

          r0(:,:,:) = r1(:,:,:)
       enddo

       nmax = 2
       allocate( g0(nmax,nmax,3) )
       allocate( g1(nmax,nmax,3) )

       rgnid_dmd = mod(rgnid-1,rgn_all) + 1
       ir        = mod(rgnid_dmd-1,rgn_all_1d) + 1
       jr        = (rgnid_dmd-ir) / rgn_all_1d + 1

       g0(1,1,:) = r0(ir  ,jr  ,:)
       g0(2,1,:) = r0(ir+1,jr  ,:)
       g0(1,2,:) = r0(ir  ,jr+1,:)
       g0(2,2,:) = r0(ir+1,jr+1,:)

       do gl = ADM_rlevel+1, ADM_glevel
          nmax_prev = nmax
          nmax = 2 * (nmax-1) + 1

          deallocate( g1 )
          allocate( g1(nmax,nmax,3) )

          call decomposition( nmax_prev, & ! [IN]
                              g0(:,:,:), & ! [IN]
                              nmax,      & ! [IN]
                              g1(:,:,:)  ) ! [OUT]

          deallocate( g0 )
          allocate( g0(nmax,nmax,3) )

          g0(:,:,:) = g1(:,:,:)
       enddo

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij = suf(i,j)

          GRD_x(ij,k,l,:) = g0(i-1,j-1,:)
       enddo
       enddo

       deallocate( r0 )
       deallocate( r1 )
       deallocate( g0 )
       deallocate( g1 )
    enddo

    ij = ADM_gslf_pl

    GRD_x_pl(ij,k,ADM_NPL,GRD_XDIR) =  0.D0
    GRD_x_pl(ij,k,ADM_NPL,GRD_YDIR) =  0.D0
    GRD_x_pl(ij,k,ADM_NPL,GRD_ZDIR) =  1.D0

    GRD_x_pl(ij,k,ADM_SPL,GRD_XDIR) =  0.D0
    GRD_x_pl(ij,k,ADM_SPL,GRD_YDIR) =  0.D0
    GRD_x_pl(ij,k,ADM_SPL,GRD_ZDIR) = -1.D0

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_standard

  !-----------------------------------------------------------------------------
  !> Apply spring dynamics
  subroutine MKGRD_spring
    use mod_adm, only: &
       ADM_prc_tab,     &
       ADM_prc_me,      &
       ADM_rgn_vnum,    &
       ADM_W,           &
       ADM_glevel,      &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_KNONE,       &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GIpJo,       &
       ADM_GIpJp,       &
       ADM_GIoJp,       &
       ADM_GImJo,       &
       ADM_GImJm,       &
       ADM_GIoJm,       &
       ADM_gmin
    use mod_cnst, only: &
       PI => CNST_PI
    use mod_comm, only: &
       COMM_data_transfer
    use mod_gtl, only: &
       GTL_max, &
       GTL_min
    implicit none

    integer, parameter :: var_vindex = 8
    integer, parameter :: I_Rx   = 1
    integer, parameter :: I_Ry   = 2
    integer, parameter :: I_Rz   = 3
    integer, parameter :: I_Wx   = 4
    integer, parameter :: I_Wy   = 5
    integer, parameter :: I_Wz   = 6
    integer, parameter :: I_Fsum = 7
    integer, parameter :: I_Ek   = 8

    real(8) :: var   ( ADM_gall,   ADM_KNONE,ADM_lall,   var_vindex)
    real(8) :: var_pl( ADM_gall_pl,ADM_KNONE,ADM_lall_pl,var_vindex)

    real(8) :: lambda
    real(8) :: dbar

    real(8) :: Px(ADM_gall,0:6)
    real(8) :: Py(ADM_gall,0:6)
    real(8) :: Pz(ADM_gall,0:6)
    real(8) :: Fx(ADM_gall,0:6)
    real(8) :: Fy(ADM_gall,0:6)
    real(8) :: Fz(ADM_gall,0:6)

    real(8) :: fixed_point(3)

    real(8) :: Ax, Ay, Az
    real(8) :: Ex, Ey, Ez
    real(8) :: Fsumx, Fsumy, Fsumz
    real(8) :: Rx, Ry, Rz
    real(8) :: Wx, Wy, Wz
    real(8) :: len, d, E

    integer, parameter :: itelim = 100000
    integer            :: ite
    real(8) :: Fsum_max, Ek_max

    real(8), parameter :: dump_coef = 1.D0  !> friction coefficent in spring dynamics
    real(8), parameter :: dt        = 2.D-2 !> delta t for solution of spring dynamics
    real(8), parameter :: criteria  = 1.D-4 !> criteria of convergence

    integer :: rgnid
    integer :: ij_singular
    integer :: n, ij, k, l, m
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOSPRING ) return

    k = ADM_KNONE
    ij_singular = suf(ADM_gmin,ADM_gmin)

    var   (:,:,:,:) = 0.D0
    var_pl(:,:,:,:) = 0.D0

    lambda = 2.D0*PI / ( 10.D0*2.D0**(ADM_glevel-1) )

    dbar = MKGRD_spring_beta * lambda

    var   (:,:,:,I_Rx:I_Rz) = GRD_x   (:,:,:,GRD_XDIR:GRD_ZDIR)
    var_pl(:,:,:,I_Rx:I_Rz) = GRD_x_pl(:,:,:,GRD_XDIR:GRD_ZDIR)

    write(ADM_LOG_FID,*) '*** Apply grid modification with spring dynamics'
    write(ADM_LOG_FID,*) '*** spring factor beta  = ', MKGRD_spring_beta
    write(ADM_LOG_FID,*) '*** length lambda       = ', lambda
    write(ADM_LOG_FID,*) '*** delta t             = ', dt
    write(ADM_LOG_FID,*) '*** conversion criteria = ', criteria
    write(ADM_LOG_FID,*) '*** dumping coefficient = ', dump_coef
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,'(3(A16))') 'itelation', 'max. Kinetic E', 'max. forcing'

    !--- Solving spring dynamics
    do ite = 1, itelim

       do l = 1, ADM_lall
          rgnid = ADM_prc_tab(l,ADM_prc_me)

          do n = 1, ADM_IooJoo_nmax
             ij = ADM_IooJoo(n,ADM_GIoJo)

             Px(ij,0) = var(ADM_IooJoo(n,ADM_GIoJo),k,l,I_Rx)
             Px(ij,1) = var(ADM_IooJoo(n,ADM_GIpJo),k,l,I_Rx)
             Px(ij,2) = var(ADM_IooJoo(n,ADM_GIpJp),k,l,I_Rx)
             Px(ij,3) = var(ADM_IooJoo(n,ADM_GIoJp),k,l,I_Rx)
             Px(ij,4) = var(ADM_IooJoo(n,ADM_GImJo),k,l,I_Rx)
             Px(ij,5) = var(ADM_IooJoo(n,ADM_GImJm),k,l,I_Rx)
             Px(ij,6) = var(ADM_IooJoo(n,ADM_GIoJm),k,l,I_Rx)

             Py(ij,0) = var(ADM_IooJoo(n,ADM_GIoJo),k,l,I_Ry)
             Py(ij,1) = var(ADM_IooJoo(n,ADM_GIpJo),k,l,I_Ry)
             Py(ij,2) = var(ADM_IooJoo(n,ADM_GIpJp),k,l,I_Ry)
             Py(ij,3) = var(ADM_IooJoo(n,ADM_GIoJp),k,l,I_Ry)
             Py(ij,4) = var(ADM_IooJoo(n,ADM_GImJo),k,l,I_Ry)
             Py(ij,5) = var(ADM_IooJoo(n,ADM_GImJm),k,l,I_Ry)
             Py(ij,6) = var(ADM_IooJoo(n,ADM_GIoJm),k,l,I_Ry)

             Pz(ij,0) = var(ADM_IooJoo(n,ADM_GIoJo),k,l,I_Rz)
             Pz(ij,1) = var(ADM_IooJoo(n,ADM_GIpJo),k,l,I_Rz)
             Pz(ij,2) = var(ADM_IooJoo(n,ADM_GIpJp),k,l,I_Rz)
             Pz(ij,3) = var(ADM_IooJoo(n,ADM_GIoJp),k,l,I_Rz)
             Pz(ij,4) = var(ADM_IooJoo(n,ADM_GImJo),k,l,I_Rz)
             Pz(ij,5) = var(ADM_IooJoo(n,ADM_GImJm),k,l,I_Rz)
             Pz(ij,6) = var(ADM_IooJoo(n,ADM_GIoJm),k,l,I_Rz)
          enddo

          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
             Px(ij_singular,6) = Px(ij_singular,1)
             Py(ij_singular,6) = Py(ij_singular,1)
             Pz(ij_singular,6) = Pz(ij_singular,1)
          endif

          do m = 1, 6
             do n = 1, ADM_IooJoo_nmax
                ij = ADM_IooJoo(n,ADM_GIoJo)

                ! A = P0 X Pm
                Ax = Py(ij,0) * Pz(ij,m) - Pz(ij,0) * Py(ij,m)
                Ay = Pz(ij,0) * Px(ij,m) - Px(ij,0) * Pz(ij,m)
                Az = Px(ij,0) * Py(ij,m) - Py(ij,0) * Px(ij,m)

                ! e0 = ( P0 X Pm ) X P0
                Ex = Ay * Pz(ij,0) - Az * Py(ij,0)
                Ey = Az * Px(ij,0) - Ax * Pz(ij,0)
                Ez = Ax * Py(ij,0) - Ay * Px(ij,0)

                ! normalize
                len = sqrt( Ex*Ex + Ey*Ey + Ez*Ez )

                ! d = P0 * Pm
                d = acos( Px(ij,0) * Px(ij,m) &
                        + Py(ij,0) * Py(ij,m) &
                        + Pz(ij,0) * Pz(ij,m) )

                Fx(ij,m) = (d-dbar) * Ex / len
                Fy(ij,m) = (d-dbar) * Ey / len
                Fz(ij,m) = (d-dbar) * Ez / len
             enddo
          enddo

          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
             Fx(ij_singular,6) = 0.D0
             Fy(ij_singular,6) = 0.D0
             Fz(ij_singular,6) = 0.D0
          endif

          do n = 1, ADM_IooJoo_nmax
             ij = ADM_IooJoo(n,ADM_GIoJo)

             Fsumx = Fx(ij,1)+Fx(ij,2)+Fx(ij,3)+Fx(ij,4)+Fx(ij,5)+Fx(ij,6)
             Fsumy = Fy(ij,1)+Fy(ij,2)+Fy(ij,3)+Fy(ij,4)+Fy(ij,5)+Fy(ij,6)
             Fsumz = Fz(ij,1)+Fz(ij,2)+Fz(ij,3)+Fz(ij,4)+Fz(ij,5)+Fz(ij,6)

             ! check dw0/dt
             var(ij,k,l,I_Fsum) = sqrt( Fsumx*Fsumx + Fsumy*Fsumy + Fsumz*Fsumz ) / lambda

             Fx(ij,0) = Fsumx - dump_coef * var(ij,k,l,I_Wx)
             Fy(ij,0) = Fsumy - dump_coef * var(ij,k,l,I_Wy)
             Fz(ij,0) = Fsumz - dump_coef * var(ij,k,l,I_Wz)
          enddo

          ! save value of fixed point
          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
             fixed_point(I_Rx) = var(ij_singular,k,l,I_Rx)
             fixed_point(I_Ry) = var(ij_singular,k,l,I_Ry)
             fixed_point(I_Rz) = var(ij_singular,k,l,I_Rz)
          endif

          ! update r0
          do n = 1, ADM_IooJoo_nmax
             ij = ADM_IooJoo(n,ADM_GIoJo)

             Rx = var(ij,k,l,I_Rx) + var(ij,k,l,I_Wx) * dt
             Ry = var(ij,k,l,I_Ry) + var(ij,k,l,I_Wy) * dt
             Rz = var(ij,k,l,I_Rz) + var(ij,k,l,I_Wz) * dt

             ! normalize
             len = sqrt( Rx*Rx + Ry*Ry + Rz*Rz )

             var(ij,k,l,I_Rx) = Rx / len
             var(ij,k,l,I_Ry) = Ry / len
             var(ij,k,l,I_Rz) = Rz / len
          enddo

          ! update w0
          do n = 1, ADM_IooJoo_nmax
             ij = ADM_IooJoo(n,ADM_GIoJo)

             Wx = var(ij,k,l,I_Wx) + Fx(ij,0) * dt
             Wy = var(ij,k,l,I_Wy) + Fy(ij,0) * dt
             Wz = var(ij,k,l,I_Wz) + Fz(ij,0) * dt

             ! horizontalize
             E = var(ij,k,l,I_Rx) * Wx &
               + var(ij,k,l,I_Ry) * Wy &
               + var(ij,k,l,I_Rz) * Wz

             var(ij,k,l,I_Wx) = Wx - E * var(ij,k,l,I_Rx)
             var(ij,k,l,I_Wy) = Wy - E * var(ij,k,l,I_Ry)
             var(ij,k,l,I_Wz) = Wz - E * var(ij,k,l,I_Rz)

             ! kinetic energy
             var(ij,k,l,I_Ek) = 0.5D0 * ( var(ij,k,l,I_Wx)*var(ij,k,l,I_Wx) &
                                        + var(ij,k,l,I_Wy)*var(ij,k,l,I_Wy) &
                                        + var(ij,k,l,I_Wz)*var(ij,k,l,I_Wz) )
          enddo

          ! restore value of fixed point
          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
             var(ij_singular,k,l,:)    = 0.D0
             var(ij_singular,k,l,I_Rx) = fixed_point(I_Rx)
             var(ij_singular,k,l,I_Ry) = fixed_point(I_Ry)
             var(ij_singular,k,l,I_Rz) = fixed_point(I_Rz)
          endif

       enddo ! l loop

       call COMM_data_transfer( var(:,:,:,:), var_pl(:,:,:,:) )

       Fsum_max = GTL_max( var(:,:,:,I_Fsum), var_pl(:,:,:,I_Fsum), 1, 1, 1 )
       Ek_max   = GTL_max( var(:,:,:,I_Ek),   var_pl(:,:,:,I_Ek)  , 1, 1, 1 )

       write(ADM_LOG_FID,'(I16,4(E16.8))') ite, Ek_max, Fsum_max

       if( Fsum_max < criteria ) exit

    enddo ! itelation loop

    GRD_x   (:,:,:,GRD_XDIR:GRD_ZDIR) = var   (:,:,:,I_Rx:I_Rz)
    GRD_x_pl(:,:,:,GRD_XDIR:GRD_ZDIR) = var_pl(:,:,:,I_Rx:I_Rz)

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_spring

  !-----------------------------------------------------------------------------
  !> Apply rotation before stretching, for 1-diamond grid system
  subroutine MKGRD_prerotate
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_lall_pl
    use mod_cnst, only: &
       PI => CNST_PI
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(8) :: g(3)
    real(8) :: angle_y, angle_z, angle_tilt
    real(8) :: alpha2

    real(8) :: d2r
    integer :: ij, k, l
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOPREROTATE ) return

    k = ADM_KNONE

    d2r        = PI / 180.D0
    alpha2     = 2.D0 * PI / 5.D0
    angle_z    = alpha2 / 2.D0
    angle_y    = 0.25D0*PI * ( 3.D0 - sqrt(3.D0) )
    angle_tilt = MKGRD_prerotation_tilt * d2r

    write(ADM_LOG_FID,*) '*** Apply pre-rotation'
    write(ADM_LOG_FID,*) '*** Diamond tilting factor = ', MKGRD_prerotation_tilt
    write(ADM_LOG_FID,*) '*** angle_z   (deg) = ', angle_z    / d2r
    write(ADM_LOG_FID,*) '*** angle_y   (deg) = ', angle_y    / d2r
    write(ADM_LOG_FID,*) '*** angle_tilt(deg) = ', angle_tilt / d2r

    do l = 1, ADM_lall
       do ij = 1, ADM_gall
          g(:) = GRD_x(ij,k,l,:)

          ! align lowermost vertex of diamond to x-z coordinate plane
          call MISC_3dvec_rotation( g(:),    & ! [INOUT]
                                    angle_z, & ! [IN]
                                    I_Zaxis  ) ! [IN]
          ! rotate around y-axis, for fitting the center of diamond to north pole
          call MISC_3dvec_rotation( g(:),    & ! [INOUT]
                                    angle_y, & ! [IN]
                                    I_Yaxis  ) ! [IN]
          ! rotate the diamond around z-axis
          call MISC_3dvec_rotation( g(:),       & ! [INOUT]
                                    angle_tilt, & ! [IN]
                                    I_Zaxis     ) ! [IN]

          GRD_x(ij,k,l,:) = g(:)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l  = 1, ADM_lall_pl
       do ij = 1, ADM_gall_pl
          g(:) = GRD_x_pl(ij,k,l,:)

          ! align lowermost vertex of diamond to x-z coordinate plane
          call MISC_3dvec_rotation( g(:),    & ! [INOUT]
                                    angle_z, & ! [IN]
                                    I_Zaxis  ) ! [IN]
          ! rotate around y-axis, for fitting the center of diamond to north pole
          call MISC_3dvec_rotation( g(:),    & ! [INOUT]
                                    angle_y, & ! [IN]
                                    I_Yaxis  ) ! [IN]
          ! rotate the diamond around z-axis
          call MISC_3dvec_rotation( g(:),       & ! [INOUT]
                                    angle_tilt, & ! [IN]
                                    I_Zaxis     ) ! [IN]

          GRD_x_pl(ij,k,l,:) = g(:)
       enddo
       enddo
    endif

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_prerotate

  !-----------------------------------------------------------------------------
  !> Apply stretching to grid system
  subroutine MKGRD_stretch
    use mod_misc, only: &
       MISC_get_latlon, &
       MISC_get_cartesian
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_lall_pl
    use mod_cnst, only: &
       PI => CNST_PI
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(8) :: lat, lon, lat_trans

    real(8), parameter :: criteria = 1.D-10

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOSTRETCH ) return

    write(ADM_LOG_FID,*) '*** Apply stretch'
    write(ADM_LOG_FID,*) '*** Stretch factor = ', MKGRD_stretch_alpha

    k = ADM_KNONE

    do l = 1, ADM_lall
       do ij = 1, ADM_gall

          call MISC_get_latlon( lat,                    & ! [OUT]
                                lon,                    & ! [OUT]
                                GRD_x(ij,k,l,GRD_XDIR), & ! [IN]
                                GRD_x(ij,k,l,GRD_YDIR), & ! [IN]
                                GRD_x(ij,k,l,GRD_ZDIR)  ) ! [IN]

          if ( 0.5D0*PI-abs(lat) > criteria ) then
             lat_trans = asin( ( MKGRD_stretch_alpha*(1.D0+sin(lat)) / (1.D0-sin(lat)) - 1.D0 ) &
                             / ( MKGRD_stretch_alpha*(1.D0+sin(lat)) / (1.D0-sin(lat)) + 1.D0 ) )
          else
             lat_trans = lat
          endif

          call MISC_get_cartesian( GRD_x(ij,k,l,GRD_XDIR), & ! [OUT]
                                   GRD_x(ij,k,l,GRD_YDIR), & ! [OUT]
                                   GRD_x(ij,k,l,GRD_ZDIR), & ! [OUT]
                                   lat_trans,              & ! [IN]
                                   lon,                    & ! [IN]
                                   1.D0                    ) ! [IN]
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l  = 1, ADM_lall_pl
       do ij = 1, ADM_gall_pl

          call MISC_get_latlon( lat,                       & ! [OUT]
                                lon,                       & ! [OUT]
                                GRD_x_pl(ij,k,l,GRD_XDIR), & ! [IN]
                                GRD_x_pl(ij,k,l,GRD_YDIR), & ! [IN]
                                GRD_x_pl(ij,k,l,GRD_ZDIR)  ) ! [IN]

          if ( 0.5D0*PI-abs(lat) > criteria ) then
             lat_trans = asin( ( MKGRD_stretch_alpha*(1.D0+sin(lat)) / (1.D0-sin(lat)) - 1.D0 ) &
                             / ( MKGRD_stretch_alpha*(1.D0+sin(lat)) / (1.D0-sin(lat)) + 1.D0 ) )
          else
             lat_trans = lat
          endif

          call MISC_get_cartesian( GRD_x_pl(ij,k,l,GRD_XDIR), & ! [OUT]
                                   GRD_x_pl(ij,k,l,GRD_YDIR), & ! [OUT]
                                   GRD_x_pl(ij,k,l,GRD_ZDIR), & ! [OUT]
                                   lat_trans,                 & ! [IN]
                                   lon,                       & ! [IN]
                                   1.D0                       ) ! [IN]
       enddo
       enddo
    endif

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_stretch

  !-----------------------------------------------------------------------------
  !> Apply shrinkng to grid system
  subroutine MKGRD_shrink
    use mod_misc, only: &
       MISC_get_latlon, &
       MISC_get_cartesian
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_lall_pl
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(8) :: o(3), g(3), len

    integer :: ij, k, l, ite
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOSHRINK ) return

    write(ADM_LOG_FID,*) '*** Apply shrink'
    write(ADM_LOG_FID,*) '*** Shrink level = ', MKGRD_shrink_level

    k = ADM_KNONE

    o(GRD_XDIR) = 0.D0
    o(GRD_YDIR) = 0.D0

    do ite = 1, MKGRD_shrink_level
       do l  = 1, ADM_lall
       do ij = 1, ADM_gall
          o(GRD_ZDIR) = sign(1.D0,GRD_x(ij,k,l,GRD_ZDIR))

          g(GRD_XDIR) = GRD_x(ij,k,l,GRD_XDIR) + o(GRD_XDIR)
          g(GRD_YDIR) = GRD_x(ij,k,l,GRD_YDIR) + o(GRD_YDIR)
          g(GRD_ZDIR) = GRD_x(ij,k,l,GRD_ZDIR) + o(GRD_ZDIR)

          len = ( g(GRD_XDIR)*g(GRD_XDIR) &
                + g(GRD_YDIR)*g(GRD_YDIR) &
                + g(GRD_ZDIR)*g(GRD_ZDIR) )

          GRD_x(ij,k,l,GRD_XDIR) = g(GRD_XDIR) / len
          GRD_x(ij,k,l,GRD_YDIR) = g(GRD_YDIR) / len
          GRD_x(ij,k,l,GRD_ZDIR) = g(GRD_ZDIR) / len
       enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
    do ite = 1, MKGRD_shrink_level-1
       do l  = 1, ADM_lall_pl
       do ij = 1, ADM_gall_pl
          o(GRD_ZDIR) = sign(1.D0,GRD_x_pl(ij,k,l,GRD_ZDIR))

          g(GRD_XDIR) = GRD_x_pl(ij,k,l,GRD_XDIR) + o(GRD_XDIR)
          g(GRD_YDIR) = GRD_x_pl(ij,k,l,GRD_YDIR) + o(GRD_YDIR)
          g(GRD_ZDIR) = GRD_x_pl(ij,k,l,GRD_ZDIR) + o(GRD_ZDIR)

          len = ( g(GRD_XDIR)*g(GRD_XDIR) &
                + g(GRD_YDIR)*g(GRD_YDIR) &
                + g(GRD_ZDIR)*g(GRD_ZDIR) )

          GRD_x_pl(ij,k,l,GRD_XDIR) = g(GRD_XDIR) / len
          GRD_x_pl(ij,k,l,GRD_YDIR) = g(GRD_YDIR) / len
          GRD_x_pl(ij,k,l,GRD_ZDIR) = g(GRD_ZDIR) / len
       enddo
       enddo
    enddo
    endif

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_shrink

  !-----------------------------------------------------------------------------
  !> Apply rotation to grid system
  subroutine MKGRD_rotate
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_lall_pl
    use mod_cnst, only: &
       PI => CNST_PI
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(8) :: g(3)
    real(8) :: angle_y, angle_z

    real(8) :: d2r
    integer :: ij, k, l
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOROTATE ) return

    write(ADM_LOG_FID,*) '*** Apply rotation'
    write(ADM_LOG_FID,*) '*** North pole -> Longitude(deg) = ', MKGRD_rotation_lon
    write(ADM_LOG_FID,*) '*** North pole -> Latitude (deg) = ', MKGRD_rotation_lat

    k = ADM_KNONE

    d2r = PI / 180.D0
    angle_y = ( MKGRD_rotation_lat - 90.D0 ) * d2r
    angle_z = - MKGRD_rotation_lon * d2r

    do l = 1, ADM_lall
       do ij = 1, ADM_gall
          g(:) = GRD_x(ij,k,l,:)

          ! rotate around y-axis
          call MISC_3dvec_rotation( g(:),    & ! [INOUT]
                                    angle_y, & ! [IN]
                                    I_Yaxis  ) ! [IN]
          ! rotate around z-axis
          call MISC_3dvec_rotation( g(:),    & ! [INOUT]
                                    angle_z, & ! [IN]
                                    I_Zaxis  ) ! [IN]

          GRD_x(ij,k,l,:) = g(:)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l  = 1, ADM_lall_pl
       do ij = 1, ADM_gall_pl
          g(:) = GRD_x_pl(ij,k,l,:)

          ! rotate around y-axis
          call MISC_3dvec_rotation( g(:),    & ! [INOUT]
                                    angle_y, & ! [IN]
                                    I_Yaxis  ) ! [IN]
          ! rotate around z-axis
          call MISC_3dvec_rotation( g(:),    & ! [INOUT]
                                    angle_z, & ! [IN]
                                    I_Zaxis  ) ! [IN]

          GRD_x_pl(ij,k,l,:) = g(:)
       enddo
       enddo
    endif

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_rotate

  !-----------------------------------------------------------------------------
  !> Arrange gravitational center
  subroutine MKGRD_gravcenter
    use mod_comm, only: &
       COMM_data_transfer
    implicit none
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** Calc gravitational center'

    write(ADM_LOG_FID,*) '*** center -> vertex'
    call MKGRD_center2vertex

    write(ADM_LOG_FID,*) '*** vertex -> center'
    call MKGRD_vertex2center

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_gravcenter

  !-----------------------------------------------------------------------------
  !> Diagnose grid property
  subroutine MKGRD_diagnosis
    use mod_misc, only: &
       MISC_3dvec_cross, &
       MISC_3dvec_dot,   &
       MISC_3dvec_abs
    use mod_adm, only: &
       ADM_prc_tab,  &
       ADM_prc_me,   &
       ADM_rgn_vnum, &
       ADM_W,        &
       ADM_glevel,   &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_KNONE,    &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_TI,       &
       ADM_TJ,       &
       ADM_gmax,     &
       ADM_gmin
    use mod_cnst, only: &
       PI     => CNST_PI,     &
       RADIUS => CNST_ERADIUS
    use mod_gmtr, only: &
       GMTR_P_AREA, &
       GMTR_P_var,  &
       GMTR_P_var_pl
    use mod_gtl, only: &
       GTL_global_sum_srf, &
       GTL_max,            &
       GTL_min
    implicit none

    real(8) :: angle    (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(8) :: angle_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(8) :: length   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(8) :: length_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(8) :: sqarea   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(8) :: sqarea_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(8) :: dummy    (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(8) :: dummy_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)

    real(8) :: len(6), ang(6)
    real(8) :: p(dir_vindex,0:7)
    real(8) :: nvlenC, nvlenS, nv(3)

    real(8) :: nlen, len_tot
    real(8) :: l_mean, area, temp
    real(8) :: sqarea_avg, sqarea_max, sqarea_min
    real(8) :: angle_max,  length_max, length_avg

    real(8) :: global_area
    integer :: global_grid

    integer :: rgnid
    integer :: i, j, ij, k, l, m
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** Diagnose grid property'

    k = ADM_KNONE

    angle    (:,:,:) = 0.D0
    angle_pl (:,:,:) = 0.D0
    length   (:,:,:) = 0.D0
    length_pl(:,:,:) = 0.D0

    nlen    = 0.D0
    len_tot = 0.D0

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij = suf(i,j)

          if (       ADM_rgn_vnum(ADM_W,rgnid) == 3 &
               .AND. i == ADM_gmin                  &
               .AND. j == ADM_gmin                  ) then ! Pentagon

             p(:,0) = GRD_xt(suf(i,  j-1),k,l,ADM_TJ,:)
             p(:,1) = GRD_xt(suf(i,  j  ),k,l,ADM_TI,:)
             p(:,2) = GRD_xt(suf(i,  j  ),k,l,ADM_TJ,:)
             p(:,3) = GRD_xt(suf(i-1,j  ),k,l,ADM_TI,:)
             p(:,4) = GRD_xt(suf(i-1,j-1),k,l,ADM_TJ,:)
             p(:,5) = GRD_xt(suf(i,  j-1),k,l,ADM_TJ,:)
             p(:,6) = GRD_xt(suf(i,  j  ),k,l,ADM_TI,:)

             len(:) = 0.D0
             ang(:) = 0.D0
             do m = 1, 5
                ! vector length of Pm->Pm-1, Pm->Pm+1 
                call MISC_3dvec_dot( len(m), p(:,m), p(:,m-1), p(:,m), p(:,m-1) )
                len(m) = sqrt( len(m) )

                ! angle of Pm-1->Pm->Pm+1
                call MISC_3dvec_dot( nvlenC, p(:,m), p(:,m-1), p(:,m), p(:,m+1) )
                call MISC_3dvec_cross( nv(:), p(:,m), p(:,m-1), p(:,m), p(:,m+1) )
                call MISC_3dvec_abs( nvlenS, nv(:) )

                ang(m) = atan2( nvlenS, nvlenC )
             enddo

             ! maximum/minimum ratio of angle between the cell vertexes
             angle(ij,k,l) = maxval( ang(1:5) ) / minval( ang(1:5) ) - 1.D0

             ! l_mean: side length of regular pentagon =sqrt(area/1.7204774005)
             area   = GMTR_P_var(ij,k,l,GMTR_P_AREA)
             l_mean = sqrt( 4.D0 / sqrt( 25.D0 + 10.D0*sqrt(5.D0)) * area )

             temp = 0.D0
             do m = 1, 5
                nlen    = nlen + 1.D0
                len_tot = len_tot + len(m)

                temp = temp + (len(m)-l_mean) * (len(m)-l_mean)
             enddo
             ! distortion of side length from l_mean
             length(ij,k,l) = sqrt( temp/5.D0 ) / l_mean

          else ! Hexagon

             p(:,0) = GRD_xt(suf(i,  j-1),k,l,ADM_TJ,:)
             p(:,1) = GRD_xt(suf(i,  j  ),k,l,ADM_TI,:)
             p(:,2) = GRD_xt(suf(i,  j  ),k,l,ADM_TJ,:)
             p(:,3) = GRD_xt(suf(i-1,j  ),k,l,ADM_TI,:)
             p(:,4) = GRD_xt(suf(i-1,j-1),k,l,ADM_TJ,:)
             p(:,5) = GRD_xt(suf(i-1,j-1),k,l,ADM_TI,:)
             p(:,6) = GRD_xt(suf(i,  j-1),k,l,ADM_TJ,:)
             p(:,7) = GRD_xt(suf(i,  j  ),k,l,ADM_TI,:)

             len(:) = 0.D0
             ang(:) = 0.D0
             do m = 1, 6
                ! vector length of Pm->Pm-1, Pm->Pm+1 
                call MISC_3dvec_dot( len(m), p(:,m), p(:,m-1), p(:,m), p(:,m-1) )
                len(m) = sqrt( len(m) )

                ! angle of Pm-1->Pm->Pm+1
                call MISC_3dvec_dot( nvlenC, p(:,m), p(:,m-1), p(:,m), p(:,m+1) )
                call MISC_3dvec_cross( nv(:), p(:,m), p(:,m-1), p(:,m), p(:,m+1) )
                call MISC_3dvec_abs( nvlenS, nv(:) )

                ang(m) = atan2( nvlenS, nvlenC )
             enddo

             ! maximum/minimum ratio of angle between the cell vertexes
             angle(ij,k,l) = maxval( ang(:) ) / minval( ang(:) ) - 1.D0

             ! l_mean: side length of equilateral triangle
             area   = GMTR_P_var(ij,k,l,GMTR_P_AREA)
             l_mean = sqrt( 4.D0 / sqrt(3.D0) / 6.D0 * area )

             temp = 0.D0
             do m = 1, 6
                nlen = nlen + 1.D0
                len_tot = len_tot + len(m)

                temp = temp + (len(m)-l_mean)*(len(m)-l_mean)
             enddo
             ! distortion of side length from l_mean
             length(ij,k,l) = sqrt( temp/6.D0 ) / l_mean

          endif
       enddo
       enddo

    enddo

    dummy    (:,:,:) = 1.D0
    dummy_pl (:,:,:) = 1.D0
    global_area = GTL_global_sum_srf( dummy(:,:,:), dummy_pl(:,:,:) )
    global_grid = 10*4**ADM_glevel + 2
    sqarea_avg = sqrt( global_area / real(global_grid,kind=8) )

    sqarea   (:,:,:) = sqrt( GMTR_P_var   (:,:,:,GMTR_P_AREA) )
    sqarea_pl(:,:,:) = sqrt( GMTR_P_var_pl(:,:,:,GMTR_P_AREA) )
    sqarea_max = GTL_max ( sqarea(:,:,:), sqarea_pl(:,:,:), 1, 1, 1 )
    sqarea_min = GTL_min ( sqarea(:,:,:), sqarea_pl(:,:,:), 1, 1, 1 )

    length_avg = len_tot / nlen
    length_max = GTL_max( length(:,:,:), length_pl(:,:,:), 1, 1, 1 )
    angle_max  = GTL_max( angle (:,:,:), angle_pl (:,:,:), 1, 1, 1 )

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '------ Diagnosis result ---'
    write(ADM_LOG_FID,*) '--- ideal  global surface area  = ', 4.D0*PI*RADIUS*RADIUS*1.D-6,' [km2]'
    write(ADM_LOG_FID,*) '--- actual global surface area  = ', global_area*1.D-6,' [km2]'
    write(ADM_LOG_FID,*) '--- global total number of grid = ', global_grid
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '--- average grid interval       = ', sqarea_avg * 1.D-3,' [km]'
    write(ADM_LOG_FID,*) '--- max grid interval           = ', sqarea_max * 1.D-3,' [km]'
    write(ADM_LOG_FID,*) '--- min grid interval           = ', sqarea_min * 1.D-3,' [km]'
    write(ADM_LOG_FID,*) '--- ratio max/min grid interval = ', sqarea_max / sqarea_min
    write(ADM_LOG_FID,*) '--- average length of arc(side) = ', length_avg * 1.D-3,' [km]'
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '--- max length distortion       = ', length_max * 1.D-3,' [km]'
    write(ADM_LOG_FID,*) '--- max angle distortion        = ', angle_max*180.D0/PI,' [deg]'

    return
  end subroutine MKGRD_diagnosis

  !---------------------------------------------------------------------------==
  subroutine decomposition( &
      n0, &
      g0, &
      n1, &
      g1  )
    implicit none

    integer, intent(in)  :: n0
    real(8), intent(in)  :: g0(n0,n0,3)
    integer, intent(in)  :: n1
    real(8), intent(out) :: g1(n1,n1,3)

    real(8) :: r
    integer :: i, j, inew, jnew
    !---------------------------------------------------------------------------

    do i = 1, n0
    do j = 1, n0
       inew = 2 * i - 1
       jnew = 2 * j - 1

       g1(inew,jnew,:) = g0(i,j,:)

       if ( i < n0 ) then
          g1(inew+1,jnew  ,:) = g0(i+1,j  ,:) + g0(i,j,:)
       endif
       if ( j < n0 ) then
          g1(inew  ,jnew+1,:) = g0(i  ,j+1,:) + g0(i,j,:)
       endif
       if ( i < n0 .AND. j < n0 ) then
          g1(inew+1,jnew+1,:) = g0(i+1,j+1,:) + g0(i,j,:)
       endif
    enddo
    enddo

    do i = 1, n1
    do j = 1, n1
       r = sqrt( g1(i,j,1)*g1(i,j,1) &
               + g1(i,j,2)*g1(i,j,2) &
               + g1(i,j,3)*g1(i,j,3) )

       g1(i,j,1) = g1(i,j,1) / r
       g1(i,j,2) = g1(i,j,2) / r
       g1(i,j,3) = g1(i,j,3) / r
    enddo
    enddo

    return
  end subroutine decomposition

  !-----------------------------------------------------------------------------
  !> suffix calculation
  !> @return suf
  function suf(i,j) result(suffix)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: suffix
    integer :: i, j
    !---------------------------------------------------------------------------

    suffix = ADM_gall_1d * (j-1) + i

  end function suf

  !-----------------------------------------------------------------------------
  !> Apply rotation matrix
  subroutine MISC_3dvec_rotation( &
      a,     &
      angle, &
      iaxis  )
    implicit none

    real(8), intent(inout) :: a(3)
    real(8), intent(in)    :: angle
    integer, intent(in)    :: iaxis

    real(8) :: m(3,3), b(3)
    !---------------------------------------------------------------------------

    if ( iaxis == I_Xaxis ) then
       m(1,1) =        1.D0
       m(1,2) =        0.D0
       m(1,3) =        0.D0

       m(2,1) =        0.D0
       m(2,2) =  cos(angle)
       m(2,3) =  sin(angle)

       m(3,1) =        0.D0
       m(3,2) = -sin(angle)
       m(3,3) =  cos(angle)
    elseif( iaxis == I_Yaxis ) then
       m(1,1) =  cos(angle)
       m(1,2) =        0.D0
       m(1,3) = -sin(angle)

       m(2,1) =        0.D0
       m(2,2) =        1.D0
       m(2,3) =        0.D0

       m(3,1) =  sin(angle)
       m(3,2) =        0.D0
       m(3,3) =  cos(angle)
    elseif( iaxis == I_Zaxis ) then
       m(1,1) =  cos(angle)
       m(1,2) =  sin(angle)
       m(1,3) =        0.D0

       m(2,1) = -sin(angle)
       m(2,2) =  cos(angle)
       m(2,3) =        0.D0

       m(3,1) =        0.D0
       m(3,2) =        0.D0
       m(3,3) =        1.D0
    else
       return
    endif

    b(1) = m(1,1) * a(1) + m(1,2) * a(2) + m(1,3) * a(3)
    b(2) = m(2,1) * a(1) + m(2,2) * a(2) + m(2,3) * a(3)
    b(3) = m(3,1) * a(1) + m(3,2) * a(2) + m(3,3) * a(3)

    a(:) = b(:)

    return
  end subroutine MISC_3dvec_rotation

  !-----------------------------------------------------------------------------
  !> gnomonic projection
  subroutine MISC_latlon2gnom( &
      x,          &
      y,          &
      lat,        &
      lon,        &
      lat_center, &
      lon_center  )
    implicit none

    real(8), intent(out) :: x          !> gnomonic, x
    real(8), intent(out) :: y          !> gnomonic, y
    real(8), intent(in)  :: lat        !> spheric, latitude
    real(8), intent(in)  :: lon        !> spheric, longitude
    real(8), intent(in)  :: lat_center !> projection center, latitude
    real(8), intent(in)  :: lon_center !> projection center, longitude

    real(8) :: cosc
    !---------------------------------------------------------------------------

    cosc = sin(lat_center) * sin(lat) &
         + cos(lat_center) * cos(lat) * cos(lon-lon_center)

    x = ( cos(lat) * sin(lon-lon_center) ) / cosc
    y = ( cos(lat_center) * sin(lat)                       &
        - sin(lat_center) * cos(lat) * cos(lon-lon_center) ) / cosc

  end subroutine MISC_latlon2gnom

  !-----------------------------------------------------------------------------
  !> gnomonic projection (inverse)
  subroutine MISC_gnom2latlon( &
      lat,        &
      lon,        &
      x,          &
      y,          &
      lat_center, &
      lon_center  )
    implicit none

    real(8), intent(out) :: lat        !> spheric, latitude
    real(8), intent(out) :: lon        !> spheric, longitude
    real(8), intent(in)  :: x          !> gnomonic, x
    real(8), intent(in)  :: y          !> gnomonic, y
    real(8), intent(in)  :: lat_center !> projection center, latitude
    real(8), intent(in)  :: lon_center !> projection center, longitude

    real(8) :: rho, c
    !---------------------------------------------------------------------------

    rho = sqrt( x*x + y*y )

    if ( rho == 0 ) then ! singular point
       lat = lat_center
       lon = lon_center
       return
    endif

    c = atan( rho )

    lon = lon_center + atan2( x*sin(c), ( rho*cos(lat_center)*cos(c) - y*sin(lat_center)*sin(c) ) )
    lat = asin( cos(c)*sin(lat_center) + y*sin(c)*cos(lat_center) / rho )

    return
  end subroutine MISC_gnom2latlon

  !-----------------------------------------------------------------------------
  !> Make center grid -> vertex grid
  subroutine MKGRD_center2vertex
    use mod_misc, only: &
       MISC_3dvec_cross, &
       MISC_3dvec_dot,   &
       MISC_3dvec_abs
    use mod_adm, only : &
      ADM_W,          &
      ADM_TI,         &
      ADM_TJ,         &
      ADM_KNONE,      &
      ADM_prc_tab,    &
      ADM_PRC_PL,     &
      ADM_prc_me,     &
      ADM_rgn_vnum,   &
      ADM_lall,       &
      ADM_gall,       &
      ADM_gmax,       &
      ADM_gmin,       &
      ADM_lall_pl,    &
      ADM_gall_pl,    &
      ADM_GSLF_PL,    &
      ADM_GMAX_PL,    &
      ADM_GMIN_PL,    &
      ADM_ImoJmo_nmax, &
      ADM_ImoJmo,      &
      ADM_GIoJo,       &
      ADM_GIoJp,       &
      ADM_GIpJp,       &
      ADM_GIpJo
    implicit none

    real(8) :: v    (dir_vindex,ADM_gall   ,4,ADM_TI:ADM_TJ)
    real(8) :: v_pl (dir_vindex,ADM_gall_pl,4)
    real(8) :: w    (dir_vindex,ADM_gall   ,3)
    real(8) :: w_pl (dir_vindex,ADM_gall_pl,3)
    real(8) :: gc   (dir_vindex,ADM_gall   )
    real(8) :: gc_pl(dir_vindex,ADM_gall_pl)

    real(8), parameter :: o(3) = 0.D0

    real(8) :: w_lenS, w_lenC, gc_len

    integer :: rgnid
    integer :: oo, po, pp, op
    integer :: k, l, m, n, t
    !---------------------------------------------------------------------------

    k  = ADM_KNONE

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do n = 1, ADM_ImoJmo_nmax
          oo = ADM_ImoJmo(n,ADM_GIoJo)
          po = ADM_ImoJmo(n,ADM_GIpJo)
          pp = ADM_ImoJmo(n,ADM_GIpJp)
          op = ADM_ImoJmo(n,ADM_GIoJp)

          v(GRD_XDIR,oo,1,ADM_TI) = GRD_x(oo,k,l,GRD_XDIR)
          v(GRD_XDIR,oo,2,ADM_TI) = GRD_x(po,k,l,GRD_XDIR)
          v(GRD_XDIR,oo,3,ADM_TI) = GRD_x(pp,k,l,GRD_XDIR)
          v(GRD_XDIR,oo,4,ADM_TI) = GRD_x(oo,k,l,GRD_XDIR)

          v(GRD_YDIR,oo,1,ADM_TI) = GRD_x(oo,k,l,GRD_YDIR)
          v(GRD_YDIR,oo,2,ADM_TI) = GRD_x(po,k,l,GRD_YDIR)
          v(GRD_YDIR,oo,3,ADM_TI) = GRD_x(pp,k,l,GRD_YDIR)
          v(GRD_YDIR,oo,4,ADM_TI) = GRD_x(oo,k,l,GRD_YDIR)

          v(GRD_ZDIR,oo,1,ADM_TI) = GRD_x(oo,k,l,GRD_ZDIR)
          v(GRD_ZDIR,oo,2,ADM_TI) = GRD_x(po,k,l,GRD_ZDIR)
          v(GRD_ZDIR,oo,3,ADM_TI) = GRD_x(pp,k,l,GRD_ZDIR)
          v(GRD_ZDIR,oo,4,ADM_TI) = GRD_x(oo,k,l,GRD_ZDIR)

          v(GRD_XDIR,oo,1,ADM_TJ) = GRD_x(oo,k,l,GRD_XDIR)
          v(GRD_XDIR,oo,2,ADM_TJ) = GRD_x(pp,k,l,GRD_XDIR)
          v(GRD_XDIR,oo,3,ADM_TJ) = GRD_x(op,k,l,GRD_XDIR)
          v(GRD_XDIR,oo,4,ADM_TJ) = GRD_x(oo,k,l,GRD_XDIR)

          v(GRD_YDIR,oo,1,ADM_TJ) = GRD_x(oo,k,l,GRD_YDIR)
          v(GRD_YDIR,oo,2,ADM_TJ) = GRD_x(pp,k,l,GRD_YDIR)
          v(GRD_YDIR,oo,3,ADM_TJ) = GRD_x(op,k,l,GRD_YDIR)
          v(GRD_YDIR,oo,4,ADM_TJ) = GRD_x(oo,k,l,GRD_YDIR)

          v(GRD_ZDIR,oo,1,ADM_TJ) = GRD_x(oo,k,l,GRD_ZDIR)
          v(GRD_ZDIR,oo,2,ADM_TJ) = GRD_x(pp,k,l,GRD_ZDIR)
          v(GRD_ZDIR,oo,3,ADM_TJ) = GRD_x(op,k,l,GRD_ZDIR)
          v(GRD_ZDIR,oo,4,ADM_TJ) = GRD_x(oo,k,l,GRD_ZDIR)
       enddo

       !--- execetion for the north and south.
       v(:,suf(ADM_gmax,ADM_gmin-1),1,ADM_TI) = v(:,suf(ADM_gmax,ADM_gmin-1),1,ADM_TJ)
       v(:,suf(ADM_gmax,ADM_gmin-1),2,ADM_TI) = v(:,suf(ADM_gmax,ADM_gmin-1),2,ADM_TJ)
       v(:,suf(ADM_gmax,ADM_gmin-1),3,ADM_TI) = v(:,suf(ADM_gmax,ADM_gmin-1),3,ADM_TJ)
       v(:,suf(ADM_gmax,ADM_gmin-1),4,ADM_TI) = v(:,suf(ADM_gmax,ADM_gmin-1),4,ADM_TJ)

       v(:,suf(ADM_gmin-1,ADM_gmax),1,ADM_TJ) = v(:,suf(ADM_gmin-1,ADM_gmax),1,ADM_TI)
       v(:,suf(ADM_gmin-1,ADM_gmax),2,ADM_TJ) = v(:,suf(ADM_gmin-1,ADM_gmax),2,ADM_TI)
       v(:,suf(ADM_gmin-1,ADM_gmax),3,ADM_TJ) = v(:,suf(ADM_gmin-1,ADM_gmax),3,ADM_TI)
       v(:,suf(ADM_gmin-1,ADM_gmax),4,ADM_TJ) = v(:,suf(ADM_gmin-1,ADM_gmax),4,ADM_TI)

       !--- exception for the west
       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          oo = suf(ADM_gmin-1,ADM_gmin-1)
          po = suf(ADM_gmin  ,ADM_gmin-1)

          v(:,oo,1,ADM_TI) = v(:,po,1,ADM_TJ)
          v(:,oo,2,ADM_TI) = v(:,po,2,ADM_TJ)
          v(:,oo,3,ADM_TI) = v(:,po,3,ADM_TJ)
          v(:,oo,4,ADM_TI) = v(:,po,4,ADM_TJ)
       endif

       do t = ADM_TI, ADM_TJ
          do m = 1, 3
             do n = 1, ADM_ImoJmo_nmax
                oo = ADM_ImoJmo(n,ADM_GIoJo)

                call MISC_3dvec_dot  ( w_lenC,    o(:), v(:,oo,m,t), o(:), v(:,oo,m+1,t) )
                call MISC_3dvec_cross( w(:,oo,m), o(:), v(:,oo,m,t), o(:), v(:,oo,m+1,t) )
                call MISC_3dvec_abs  ( w_lenS, w(:,oo,m) )

                w(:,oo,m) = w(:,oo,m) / w_lenS * atan2( w_lenS, w_lenC )
             enddo
          enddo

          do n = 1, ADM_ImoJmo_nmax
             oo = ADM_ImoJmo(n,ADM_GIoJo)

             gc(:,oo) = w(:,oo,1) &
                      + w(:,oo,2) &
                      + w(:,oo,3)

             call MISC_3dvec_abs( gc_len, gc(:,oo) )

             GRD_xt(oo,k,l,t,GRD_XDIR) = gc(GRD_XDIR,oo) / gc_len
             GRD_xt(oo,k,l,t,GRD_YDIR) = gc(GRD_YDIR,oo) / gc_len
             GRD_xt(oo,k,l,t,GRD_ZDIR) = gc(GRD_ZDIR,oo) / gc_len
          enddo
       enddo

    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          do oo = ADM_GMIN_PL, ADM_GMAX_PL-1
             v_pl(:,oo,1) = GRD_x_pl(ADM_GSLF_PL,k,l,:)
             v_pl(:,oo,2) = GRD_x_pl(oo+1,       k,l,:)
             v_pl(:,oo,3) = GRD_x_pl(oo  ,       k,l,:)
             v_pl(:,oo,4) = GRD_x_pl(ADM_GSLF_PL,k,l,:)
          enddo
           v_pl(:,ADM_GMAX_PL,1) = GRD_x_pl(ADM_GSLF_PL,k,l,:)
           v_pl(:,ADM_GMAX_PL,2) = GRD_x_pl(ADM_GMIN_PL,k,l,:)
           v_pl(:,ADM_GMAX_PL,3) = GRD_x_pl(ADM_GMAX_PL,k,l,:)
           v_pl(:,ADM_GMAX_PL,4) = GRD_x_pl(ADM_GSLF_PL,k,l,:)

          do n = ADM_GMIN_PL, ADM_GMAX_PL
             do m = 1, 3
                call MISC_3dvec_dot  ( w_lenC,      o(:), v_pl(:,n,m), o(:), v_pl(:,n,m+1) )
                call MISC_3dvec_cross( w_pl(:,n,m), o(:), v_pl(:,n,m), o(:), v_pl(:,n,m+1) )
                call MISC_3dvec_abs  ( w_lenS, w_pl(:,n,m) )

                w_pl(:,n,m) = w_pl(:,n,m) / w_lenS * atan2( w_lenS, w_lenC )
             enddo

             gc_pl(:,n) = w_pl(:,n,1) &
                        + w_pl(:,n,2) &
                        + w_pl(:,n,3)

             call MISC_3dvec_abs( gc_len, gc_pl(:,n) )

             GRD_xt_pl(n,k,l,GRD_XDIR) = gc_pl(GRD_XDIR,n) / gc_len
             GRD_xt_pl(n,k,l,GRD_YDIR) = gc_pl(GRD_YDIR,n) / gc_len
             GRD_xt_pl(n,k,l,GRD_ZDIR) = gc_pl(GRD_ZDIR,n) / gc_len
          enddo
       enddo
    endif

    return
  end subroutine MKGRD_center2vertex

  !-----------------------------------------------------------------------------
  !> Make vertex grid -> center grid
  subroutine MKGRD_vertex2center
    use mod_misc, only: &
       MISC_3dvec_cross, &
       MISC_3dvec_dot,   &
       MISC_3dvec_abs
    use mod_adm, only : &
      ADM_W,          &
      ADM_TI,         &
      ADM_TJ,         &
      ADM_KNONE,      &
      ADM_prc_tab,    &
      ADM_PRC_PL,     &
      ADM_prc_me,     &
      ADM_rgn_vnum,   &
      ADM_lall,       &
      ADM_gall,       &
      ADM_gmin,       &
      ADM_lall_pl,    &
      ADM_gall_pl,    &
      ADM_GSLF_PL,    &
      ADM_GMAX_PL,    &
      ADM_GMIN_PL,    &
      ADM_IooJoo_nmax, &
      ADM_IooJoo,      &
      ADM_GIoJo,       &
      ADM_GIoJm,       &
      ADM_GImJm,       &
      ADM_GImJo
    use mod_grd, only : &
      GRD_XDIR,       &
      GRD_YDIR,       &
      GRD_ZDIR,       &
      GRD_x,          &
      GRD_x_pl,       &
      GRD_xt,         &
      GRD_xt_pl
    implicit none

    real(8) :: v    (dir_vindex,ADM_gall   ,7)
    real(8) :: v_pl (dir_vindex,ADM_gall_pl,6)
    real(8) :: w    (dir_vindex,ADM_gall   ,6)
    real(8) :: w_pl (dir_vindex,ADM_gall_pl,5)
    real(8) :: gc   (dir_vindex,ADM_gall   )
    real(8) :: gc_pl(dir_vindex,ADM_gall_pl)

    real(8), parameter :: o(3) = 0.D0

    real(8) :: w_lenC, w_lenS, gc_len

    integer :: rgnid
    integer :: oo, mo, mm, om
    integer :: k, l, m, n
    !---------------------------------------------------------------------------

    k  = ADM_KNONE

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do n = 1, ADM_IooJoo_nmax
          oo = ADM_IooJoo(n,ADM_GIoJo)
          mo = ADM_IooJoo(n,ADM_GImJo)
          mm = ADM_IooJoo(n,ADM_GImJm)
          om = ADM_IooJoo(n,ADM_GIoJm)

          v(GRD_XDIR,oo,1) = GRD_xt(om,k,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,oo,2) = GRD_xt(oo,k,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,oo,3) = GRD_xt(oo,k,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,oo,4) = GRD_xt(mo,k,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,oo,5) = GRD_xt(mm,k,l,ADM_TJ,GRD_XDIR)
          v(GRD_XDIR,oo,6) = GRD_xt(mm,k,l,ADM_TI,GRD_XDIR)
          v(GRD_XDIR,oo,7) = GRD_xt(om,k,l,ADM_TJ,GRD_XDIR)

          v(GRD_YDIR,oo,1) = GRD_xt(om,k,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,oo,2) = GRD_xt(oo,k,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,oo,3) = GRD_xt(oo,k,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,oo,4) = GRD_xt(mo,k,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,oo,5) = GRD_xt(mm,k,l,ADM_TJ,GRD_YDIR)
          v(GRD_YDIR,oo,6) = GRD_xt(mm,k,l,ADM_TI,GRD_YDIR)
          v(GRD_YDIR,oo,7) = GRD_xt(om,k,l,ADM_TJ,GRD_YDIR)

          v(GRD_ZDIR,oo,1) = GRD_xt(om,k,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,oo,2) = GRD_xt(oo,k,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,oo,3) = GRD_xt(oo,k,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,oo,4) = GRD_xt(mo,k,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,oo,5) = GRD_xt(mm,k,l,ADM_TJ,GRD_ZDIR)
          v(GRD_ZDIR,oo,6) = GRD_xt(mm,k,l,ADM_TI,GRD_ZDIR)
          v(GRD_ZDIR,oo,7) = GRD_xt(om,k,l,ADM_TJ,GRD_ZDIR)
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          oo = suf(ADM_gmin,ADM_gmin)
          v(:,oo,6) = v(:,oo,1)
          v(:,oo,7) = v(:,oo,2)
       endif

       do m = 1, 6
          do n = 1, ADM_IooJoo_nmax
             oo = ADM_IooJoo(n,ADM_GIoJo)

             call MISC_3dvec_dot  ( w_lenC,    o(:), v(:,oo,m), o(:), v(:,oo,m+1) )
             call MISC_3dvec_cross( w(:,oo,m), o(:), v(:,oo,m), o(:), v(:,oo,m+1) )
             call MISC_3dvec_abs  ( w_lenS, w(:,oo,m) )

             w(:,oo,m) = w(:,oo,m) / w_lenS * atan2( w_lenS, w_lenC )
          enddo
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          w(:,suf(ADM_gmin,ADM_gmin),6) = 0.D0
       endif

       do n = 1, ADM_IooJoo_nmax
          oo = ADM_IooJoo(n,ADM_GIoJo)

          gc(:,oo) = w(:,oo,1) &
                   + w(:,oo,2) &
                   + w(:,oo,3) &
                   + w(:,oo,4) &
                   + w(:,oo,5) &
                   + w(:,oo,6)

          call MISC_3dvec_abs( gc_len, gc(:,oo) )

          GRD_x(oo,k,l,:) = gc(:,oo) / gc_len
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1,ADM_lall_pl
          oo = ADM_GSLF_PL

          do n = ADM_GMIN_PL, ADM_GMAX_PL

             v_pl(:,oo,n-ADM_GMIN_PL+1) = GRD_xt_pl(n,k,l,:)

          enddo
          v_pl(:,oo,6) = v_pl(:,oo,1)

          do m = 1, 5
             call MISC_3dvec_dot  ( w_lenC,       o(:), v_pl(:,oo,m), o(:), v_pl(:,oo,m+1) )
             call MISC_3dvec_cross( w_pl(:,oo,m), o(:), v_pl(:,oo,m), o(:), v_pl(:,oo,m+1) )
             call MISC_3dvec_abs  ( w_lenS, w_pl(:,oo,m) )

             w_pl(:,oo,m) = w_pl(:,oo,m) / w_lenS * atan2( w_lenS, w_lenC )
          enddo

          gc_pl(:,oo) = w_pl(:,oo,1) &
                      + w_pl(:,oo,2) &
                      + w_pl(:,oo,3) &
                      + w_pl(:,oo,4) &
                      + w_pl(:,oo,5)

          call MISC_3dvec_abs( gc_len, gc_pl(:,oo) )

          GRD_x_pl(oo,k,l,:) = -gc_pl(:,oo) / gc_len
       enddo
    endif

    return
  end subroutine MKGRD_vertex2center

end module mod_mkgrd
!-------------------------------------------------------------------------------
