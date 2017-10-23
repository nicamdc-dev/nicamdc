!-------------------------------------------------------------------------------
!> Module administration
!!
!! @par Description
!!         This module is for the management of process and region
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_adm
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: ADM_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !====== Basic definition & information ======

  integer,  public, parameter :: I_l   = 1   !< local region
  integer,  public, parameter :: I_prc = 2   !< process

  integer,  public, parameter :: I_RGNID = 1 !< region id
  integer,  public, parameter :: I_DIR   = 2 !< direction

  ! Identifiers of directions of region edges
  integer,  public, parameter :: I_SW = 1    !< south west
  integer,  public, parameter :: I_NW = 2    !< north west
  integer,  public, parameter :: I_NE = 3    !< north east
  integer,  public, parameter :: I_SE = 4    !< south east

  ! Identifiers of directions of region vertices
  integer,  public, parameter :: I_W = 1     !< west
  integer,  public, parameter :: I_N = 2     !< north
  integer,  public, parameter :: I_E = 3     !< east
  integer,  public, parameter :: I_S = 4     !< south

  ! Identifier of poles (north pole or south pole)
  integer,  public, parameter :: I_NPL = 1   !< north pole
  integer,  public, parameter :: I_SPL = 2   !< south pole

  ! Identifier of triangle element (i-axis-side or j-axis side)
  integer,  public, parameter :: ADM_TI = 1
  integer,  public, parameter :: ADM_TJ = 2

  ! Identifier of arc element (i-axis-side, ij-axis side, or j-axis side)
  integer,  public, parameter :: ADM_AI  = 1
  integer,  public, parameter :: ADM_AIJ = 2
  integer,  public, parameter :: ADM_AJ  = 3

  ! Identifier of 1 variable
  integer,  public, parameter :: ADM_KNONE = 1

  ! dimension of the spacial vector
  integer,  public, parameter :: ADM_nxyz = 3

  !#############################################################################
  ! Global Parameters
  !#############################################################################
  integer, public            :: GLOBAL_glevel                ! grid     division level
  integer, public            :: GLOBAL_rlevel                ! region   division level
  integer, public            :: GLOBAL_elevel                ! ensemble division level
  integer, public            :: GLOBAL_COMM_WORLD            ! communication world per member
  integer, public            :: GLOBAL_prc_all               ! number of MPI process
  integer, public            :: GLOBAL_prc_me                ! my process ID
  integer, public, parameter :: GLOBAL_prc_master  = 1       ! master process ID

  logical, public            :: GLOBAL_ensemble_on = .false. ! ensemble run is on?

  character(len=H_LONG), public              :: GLOBAL_prefix_dir    = '' !< prefix,    directory name
  character(len=H_LONG), public              :: GLOBAL_extension_ens = '' !< extension, ensemble name
  character(len=H_LONG), public, allocatable :: GLOBAL_extension_rgn(:)   !< extension, region name

#ifdef _FIXEDINDEX_
  include 'inc_index.h'
#else
  !#############################################################################
  ! Basic Index Parameters
  !#############################################################################

  ! main parameter
  integer,  public            :: ADM_glevel           ! grid   division level
  integer,  public            :: ADM_rlevel           ! region division level
  integer,  public            :: ADM_vlayer           ! number of vertical layer
  integer,  public            :: ADM_DMD              ! number of diamond

  ! region
  integer,  public            :: ADM_rgn_nmax         ! number of regular region
  integer,  public            :: ADM_lall             ! number of regular region per process
  integer,  public, parameter :: ADM_rgn_nmax_pl =  2 ! number of pole    region
  integer,  public, parameter :: ADM_lall_pl     =  2 ! number of pole    region per process

  ! horizontal grid
  integer,  public            :: ADM_gall             ! number of horizontal grid per regular region
  integer,  public            :: ADM_gall_in          ! number of horizontal grid (inner part)
  integer,  public            :: ADM_gall_1d          ! number of horizontal grid (1D)
  integer,  public            :: ADM_gmin             ! start index of 1D horizontal grid
  integer,  public            :: ADM_gmax             ! end   index of 1D horizontal grid

  integer,  public            :: ADM_vlink       = -1 ! maximum number of vertex linkage, ICO:5, PSP:6, LCP, MLCP:k
  integer,  public            :: ADM_gall_pl          ! number of horizontal grid for pole region
  integer,  public, parameter :: ADM_gslf_pl     =  1 ! index for pole point
  integer,  public, parameter :: ADM_gmin_pl     =  2 ! start index of grid around the pole point
  integer,  public            :: ADM_gmax_pl          ! end   index of grid around the pole point

  ! vertical grid
  integer,  public            :: ADM_kall             ! number of vertical grid
  integer,  public            :: ADM_kmin             ! start index of vertical grid
  integer,  public            :: ADM_kmax             ! end   index of vertical grid
#endif

  !
  !====== Information for processes ======
  !

  integer,  public            :: ADM_prc_me               ! my process ID
  integer,  public, parameter :: ADM_prc_master = 1       ! master process ID
  integer,  public            :: ADM_prc_pl               ! process ID which manages the pole regions
  logical,  public            :: ADM_have_pl              ! this ID manages pole region?

  !
  !====== Information for processes-region relationship ======
  !
  integer,  public, parameter   :: RGNMNG_llim = 2560          !< maximum number of region per process.

  integer,  public, allocatable :: RGNMNG_edge_tab   (:,:,:)   !< region link information (for 4 edges)

  integer,  public, allocatable :: RGNMNG_vert_num   (:,:)     !< number of region around the vertex (4 vertexes)
  integer,  public, allocatable :: RGNMNG_vert_tab   (:,:,:,:) !< region link information (for 4 vertexes)
  integer,  public, allocatable :: RGNMNG_vert_tab_pl(:,:,:)   !< region link information (for 4 vertexes)

  integer,  public, allocatable :: RGNMNG_lnum(:)              !< number of region in each process
  integer,  public, allocatable :: RGNMNG_lp2r(:,:)            !< l,prc       => rgn
  integer,  public, allocatable :: RGNMNG_r2lp(:,:)            !< rgn         => l,prc
  integer,  public, allocatable :: RGNMNG_l2r (:)              !< l,prc_me    => rgn

  integer,  public              :: RGNMNG_r2p_pl(I_NPL:I_SPL)  ! process ID which have the pole regions
  integer,  public              :: RGNMNG_rgn4pl(I_NPL:I_SPL)  !< region, having pole data in the halo

  !
  !====== Information for regions ======
  !
  integer,  public              :: ADM_l_me               ! Present Local region number

  logical,  public, allocatable :: ADM_have_sgp(:)        ! region have singlar point?

  !
  !====== Information for grids ======
  !
  character(len=H_SHORT), public :: ADM_HGRID_SYSTEM = 'ICO' ! [XTMS] Horizontal Grid type
                                                     ! 'ICO'      icosahedral
                                                     ! 'ICO-XTMS' icosahedral but XTMS is used in oprt
                                                     ! 'LCP'      Lambert-cornial (including PSP)
                                                     ! 'MLCP'     Mercator+Lambert-cornial
                                                     ! 'MLCP-OLD' OLD vergion (only for s=1)

  integer,                public :: ADM_XTMS_MLCP_S  = 1 ! [XTMS] Number of segment for MLCP

  !
  !====== Information for ensemble ======
  !
  integer, public              :: ENS_COMM_WORLD     ! communication world for ensemble
  integer, public              :: ENS_prc_all        ! number of MPI process in the ensemble
  integer, public              :: ENS_prc_me         ! my process ID in the ensemble
  integer, public, parameter   :: ENS_prc_master = 1 ! master process ID in the emsemble
  logical, public              :: ENS_have_pl        ! this ID manages pole grid?

  integer, public              :: ENS_elevel         ! ensemble division level
  integer, public              :: ENS_div_nmax       ! number of divided part
  integer, public              :: ENS_gall           ! number of horizontal grid
  integer, public, parameter   :: ENS_gall_pl = 2    ! number of horizontal grid, pole point
  integer, public              :: ENS_eall           ! number of ensemble

  !
  !====== Information for processes-division relationship ======
  !
  integer, public              :: ENS_ediv           ! number of divided part per process
  integer, public              :: ENS_gdiv           ! number of horizontal grid per divided part
  integer, public              :: ENS_gdiv_1d        ! number of horizontal grid(1D) per divided part

  integer, public, allocatable :: ENS_d2e(:)         ! division ID -> ensemble member ID
  integer, public, allocatable :: ENS_d2l(:)         ! division ID -> local l
  integer, public, allocatable :: ENS_imin(:)        ! start i-index of 1D horizontal grid for devide
  integer, public, allocatable :: ENS_imax(:)        ! end   i-index of 1D horizontal grid for divide
  integer, public, allocatable :: ENS_jmin(:)        ! start j-index of 1D horizontal grid for devide
  integer, public, allocatable :: ENS_jmax(:)        ! end   j-index of 1D horizontal grid for divide

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: RGNMNG_setup
  private :: RGNMNG_input
  private :: RGNMNG_output
  private :: RGNMNG_generate_ico
  private :: RGNMNG_vertex_walkaround
  private :: output_info

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private            :: RGNMNG_in_fname    = ''  !< input  file name for region management file
  character(len=H_LONG), private            :: RGNMNG_out_fname   = ''  !< output file name for region management file

  character(len=2),      private, parameter :: RGNMNG_edgename(4) = (/'SW','NW','NE','SE'/) !< name table
  character(len=2),      private, parameter :: RGNMNG_vertname(4) = (/'W ','N ','E ','S '/) !< name table

  logical, private :: debug = .false. !< debug option

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ADM_setup
    use mod_process, only: &
       PRC_MPIstop, &
       PRC_myrank,  &
       PRC_nprocs
    implicit none

    integer               :: glevel      = -1 !> grid division level
    integer               :: rlevel      = -1 !> region division level
    integer               :: vlayer      =  1 !> number of inner vertical layer
    character(LEN=H_LONG) :: rgnmngfname = '' !> region management file name

    namelist / ADMPARAM / &
        glevel,           &
        rlevel,           &
        vlayer,           &
        rgnmngfname,      &
        ADM_HGRID_SYSTEM, &
#ifndef _FIXEDINDEX_
        ADM_vlink,        &
#endif
        ADM_XTMS_MLCP_S,  &
        debug

    integer  :: nmax, dmd
    integer  :: l, rgnid
    integer  :: ierr
    !---------------------------------------------------------------------------

    ADM_prc_me = PRC_myrank + 1
    ADM_prc_pl = 1

#ifdef _FIXEDINDEX_
    if ( PRC_nprocs /= PRC_nprocs ) then
       write(*,*) 'xxx Fixed prc_all is not match (fixed,requested): ', PRC_nprocs, PRC_nprocs
       stop
    endif
#endif

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[adm]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=ADMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(*,*) 'xxx Not found namelist ADMPARAM! STOP.'
       call PRC_MPIstop
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=ADMPARAM)

    ! Error if glevel & rlevel are not defined
    if ( glevel < 1 ) then
       write(*,*) 'xxx [ADM_setup] glevel is not appropriate :', glevel
       call PRC_MPIstop
    endif
    if ( rlevel < 0 ) then
       write(*,*) 'xxx [ADM_setup] rlevel is not appropriate :', rlevel
       call PRC_MPIstop
    endif

    ! Error if glevel & rlevel are different from global parameter
    if ( GLOBAL_ensemble_on ) then
       if ( glevel /= GLOBAL_glevel ) then
          write(*,*) 'xxx glevel is not equal (global,local) :', GLOBAL_glevel, glevel
          call PRC_MPIstop
       endif
       if ( rlevel /= GLOBAL_rlevel ) then
          write(*,*) 'xxx rlevel is not equal (global,local) :', GLOBAL_rlevel, rlevel
          call PRC_MPIstop
       endif
    endif

#ifdef _FIXEDINDEX_
    if ( ADM_HGRID_SYSTEM == 'ICO' ) then
       dmd        = 10
    elseif( ADM_HGRID_SYSTEM == 'PERIODIC-1DMD' ) then ! T.Ohno 110721
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' ) then ! M.Hara 110721
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == 'ICO-XTMS' ) then
       dmd        = 10
    else
       write(*,*) 'xxx [ADM_setup] Not appropriate param for ADM_HGRID_SYSTEM. STOP.', trim(ADM_HGRID_SYSTEM)
       call PRC_MPIstop
    endif
#else
    if ( ADM_HGRID_SYSTEM == 'ICO' ) then
       ADM_vlink  = 5
       dmd        = 10
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'LCP' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = 4 * ADM_vlink
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'MLCP-OLD' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = 2 * ADM_vlink
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'MLCP' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = (1+ADM_XTMS_MLCP_S)  * ADM_vlink
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'PERIODIC-1DMD' ) then ! T.Ohno 110721
       ADM_vlink  = 5
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' ) then ! M.Hara 110721
       ADM_vlink  = 5
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == 'ICO-XTMS' ) then
       ADM_vlink  = 5
       dmd        = 10
       ADM_prc_pl = 1
    else
       write(*,*) 'xxx [ADM_setup] Not appropriate param for ADM_HGRID_SYSTEM. STOP.', trim(ADM_HGRID_SYSTEM)
       call PRC_MPIstop
    endif

    ADM_gall_pl = ADM_vlink + 1
    ADM_gmax_pl = ADM_vlink + 1
#endif

#ifdef _FIXEDINDEX_
    if ( ADM_glevel /= glevel ) then
       write(*,*) 'xxx [ADM_setup] Fixed glevel is not match (fixed,requested): ', ADM_glevel, glevel
       call PRC_MPIstop
    endif
    if ( ADM_rlevel /= rlevel ) then
       write(*,*) 'xxx [ADM_setup] Fixed rlevel is not match (fixed,requested): ', ADM_rlevel, rlevel
       call PRC_MPIstop
    endif
    if ( ADM_vlayer /= vlayer ) then
       write(*,*) 'xxx [ADM_setup] Fixed vlayer is not match (fixed,requested): ', ADM_vlayer, vlayer
       call PRC_MPIstop
    endif
    if ( ADM_DMD /= dmd ) then
       write(*,*) 'xxx [ADM_setup] Fixed dmd is not match (fixed,requested): ', ADM_DMD, dmd
       call PRC_MPIstop
    endif
#else
    ADM_glevel   = glevel
    ADM_rlevel   = rlevel
    ADM_vlayer   = vlayer
    ADM_DMD      = dmd

    ADM_rgn_nmax = 2**ADM_rlevel * 2**ADM_rlevel * ADM_DMD
    ADM_lall     = ADM_rgn_nmax / PRC_nprocs

    nmax         = 2**(ADM_glevel-ADM_rlevel)
    ADM_gall_1d  = 1 + nmax + 1
    ADM_gmin     = 1 + 1
    ADM_gmax     = 1 + nmax

    ADM_gall     = ( 1+nmax+1 ) * ( 1+nmax+1 )
    ADM_gall_in  = (   nmax+1 ) * (   nmax+1 )

    if ( ADM_vlayer == 1 ) then
       ADM_kall = 1
       ADM_kmin = 1
       ADM_kmax = 1
    else
       ADM_kall = 1 + ADM_vlayer + 1
       ADM_kmin = 1 + 1
       ADM_kmax = 1 + ADM_vlayer
    endif
#endif



    call RGNMNG_setup( rgnmngfname )



    allocate( GLOBAL_extension_rgn(ADM_lall) )
    allocate( ADM_have_sgp        (ADM_lall) )
    GLOBAL_extension_rgn(:) = ''
    ADM_have_sgp        (:) = .false.

    do l = 1, ADM_lall
       rgnid = RGNMNG_lp2r(l,ADM_prc_me)

       write(GLOBAL_extension_rgn(l),'(A,I5.5)') '.rgn', rgnid-1

       if ( RGNMNG_vert_num(I_W,rgnid) == 3 ) then
          ADM_have_sgp(l) = .true.
       endif
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       ADM_have_pl = .true.
    else
       ADM_have_pl = .false.
    endif

    ADM_l_me = 0

    call output_info

    return
  end subroutine ADM_setup

  !-----------------------------------------------------------------------------
  subroutine RGNMNG_setup( &
       rgnmngfname )
    use mod_process, only: &
       PRC_nprocs, &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: rgnmngfname

    namelist / RGNMNGPARAM / &
        RGNMNG_in_fname,  &
        RGNMNG_out_fname, &
        debug

    integer  :: l, p, r
    integer  :: ierr
    !---------------------------------------------------------------------------

    RGNMNG_in_fname = trim(rgnmngfname)

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[rgnmng]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=RGNMNGPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** RGNMNGPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist RGNMNGPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=RGNMNGPARAM)

    if ( ADM_lall > RGNMNG_llim ) then
       write(*,*) 'xxx limit exceed! local region: ', ADM_lall, RGNMNG_llim
       call PRC_MPIstop
    endif

    ! Global information (Each process has all the information)
    allocate( RGNMNG_edge_tab(I_RGNID:I_DIR,I_SW:I_SE,ADM_rgn_nmax) )
    allocate( RGNMNG_lnum    (PRC_nprocs) )
    allocate( RGNMNG_lp2r    (RGNMNG_llim,PRC_nprocs) )

    if ( RGNMNG_in_fname /= '' ) then
       call RGNMNG_input( RGNMNG_in_fname,        & ! [IN]
                          ADM_rgn_nmax,           & ! [IN]
                          PRC_nprocs,             & ! [IN]
                          RGNMNG_llim,            & ! [IN]
                          RGNMNG_edge_tab(:,:,:), & ! [OUT]
                          RGNMNG_lnum    (:),     & ! [OUT]
                          RGNMNG_lp2r    (:,:)    ) ! [OUT]
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** input file is not specified.'

       if ( ADM_HGRID_SYSTEM == 'ICO' ) then
          call RGNMNG_generate_ico( ADM_rlevel,             & ! [IN]
                                    ADM_rgn_nmax,           & ! [IN]
                                    PRC_nprocs,             & ! [IN]
                                    RGNMNG_llim,            & ! [IN]
                                    RGNMNG_edge_tab(:,:,:), & ! [OUT]
                                    RGNMNG_lnum    (:),     & ! [OUT]
                                    RGNMNG_lp2r    (:,:)    ) ! [OUT]
       else
          write(*,*) 'xxx auto generation of mnginfo is not supported for ADM_HGRID_SYSTEM = ', ADM_HGRID_SYSTEM
          call PRC_MPIstop
       endif
    endif

    if ( RGNMNG_out_fname /= '' ) then
       call RGNMNG_output( RGNMNG_out_fname,       & ! [IN]
                           ADM_rgn_nmax,           & ! [IN]
                           PRC_nprocs,             & ! [IN]
                           RGNMNG_llim,            & ! [IN]
                           RGNMNG_edge_tab(:,:,:), & ! [IN]
                           RGNMNG_lnum    (:),     & ! [IN]
                           RGNMNG_lp2r    (:,:)    ) ! [IN]
    endif

    if ( RGNMNG_lnum(ADM_prc_me) /= ADM_lall ) then
       write(*,*) 'xxx local region number is not match: ', RGNMNG_lnum(ADM_prc_me), ADM_lall
       call PRC_MPIstop
    endif

    if ( ADM_lall > RGNMNG_llim ) then
       write(*,*) 'xxx limit exceed! local region: ', ADM_lall, RGNMNG_llim
       call PRC_MPIstop
    endif

    !--- additional table (reversal,local)
    allocate( RGNMNG_r2lp(2,ADM_rgn_nmax) )
    allocate( RGNMNG_l2r (ADM_lall) )

    do p = 1, PRC_nprocs
    do l = 1, ADM_lall
       RGNMNG_r2lp(I_l,  RGNMNG_lp2r(l,p)) = l
       RGNMNG_r2lp(I_prc,RGNMNG_lp2r(l,p)) = p
    enddo
    enddo

    do l = 1, ADM_lall
       RGNMNG_l2r(l) = RGNMNG_lp2r(l,ADM_prc_me)
    enddo

    !--- region connection chains around the diamond vertexes
    allocate( RGNMNG_vert_num   ( I_W:I_S,ADM_rgn_nmax ) )
    allocate( RGNMNG_vert_tab   ( I_RGNID:I_DIR,   &
                                  I_W:I_S,         &
                                  ADM_rgn_nmax,    &
                                  ADM_vlink        ) )
    allocate( RGNMNG_vert_tab_pl( I_RGNID:I_DIR,   &
                                  ADM_rgn_nmax_pl, &
                                  ADM_vlink        ) )

    call RGNMNG_vertex_walkaround( ADM_rgn_nmax,             & ! [IN]
                                   ADM_rgn_nmax_pl,          & ! [IN]
                                   ADM_vlink,                & ! [IN]
                                   RGNMNG_edge_tab(:,:,:),   & ! [IN]
                                   RGNMNG_vert_num(:,:),     & ! [OUT]
                                   RGNMNG_vert_tab(:,:,:,:), & ! [OUT]
                                   RGNMNG_vert_tab_pl(:,:,:) ) ! [OUT]

    !--- tables for pole
    do r = 1, ADM_rgn_nmax
       if ( RGNMNG_vert_num(I_N,r) == ADM_vlink ) then
          RGNMNG_rgn4pl(I_NPL) = r
          exit
       endif
    enddo

    do r = 1, ADM_rgn_nmax
       if ( RGNMNG_vert_num(I_S,r) == ADM_vlink ) then
          RGNMNG_rgn4pl(I_SPL) = r
          exit
       endif
    enddo

    RGNMNG_r2p_pl(I_NPL) = ADM_prc_pl
    RGNMNG_r2p_pl(I_SPL) = ADM_prc_pl

    return
  end subroutine RGNMNG_setup

  !-----------------------------------------------------------------------------
  !> Input mnginfo file
  subroutine RGNMNG_input( &
       in_fname, &
       rall,     &
       pall,     &
       lall,     &
       edge_tab, &
       lnum,     &
       lp2r      )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in)  :: in_fname           !< input file
    integer,          intent(in)  :: rall               !< total number of region
    integer,          intent(in)  :: pall               !< total number of process
    integer,          intent(in)  :: lall               !< limit number of local region
    integer,          intent(out) :: edge_tab(2,4,rall) !< region link information (for 4 edges)
    integer,          intent(out) :: lnum(pall)         !< number of local region
    integer,          intent(out) :: lp2r(lall,pall)    !< l,prc => region

    integer  :: num_of_rgn !< number of region

    namelist / rgn_info / &
         num_of_rgn

    integer  :: rgnid      !< region ID
    integer  :: sw(2) = -1 !< south-west region info
    integer  :: nw(2) = -1 !< north-west region info
    integer  :: ne(2) = -1 !< north-east region info
    integer  :: se(2) = -1 !< south-east region info

    namelist / rgn_link_info / &
         rgnid, &
         sw,    &
         nw,    &
         ne,    &
         se

    integer  :: num_of_proc !< number of processes

    namelist / proc_info / &
         num_of_proc

    integer  :: peid            !< process ID
    integer  :: num_of_mng      !< number of regions be managed
    integer  :: mng_rgnid(RGNMNG_llim) !< managed region ID

    namelist / rgn_mng_info / &
         peid,       &
         num_of_mng, &
         mng_rgnid

    integer  :: fid, ierr
    integer  :: r, p, l
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** input mnginfo: ',trim(in_fname)

    fid = IO_get_available_fid()
    open( unit   = fid,            &
          file   = trim(in_fname), &
          form   = 'formatted',    &
          status = 'old',          &
          iostat = ierr            )

       ! ERROR if filename are not defined
       if ( ierr /= 0 ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx mnginfo file is not found! STOP. ', trim(in_fname)
          call PRC_MPIstop
       endif

       read(fid,nml=rgn_info)

       if ( num_of_rgn /= rall ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx No match for region number! STOP.'
          if( IO_L ) write(IO_FID_LOG,*) 'xxx rall= ',rall,' num_of_rgn=',num_of_rgn
          call PRC_MPIstop
       endif

       do r = 1, rall
          read(fid,nml=rgn_link_info)

          edge_tab(I_RGNID:I_DIR,I_SW,rgnid) = sw(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_NW,rgnid) = nw(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_NE,rgnid) = ne(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_SE,rgnid) = se(I_RGNID:I_DIR)
       enddo

       read(fid,nml=proc_info)
       if ( num_of_proc /= pall ) then
          if( IO_L ) write(IO_FID_LOG,*) ' xxx No match for  process number! STOP.'
          if( IO_L ) write(IO_FID_LOG,*) ' xxx prc_all= ',pall,' num_of_proc=',num_of_proc
          call PRC_MPIstop
       endif

       do p = 1, pall
          read(fid,nml=rgn_mng_info)

          lnum(peid) = num_of_mng

          if ( lnum(peid) /= ADM_lall ) then
             write(*,*) 'xxx local region number is not match: ', peid, lnum(peid), ADM_lall
             call PRC_MPIstop
          endif

          lp2r(:,peid) = -1 ! initialize
          do l = 1, lnum(peid)
             lp2r(l,peid) = mng_rgnid(l)
          enddo
       enddo

    close(fid)

    return
  end subroutine RGNMNG_input

  !-----------------------------------------------------------------------------
  !> Output mnginfo file
  subroutine RGNMNG_output( &
       out_fname, &
       rall,      &
       pall,      &
       lall,      &
       edge_tab,  &
       lnum,      &
       lp2r       )
    implicit none

    character(len=*), intent(in) :: out_fname          !< output file
    integer,          intent(in) :: rall               !< total number of region
    integer,          intent(in) :: pall               !< total number of process
    integer,          intent(in) :: lall               !< limit number of local region
    integer,          intent(in) :: edge_tab(2,4,rall) !< region link information (for 4 edges)
    integer,          intent(in) :: lnum(pall)         !< number of local region
    integer,          intent(in) :: lp2r(lall,pall)    !< l,prc => region

    integer  :: num_of_rgn !< number of region

    namelist / rgn_info / &
       num_of_rgn

    integer  :: rgnid      !< region ID
    integer  :: sw(2) = -1 !< south-west region info
    integer  :: nw(2) = -1 !< north-west region info
    integer  :: ne(2) = -1 !< north-east region info
    integer  :: se(2) = -1 !< south-east region info

    namelist / rgn_link_info / &
       rgnid, &
       sw,    &
       nw,    &
       ne,    &
       se

    integer  :: num_of_proc !< number of processes

    namelist / proc_info / &
       num_of_proc

    integer  :: peid            !< process ID
    integer  :: num_of_mng      !< number of regions be managed
    integer  :: mng_rgnid(RGNMNG_llim) !< managed region ID

    namelist / rgn_mng_info / &
       peid,       &
       num_of_mng, &
       mng_rgnid

    integer  :: fid
    integer  :: r, p, l
    !---------------------------------------------------------------------------

    if ( ADM_prc_me == ADM_prc_master ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** output mnginfo: ',trim(out_fname)

       fid = IO_get_available_fid()
       open( unit   = fid,             &
             file   = trim(out_fname), &
             form   = 'formatted'      )

          num_of_rgn = rall
          write(fid,nml=rgn_info)

          do r = 1, rall
             rgnid = r
             sw(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_SW,rgnid)
             nw(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_NW,rgnid)
             ne(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_NE,rgnid)
             se(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_SE,rgnid)

             write(fid,nml=rgn_link_info)
          enddo

          num_of_proc = pall
          write(fid,nml=proc_info)

          do p = 1, pall
             peid = p
             num_of_mng = lnum(peid)

             mng_rgnid(:) = -1
             do l = 1, lnum(peid)
                mng_rgnid(l) = lp2r(l,peid)
             enddo

             write(fid,nml=rgn_mng_info)
          enddo

       close(fid)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** output mnginfo: only in the master process'
    endif

    return
  end subroutine RGNMNG_output

  !-----------------------------------------------------------------------------
  !> Generate region management info
  Subroutine RGNMNG_generate_ico( &
       rlevel,   &
       rall,     &
       pall,     &
       lall,     &
       edge_tab, &
       lnum,     &
       lp2r      )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    integer, intent(in)  :: rlevel
    integer, intent(in)  :: rall
    integer, intent(in)  :: pall
    integer, intent(in)  :: lall
    integer, intent(out) :: edge_tab(2,4,rall)
    integer, intent(out) :: lnum(pall)
    integer, intent(out) :: lp2r(lall,pall)

    integer, parameter :: nmax_dmd = 10
    integer  :: dmd_data(4,nmax_dmd)

    integer  :: rall_1d, rall_1dmd

    integer  :: d_nb, i_nb, j_nb, rgnid_nb, direction
    integer  :: d, i, j, rgnid
    integer  :: l, p
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** generate region management table'
    if( IO_L ) write(IO_FID_LOG,*) '*** Topology: default icosahedral'

    dmd_data(:, 1) = (/  6, 5, 2,10 /)
    dmd_data(:, 2) = (/ 10, 1, 3, 9 /)
    dmd_data(:, 3) = (/  9, 2, 4, 8 /)
    dmd_data(:, 4) = (/  8, 3, 5, 7 /)
    dmd_data(:, 5) = (/  7, 4, 1, 6 /)
    dmd_data(:, 6) = (/  7, 5, 1,10 /)
    dmd_data(:, 7) = (/  8, 4, 5, 6 /)
    dmd_data(:, 8) = (/  9, 3, 4, 7 /)
    dmd_data(:, 9) = (/ 10, 2, 3, 8 /)
    dmd_data(:,10) = (/  6, 1, 2, 9 /)

    rall_1d   = 2**rlevel
    rall_1dmd = rall_1d*rall_1d

    !--- make region link table
    do d = 1, nmax_dmd
    do i = 1, rall_1d
    do j = 1, rall_1d
       rgnid = (d-1)*rall_1dmd + (j-1)*rall_1d + i

       !--- I_SW
       if ( j == 1 ) then
          if ( d <= 5 ) then
             i_nb = i
             j_nb = rall_1d
             d_nb = dmd_data(I_SW,d)
             direction = I_NE
          else
             i_nb = rall_1d
             j_nb = rall_1d+1-i
             d_nb = dmd_data(I_SW,d)
             direction = I_SE
          endif
       else
          i_nb = i
          j_nb = j-1
          d_nb = d
          direction = I_NE
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_SW,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_SW,rgnid) = direction

       !--- I_NW
       if ( i == 1 ) then
          if ( d <= 5 ) then
             i_nb = rall_1d+1-j
             j_nb = rall_1d
             d_nb = dmd_data(I_NW,d)
             direction = I_NE
          else
             i_nb = rall_1d
             j_nb = j
             d_nb = dmd_data(I_NW,d)
             direction = I_SE
          endif
       else
          i_nb = i-1
          j_nb = j
          d_nb = d
          direction = I_SE
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_NW,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_NW,rgnid) = direction

       !--- I_NE
       if ( j == rall_1d ) then
          if ( d <= 5 ) then
             i_nb = 1
             j_nb = rall_1d+1-i
             d_nb = dmd_data(I_NE,d)
             direction = I_NW
          else
             i_nb = i
             j_nb = 1
             d_nb = dmd_data(I_NE,d)
             direction = I_SW
          endif
       else
          i_nb = i
          j_nb = j+1
          d_nb = d
          direction = I_SW
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_NE,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_NE,rgnid) = direction

       !--- I_SE
       if ( i == rall_1d ) then
          if ( d <= 5 ) then
             i_nb = 1
             j_nb = j
             d_nb = dmd_data(I_SE,d)
             direction = I_NW
          else
             i_nb = rall_1d+1-j
             j_nb = 1
             d_nb = dmd_data(I_SE,d)
             direction = I_SW
          endif
       else
          i_nb = i+1
          j_nb = j
          d_nb = d
          direction = I_NW
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_SE,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_SE,rgnid) = direction

    enddo
    enddo
    enddo

    if ( mod(rall,pall) /= 0 ) then
       write(*,*) 'invalid number of process!', rall, pall
       call PRC_MPIstop
    else
       if ( rall/pall /= lall ) then
          write(*,*) 'xxx local region number is not match: ', rall/pall, lall
          call PRC_MPIstop
       endif

       do p = 1, pall
          lnum(p) = rall / pall
       enddo
    endif

    !--- make region-pe relationship
    lp2r(:,:) = -1
    rgnid = 0

    do p = 1, pall
    do l = 1, lnum(p)
       rgnid = rgnid + 1

       lp2r(l,p) = rgnid
    enddo
    enddo

    return
  end Subroutine RGNMNG_generate_ico

  !-----------------------------------------------------------------------------
  !> Search region with walking around the diamond vertexes
  subroutine RGNMNG_vertex_walkaround( &
       rall,       &
       rall_pl,    &
       vlink_nmax, &
       edge_tab,   &
       vert_num,   &
       vert_tab,   &
       vert_tab_pl )
    implicit none

    integer, intent(in)  :: rall
    integer, intent(in)  :: rall_pl
    integer, intent(in)  :: vlink_nmax
    integer, intent(in)  :: edge_tab(2,4,rall)

    integer, intent(out) :: vert_num   (4,rall)
    integer, intent(out) :: vert_tab   (2,4,rall,vlink_nmax)
    integer, intent(out) :: vert_tab_pl(2,rall_pl,vlink_nmax)

    integer  :: rgnid, dir
    integer  :: rgnid_next, dir_next

    integer  :: r, d, v
    !---------------------------------------------------------------------------

    vert_num(:,:)     = -1
    vert_tab(:,:,:,:) = -1

    do r = 1, rall
    do d = I_W, I_S

       rgnid = r
       select case(d)
       case(I_W)
          dir = I_SW
       case(I_N)
          dir = I_NW
       case(I_E)
          dir = I_NE
       case(I_S)
          dir = I_SE
       endselect

       v = 0
       do
          rgnid_next = edge_tab(I_RGNID,dir,rgnid)
          dir_next   = edge_tab(I_DIR,  dir,rgnid) - 1

          if( dir_next == 0 ) dir_next = 4
          v = v + 1
          vert_tab(I_RGNID,d,r,v) = rgnid
          vert_tab(I_DIR,  d,r,v) = dir

          rgnid = rgnid_next
          dir   = dir_next

          if( rgnid == r ) exit
       enddo
       vert_num(d,r) = v

    enddo
    enddo

    vert_tab_pl(:,:,:) = -1

    do r = 1, rall
       if ( vert_num(I_N,r) == vlink_nmax ) then
          do v = 1, vlink_nmax
             vert_tab_pl(I_RGNID,I_NPL,v) = vert_tab(I_RGNID,I_N,r,v)
             vert_tab_pl(I_DIR,  I_NPL,v) = vert_tab(I_DIR,  I_N,r,v)
          enddo
          exit
       endif
    enddo

    do r = 1, rall
       if ( vert_num(I_S,r) == vlink_nmax ) then
          do v = 1, vlink_nmax
             vert_tab_pl(I_RGNID,I_SPL,v) = vert_tab(I_RGNID,I_S,r,v)
             vert_tab_pl(I_DIR,  I_SPL,v) = vert_tab(I_DIR,  I_S,r,v)
          enddo
          exit
       endif
    enddo

    return
  end subroutine RGNMNG_vertex_walkaround

  !-----------------------------------------------------------------------------
  subroutine output_info
    use mod_process, only: &
       PRC_nprocs
    implicit none

    integer          :: rgnid, rgnid_next
    character(len=2) :: dstr, dstr_next

    integer  :: l, d, v
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '====== Process management info. ======'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Total number of process           : ', PRC_nprocs
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- My Process number = (my rank + 1) : ', ADM_prc_me
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '====== Region/Grid topology info. ======'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)' ) '--- Grid sysytem                      : ', ADM_HGRID_SYSTEM
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- #  of diamond                     : ', ADM_DMD
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '====== Region management info. ======'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Region level (RL)                 : ', ADM_rlevel
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7,3(A,I4),A)') '--- Total number of region            : ', ADM_rgn_nmax, &
                                             ' (', 2**ADM_rlevel, ' x', 2**ADM_rlevel, ' x', ADM_DMD, ' )'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- #  of region per process          : ', ADM_lall
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '--- ID of region in my process        : '
    if( IO_L ) write(IO_FID_LOG,*)           RGNMNG_lp2r(1:ADM_lall,ADM_prc_me)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Region ID, contains north pole    : ', RGNMNG_rgn4pl(I_NPL)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Region ID, contains south pole    : ', RGNMNG_rgn4pl(I_SPL)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Process rank, managing north pole : ', RGNMNG_r2p_pl(I_NPL)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Process rank, managing south pole : ', RGNMNG_r2p_pl(I_SPL)
    if( IO_L ) write(IO_FID_LOG,'(1x,A)'   ) '====== Grid management info. ======'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Grid level (GL)                   : ', ADM_glevel
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7,2(A,I4),A,I7,A)') '--- Total number of grid (horizontal) : ', &
                                                 4**(ADM_glevel-ADM_rlevel)*ADM_rgn_nmax,    &
                                                 ' (', 2**(ADM_glevel-ADM_rlevel),           &
                                                 ' x', 2**(ADM_glevel-ADM_rlevel),           &
                                                 ' x', ADM_rgn_nmax, ' )'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- Number of vertical layer          : ', ADM_kmax-ADM_kmin+1

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '====== region management information ======'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)') '--- # of region in this node : ', ADM_lall

       if( IO_L ) write(IO_FID_LOG,*) '--- (l,prc_me) => (rgn)'
       do l = 1, ADM_lall
          rgnid = RGNMNG_l2r(l)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,I6,A,I6,A)') '--- (',l,',',ADM_prc_me,') => (',rgnid,') '
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Link information'
       do l = 1, ADM_lall
          rgnid = RGNMNG_l2r(l)

          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '--- edge link: (rgn,direction)'
          do d = I_SW, I_SE
             rgnid_next = RGNMNG_edge_tab(I_RGNID,d,rgnid)
             dstr       = RGNMNG_edgename(d)
             dstr_next  = RGNMNG_edgename(RGNMNG_edge_tab(I_DIR,d,rgnid))
             if( IO_L ) write(IO_FID_LOG,'(5x,A,I6,A,A,A,I6,A,A,A)') '(',rgnid,',',dstr,') -> (', rgnid_next,',', dstr_next,')'
          enddo

          if( IO_L ) write(IO_FID_LOG,*) '--- vertex link: (rgn)'
          do d = I_W, I_S
             dstr = RGNMNG_vertname(d)
             if( IO_L ) write(IO_FID_LOG,'(5x,A,I6,A,A,A)',advance='no') '(',rgnid,',',dstr,')'
             do v = 2, RGNMNG_vert_num(d,rgnid)
                dstr = RGNMNG_vertname(RGNMNG_vert_tab(I_DIR,d,rgnid,v))
                if( IO_L ) write(IO_FID_LOG,'(A,I6,A,A,A)',advance='no') ' -> (',RGNMNG_vert_tab(I_RGNID,d,rgnid,v),',',dstr,')'
             enddo
             if( IO_L ) write(IO_FID_LOG,*)
          enddo
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Pole information (in the global scope)'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- region, having north pole data : ', RGNMNG_rgn4pl(I_NPL)
       if( IO_L ) write(IO_FID_LOG,*) '--- vertex link: (north pole)'
       do v = 2, ADM_vlink
          rgnid = RGNMNG_vert_tab_pl(I_RGNID,I_NPL,v)
          dstr  = RGNMNG_vertname(RGNMNG_vert_tab_pl(I_DIR,I_NPL,v))
          if( IO_L ) write(IO_FID_LOG,'(A,I6,A,A,A)',advance='no') ' -> (',rgnid,',',dstr,')'
       enddo
       rgnid = RGNMNG_vert_tab_pl(I_RGNID,I_NPL,1)
       dstr  = RGNMNG_vertname(RGNMNG_vert_tab_pl(I_DIR,I_NPL,1))
       if( IO_L ) write(IO_FID_LOG,'(A,I6,A,A,A)',advance='no') ' -> (',rgnid,',',dstr,')'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- process, managing north pole : ', RGNMNG_r2p_pl(I_NPL)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- region, having south pole data : ', RGNMNG_rgn4pl(I_SPL)
       if( IO_L ) write(IO_FID_LOG,*) '--- vertex link: (south pole)'
       do v = 2, ADM_vlink
          rgnid = RGNMNG_vert_tab_pl(I_RGNID,I_SPL,v)
          dstr  = RGNMNG_vertname(RGNMNG_vert_tab_pl(I_DIR,I_SPL,v))
          if( IO_L ) write(IO_FID_LOG,'(A,I6,A,A,A)',advance='no') ' -> (',rgnid,',',dstr,')'
       enddo
       rgnid = RGNMNG_vert_tab_pl(I_RGNID,I_SPL,1)
       dstr  = RGNMNG_vertname(RGNMNG_vert_tab_pl(I_DIR,I_SPL,1))
       if( IO_L ) write(IO_FID_LOG,'(A,I6,A,A,A)',advance='no') ' -> (',rgnid,',',dstr,')'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- process, managing south pole : ', RGNMNG_r2p_pl(I_SPL)
    endif

    return
  end subroutine output_info

  !-----------------------------------------------------------------------------
  subroutine ADM_make_idstr( &
       fname,    &
       besename, &
       l         )
    implicit none

    character(len=*), intent(out) :: fname    ! combined string (file name)
    character(len=*), intent(in)  :: besename ! basename
    integer,          intent(in)  :: l        ! local region number
    !---------------------------------------------------------------------------

    fname =  trim(GLOBAL_prefix_dir)      &
          // trim(besename)               &
          // trim(GLOBAL_extension_ens)   &
          // trim(GLOBAL_extension_rgn(l) )

    return
  end subroutine ADM_make_idstr

end module mod_adm
