!-------------------------------------------------------------------------------
!>
!! Administration module
!!
!! @par Description
!!         This module is for the management of process and region on
!!         the icosahedral grid configuration.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2007-10-22 (T.Mitsui)  change value of PRC_RGN_NMAX
!! @li      2008-01-30 (S.Iga)     private procedure mk_suffix is changed to public procedure
!! @li      2009-08-18 (T.Mitsui)  modify adm_proc_stop to keep out extra process from main routines.
!! @li      2010-04-26 (M.Satoh)   add ADM_l_me
!! @li      2010-06-07 (S.Iga)     new grid (Iga 2010) is implemented. (see string XTMS)
!! @li      2011-06-30 (T.Seiki)   fix undefined value (after, 07-10-22)
!! @li      2011-07-21 (T.Ohno)    2 new grid systems (1DMD-ON-SPHERE are added by Hara-san@JAMSTEC)
!! @li      2012-01-12 (H.Yashiro) add filename specification for logfile(optional)
!! @li      2012-06-11 (H.Yashiro) Milestone-project, code cleanup
!!
!<
module mod_adm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ADM_proc_init
  public :: ADM_proc_stop
  public :: ADM_setup
  public :: ADM_mk_suffix

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !------ Character length of system control
  integer, public, parameter :: ADM_NSYS = 32
  !
  !------ Maximum length of file name
  integer, public, parameter :: ADM_MAXFNAME = 128

  !
  !====== Basic definition & information ======
  !
  !------ Log file ID & Control file ID
  integer, public, parameter :: ADM_LOG_FID = 30
  integer, public, parameter :: ADM_CTL_FID = 35
  !
  !------ Identifier for single computation or parallel computation
  integer, public, parameter :: ADM_SINGLE_PRC = 0
  integer, public, parameter :: ADM_MULTI_PRC  = 1
  !
  !------ Identifiers of directions of region edges
  integer, public, parameter :: ADM_SW = 1
  integer, public, parameter :: ADM_NW = 2
  integer, public, parameter :: ADM_NE = 3
  integer, public, parameter :: ADM_SE = 4
  !
  !------ Identifiers of directions of region vertices
  integer, public, parameter :: ADM_W = 1
  integer, public, parameter :: ADM_N = 2
  integer, public, parameter :: ADM_E = 3
  integer, public, parameter :: ADM_S = 4
  !
  !--- Identifier of triangle element (i-axis-side or j-axis side)
  integer, public, parameter :: ADM_TI = 1
  integer, public, parameter :: ADM_TJ = 2
  !
  !--- Identifier of line element (i-axis-side, ij-axis side, or j-axis side)
  integer, public, parameter :: ADM_AI  = 1
  integer, public, parameter :: ADM_AIJ = 2
  integer, public, parameter :: ADM_AJ  = 3
  !
  !------ Identifier of 1 variable
  integer, public, parameter :: ADM_KNONE = 1
  integer, public, parameter :: ADM_VNONE = 1
  !
  !------ Identifier of poles (north pole or south pole)
  integer, public, parameter :: ADM_NPL = 1
  integer, public, parameter :: ADM_SPL = 2
  !
  !------ Fist colomn on the table for region and direction
  integer, public, parameter :: ADM_RID = 1
  integer, public, parameter :: ADM_DIR = 2
  !
  real(8), public, parameter :: ADM_VMISS = 1.D0

  !
  !====== Information for processes ======
  !
  !------ Communication world for NICAM
  integer, public, save      :: ADM_COMM_WORLD
  !
  !------ Master process
  integer, public, parameter :: ADM_prc_run_master = 1
  !
  !------ Total number of process
  integer, public, save      :: ADM_prc_all
  !
  !------ My process ID
  integer, public, save      :: ADM_prc_me
  !
  !------ Process ID which manages the pole regions.
  integer, public, save      :: ADM_prc_pl
  !
  !------ Process ID which have the pole regions.
  integer, public,  save     :: ADM_prc_npl
  integer, public,  save     :: ADM_prc_spl
  integer, public,  save     :: ADM_prc_nspl(ADM_NPL:ADM_SPL)

  logical, public,  save     :: ADM_have_pl

  !
  !====== Information for processes-region relationship ======
  !
  !------ Maximum number of regions managed by 1 process.
  integer, public,                parameter :: PRC_RGN_NMAX = 2560
  !
  !------ Regin managing file name
  character(len=ADM_MAXFNAME), public, save :: ADM_rgnmngfname
  !
  !------ Number of regions mangeged by each process
  integer, public, allocatable,        save :: ADM_prc_rnum(:)
  !
  !------ Table of regions managed by each process
  integer, public, allocatable,        save :: ADM_prc_tab(:,:)
  !
  !------ Table of edge link information
  integer, public, allocatable,        save :: ADM_rgn_etab(:,:,:)
  !<-----
  !<----- ADM_rgn_etab( ADM_RID:ADM_DIR, &
  !<-----               ADM_SW:ADM_SE,   &
  !<-----               ADM_rgn_nmax     )
  !<-----
  !
  !------ Table of process ID from region ID
  integer, public, allocatable,        save :: ADM_rgn2prc(:)
  !<-----
  !<----- ADM_rgn2prc(ADM_rgn_nmax)
  !<-----
  !
  !------ Maximum number of vertex linkage
  !integer, public, parameter :: ADM_VLINK_NMAX=5 ! S.Iga 100607
  integer, public,                     save :: ADM_VLINK_NMAX ! S.Iga 100607
  !
  !------ Table of n-vertex-link(?) at the region vertex
  integer, public, allocatable,        save :: ADM_rgn_vnum(:,:)
  !<-----
  !<----- ADM_rgn_vnum( ADM_W:ADM_S, &
  !<-----               ADM_rgn_nmax )
  !<-----
  !
  !------ Table of vertex link information
  integer, public, allocatable,        save :: ADM_rgn_vtab(:,:,:,:)
  !<-----
  !<----- ADM_rgn_vtab( ADM_RID:ADM_DIR, &
  !<-----               ADM_W:ADM_S,     &
  !<-----               ADM_rgn_nmax,    &
  !<-----               ADM_VLINK_NMAX   )
  !<-----
  !
  !------ Table of vertex link information for poles
  integer, public, allocatable,        save :: ADM_rgn_vtab_pl(:,:,:)
  !<-----
  !<----- ADM_rgn_vtab_pl( ADM_RID:ADM_DIR, &
  !<-----                  ADM_RGN_NMAX_PL, &
  !<-----                  ADM_VLINK_NMAX   )
  !<-----
  !
  !------ Region ID (reguler) of north pole management
  integer, public, save :: ADM_rgnid_npl_mng
  integer, public, save :: ADM_rgnid_spl_mng


  !
  !====== Information for regions ======
  !
  !------ Region division level
  integer, public, save      :: ADM_rlevel
  !
  !------ Total number of regular regions managed by all process
  integer, public, save      :: ADM_rgn_nmax
  !
  !------ Maximum number of pole regions
  integer, public, parameter :: ADM_rgn_nmax_pl = 2
  !
  !------ Local region number
  integer, public, save      :: ADM_lall
  !
  !------ Local region number for poles
  integer, public, save      :: ADM_lall_pl = ADM_rgn_nmax_pl
  !
  !------ Present Local region number ! 2010.4.26 M.Satoh
  integer, public, save      :: ADM_l_me

  logical, public, allocatable, save :: ADM_have_sgp(:) ! region have singlar point?

  !
  !====== Grid resolution informations  ======
  !
  !------ Grid division level
  integer, public, save      :: ADM_glevel
  !
  !------ Horizontal grid numbers
  integer, public, save      :: ADM_gmin
  integer, public, save      :: ADM_gmax
  integer, public, save      :: ADM_gall_1d
  integer, public, save      :: ADM_gall
  !
  !----- grid number of inner region in the diamond
  integer, public, save      :: ADM_gall_in
  !
  !------ Identifiers of grid points around poles.
  integer, public, parameter :: ADM_gslf_pl = 1
  integer, public, parameter :: ADM_gmin_pl = 2
  integer, public, save      :: ADM_gmax_pl     ! [mod] S.Iga 100607
  integer, public, save      :: ADM_gall_pl     ! [mod] S.Iga 100607
  !
  !------ Vertica grid numbers
  integer, public, save      :: ADM_vlayer
  integer, public, save      :: ADM_kmin
  integer, public, save      :: ADM_kmax
  integer, public, save      :: ADM_kall

  !
  !======  List vector for 1-dimensional array in the horiz. dir. ======
  !
  !------ Identifiers of grid points around a grid point
  integer, public, parameter :: ADM_GIJ_nmax = 7
  integer, public, parameter :: ADM_GIoJo = 1
  integer, public, parameter :: ADM_GIpJo = 2
  integer, public, parameter :: ADM_GIpJp = 3
  integer, public, parameter :: ADM_GIoJp = 4
  integer, public, parameter :: ADM_GImJo = 5
  integer, public, parameter :: ADM_GImJm = 6
  integer, public, parameter :: ADM_GIoJm = 7
  !
  !------ List vectors
  integer, public,              save :: ADM_IooJoo_nmax
  integer, public, allocatable, save :: ADM_IooJoo(:,:)
  !<-----
  !<----- ADM_IooJoo(ADM_IooJoo_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_IooJmo_nmax
  integer, public, allocatable, save :: ADM_IooJmo(:,:)
  !<-----
  !<----- ADM_IooJmo(ADM_IooJmo_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_IooJop_nmax
  integer, public, allocatable, save :: ADM_IooJop(:,:)
  !<-----
  !<----- ADM_IooJop(ADM_IooJop_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_IooJmp_nmax
  integer, public, allocatable, save :: ADM_IooJmp(:,:)
  !<-----
  !<----- ADM_IooJmp(ADM_IooJmp_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_ImoJoo_nmax
  integer, public, allocatable, save :: ADM_ImoJoo(:,:)
  !<-----
  !<----- ADM_ImoJoo(ADM_ImoJoo_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_ImoJmo_nmax
  integer, public, allocatable, save :: ADM_ImoJmo(:,:)
  !<-----
  !<----- ADM_ImoJmo(ADM_ImoJmo_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_ImoJop_nmax
  integer, public, allocatable, save :: ADM_ImoJop(:,:)
  !<-----
  !<----- ADM_ImoJop(ADM_ImoJop_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_ImoJmp_nmax
  integer, public, allocatable, save :: ADM_ImoJmp(:,:)
  !<-----
  !<----- ADM_ImoJmp(ADM_ImoJmp_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_IopJoo_nmax
  integer, public, allocatable, save :: ADM_IopJoo(:,:)
  !<-----
  !<----- ADM_IopJoo(ADM_IopJoo_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_IopJmo_nmax
  integer, public, allocatable, save :: ADM_IopJmo(:,:)
  !<-----
  !<----- ADM_IopJmo(ADM_IopJmo_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_IopJop_nmax
  integer, public, allocatable, save :: ADM_IopJop(:,:)
  !<-----
  !<----- ADM_IopJop(ADM_IopJop_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_IopJmp_nmax
  integer, public, allocatable, save :: ADM_IopJmp(:,:)
  !<-----
  !<----- ADM_IopJmp(ADM_IopJmp_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_ImpJoo_nmax
  integer, public, allocatable, save :: ADM_ImpJoo(:,:)
  !<-----
  !<----- ADM_ImpJoo(ADM_ImpJoo_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_ImpJmo_nmax
  integer, public, allocatable, save :: ADM_ImpJmo(:,:)
  !<-----
  !<----- ADM_ImpJmo(ADM_ImpJmo_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_ImpJop_nmax
  integer, public, allocatable, save :: ADM_ImpJop(:,:)
  !<-----
  !<----- ADM_ImpJop(ADM_ImpJop_nmax,ADM_GIJ_nmax)
  !<-----
  integer, public,              save :: ADM_ImpJmp_nmax
  integer, public, allocatable, save :: ADM_ImpJmp(:,:)
  !<-----
  !<----- ADM_ImpJmp(ADM_ImpJmp_nmax,ADM_GIJ_nmax)
  !<-----

  !=========== For New Grid (XTMS) start    <= S.Iga100607
  !
  !------ Horizontal Grid type
  character(len=ADM_MAXFNAME), public, save  :: ADM_HGRID_SYSTEM = 'ICO' ! icosahedral
  !                                             'ICO-XTMS' icosahedral but XTMS is used in oprt
  !                                             'LCP'      Lambert-cornial (including PSP)
  !                                             'MLCP'     Mercator+Lambert-cornial
  !                                             'MLCP-OLD' OLD vergion (only for s=1)
  !
  !------ Number of lines at each pole (maybe this is identical to ADM_VLINK_NMAX)
  integer,                     public, save  :: ADM_XTMS_K=-1 ! default
  !                                             ICO:5
  !                                             PSP:6
  !                                             LCP, MLCP:k
  !
  !------ Number of segment for MLCP
  integer,                     public, save  :: ADM_XTMS_MLCP_S= 1
  !
  !------ XTMS LEVEL (it is conveniently defined especially for mod_oprt)
  integer,                     public, save  :: ADM_XTMS_LEVEL = 0 ! original icosahedral (NICAM)
  !                                             = 1 ! XTMS level 1
  !                                             = 2 ! XTMS level 2 (to be implemented)
  !=========== For New Grid (XTMS) end    S.Iga100607 =>

  logical,                     public, save  :: ADM_debug = .false. ! [ADD] H.Yashiro 20120703

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: input_mnginfo
  private :: output_info
  private :: setup_vtab

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: GDUMMY = 1 ! Horizontal dummy(halo) cell
  integer, private, parameter :: KDUMMY = 1 ! Vertical   dummy(halo) cell

  integer, private, save      :: ADM_run_type    ! Run type (single or multi processes)

  integer, private, save      :: NMAX_DMD = -999 ! number of diamond

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine ADM_proc_init
  !>
  subroutine ADM_proc_init( rtype )
    implicit none

    integer, intent(in) :: rtype ! multi or single processor?

    integer :: my_rank
    integer :: ierr
    !---------------------------------------------------------------------------

    ADM_run_type = rtype

    if ( rtype == ADM_MULTI_PRC ) then

       call MPI_Init(ierr)
       call MPI_Comm_size(MPI_COMM_WORLD, ADM_prc_all, ierr)
       call MPI_Comm_rank(MPI_COMM_WORLD, my_rank,     ierr)

       call MPI_Comm_split(MPI_COMM_WORLD, 0, my_rank, ADM_COMM_WORLD,ierr)

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
    else
       ADM_prc_all = 1
       my_rank     = 0
    endif

    ADM_prc_me = my_rank + 1
    ADM_prc_pl = 1

    return
  end subroutine ADM_proc_init

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine ADM_proc_stop
  !>
  subroutine ADM_proc_stop
    implicit none

    character(len=ADM_NSYS) :: request
    integer                 :: ierr
    !---------------------------------------------------------------------------

    ! flush 1kbyte
    write(ADM_LOG_FID,'(32A32)') '                                '

    if ( ADM_run_type == ADM_MULTI_PRC ) then
       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) 'MPI process going to STOP...'

       request='STOP'
       call MPI_BCAST( request,              & !--- starting address
                       ADM_NSYS,             & !--- number of array
                       MPI_CHARACTER,        & !--- type
                       ADM_prc_run_master-1, & !--- source rank
                       MPI_COMM_WORLD,       & !--- world
                       ierr                  ) !--- error id

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       write(ADM_LOG_FID,*) 'MPI process has normally finished.'
       write(ADM_LOG_FID,*) '############################################################'
       call MPI_Finalize(ierr)

       close(ADM_CTL_FID)
       close(ADM_LOG_FID)

       stop
    else
       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) 'Serial process stopeed.'
       write(ADM_LOG_FID,*) '############################################################'

       close(ADM_CTL_FID)
       close(ADM_LOG_FID)

       stop
    endif

    return
  end subroutine ADM_proc_stop

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine ADM_setup
  !>
  subroutine ADM_setup( &
       param_fname, &
       msg_base     )
    use mod_misc, only: &
       MISC_make_idstr
    implicit none

    character(LEN=*), intent(in) :: param_fname ! namelist file name

    character(len=*), intent(in), optional :: msg_base ! output file for msg.pexxxxx file

    integer                     :: glevel      = -1
    integer                     :: rlevel      = -1
    integer                     :: vlayer      =  1
    character(LEN=ADM_MAXFNAME) :: rgnmngfname = ''

    namelist / ADMPARAM / &
        glevel,           & !--- grid division level
        rlevel,           & !--- region division level
        vlayer,           & !--- number of inner vertical layer
        rgnmngfname,      & !--- region management file name
        ADM_HGRID_SYSTEM, & !--- grid system (default: ico)  ! S.Iga100607
        ADM_XTMS_K,       & !--- num of lines at PL          ! S.Iga100607
        ADM_XTMS_MLCP_S,  & !--- num of segment for MLCP     ! S.Iga100607
        ADM_debug

    integer :: rgn_nmax
    integer :: nmax
    integer :: l, rgnid
    integer :: ierr

    character(LEN=ADM_MAXFNAME) :: fname
    character(LEN=ADM_MAXFNAME) :: msg
    !---------------------------------------------------------------------------

    msg = 'msg'
    if( present(msg_base) ) msg = msg_base ! [add] H.Yashiro 20110701

    !--- open message file
    call MISC_make_idstr(fname,trim(msg),'pe',ADM_prc_me)
    open( unit = ADM_LOG_FID, &
          file = trim(fname), &
          form = 'formatted'  )

    write(ADM_LOG_FID,*) '############################################################'
    write(ADM_LOG_FID,*) '#                                                          #'
    write(ADM_LOG_FID,*) '#   NICAM : Nonhydrostatic ICosahedal Atmospheric Model    #'
    write(ADM_LOG_FID,*) '#                                                          #'
    write(ADM_LOG_FID,*) '############################################################'

    !--- open control file
    open( unit   = ADM_CTL_FID,       &
          file   = trim(param_fname), &
          form   = 'formatted',       &
          status = 'old',             &
          iostat = ierr               )

    if ( ierr /= 0 ) then
       write(*,*) 'xxx Cannot open parameter control file!'
       write(*,*) 'xxx filename:', trim(param_fname)
       call ADM_proc_stop
    endif

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[adm]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=ADMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(*,          *) 'xxx Not found namelist! STOP.'
       write(ADM_LOG_FID,*) 'xxx Not found namelist! STOP.'
       call ADM_proc_stop
    elseif ( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=ADMPARAM)

    ADM_glevel      = glevel
    ADM_rlevel      = rlevel
    ADM_vlayer      = vlayer
    ADM_rgnmngfname = trim(rgnmngfname)

    ! S.Iga 100607 start =>
    if ( trim(ADM_HGRID_SYSTEM) == 'ICO' ) then
       ADM_XTMS_level = 0
       ADM_XTMS_K     = 5
       NMAX_DMD       = 10
    elseif( trim(ADM_HGRID_SYSTEM) == 'LCP' ) then
       if( ADM_XTMS_K == -1 ) ADM_XTMS_K = 6
       ADM_XTMS_level = 1
       NMAX_DMD       = 4* ADM_XTMS_K
    elseif( trim(ADM_HGRID_SYSTEM) == 'MLCP-OLD' ) then
       if( ADM_XTMS_K == -1 ) ADM_XTMS_K = 6
       ADM_XTMS_level = 1
       NMAX_DMD       = 2* ADM_XTMS_K
    elseif( trim(ADM_HGRID_SYSTEM) == 'MLCP' ) then
       if( ADM_XTMS_K == -1 ) ADM_XTMS_K = 6
       ADM_XTMS_level = 1
       NMAX_DMD       = (1+ADM_XTMS_MLCP_S)  * ADM_XTMS_K
    elseif( trim(ADM_HGRID_SYSTEM) == 'PERIODIC-1DMD' ) then ! T.Ohno 110721
       ADM_XTMS_level = 0
       ADM_XTMS_K     = 5
       NMAX_DMD       = 1
       ADM_prc_pl     = -999
    elseif( trim(ADM_HGRID_SYSTEM) == '1DMD-ON-SPHERE' ) then ! M.Hara 110721
       ADM_XTMS_level = 0
       ADM_XTMS_K     = 5
       NMAX_DMD       = 1
       ADM_prc_pl     = -999
    elseif( trim(ADM_HGRID_SYSTEM) == 'ICO-XTMS' ) then
       ADM_XTMS_level = 1
       ADM_XTMS_K     = 5
       NMAX_DMD       = 10
    else
       write(*          ,*) 'xxx Name of ADM_HGRID_SYSTEM is wrong. STOP.'
       write(ADM_LOG_FID,*) 'xxx Name of ADM_HGRID_SYSTEM is wrong. STOP.'
       call ADM_proc_stop
    endif

    ADM_VLINK_NMAX = ADM_XTMS_K
    ADM_GMAX_PL    = ADM_VLINK_NMAX + 1
    ADM_GALL_PL    = ADM_VLINK_NMAX + 1
    ! <= S.Iga 100607 end

    ! ERROR if Glevel & Rlevel are not defined
    if ( ADM_glevel < 1 ) then
       write(*          ,*) 'xxx Glevel is not appropriate, STOP. GL=', ADM_glevel
       write(ADM_LOG_FID,*) 'xxx Glevel is not appropriate, STOP. GL=', ADM_glevel
       call ADM_proc_stop
    endif
    if ( ADM_rlevel < 0 ) then
       write(*          ,*) 'xxx Rlevel is not appropriate, STOP. RL=', ADM_rlevel
       write(ADM_LOG_FID,*) 'xxx Rlevel is not appropriate, STOP. RL=', ADM_rlevel
       call ADM_proc_stop
    endif

    rgn_nmax     = 2**ADM_rlevel
    ADM_rgn_nmax = rgn_nmax * rgn_nmax * NMAX_DMD

    call input_mnginfo( ADM_rgnmngfname )

    ADM_prc_npl = ADM_prc_pl
    ADM_prc_spl = ADM_prc_pl

    ADM_prc_nspl(ADM_NPL) = ADM_prc_npl
    ADM_prc_nspl(ADM_SPL) = ADM_prc_spl

    nmax        = 2**( ADM_glevel - ADM_rlevel )
    ADM_gmin    = GDUMMY + 1
    ADM_gmax    = GDUMMY + nmax
    ADM_gall_1d = GDUMMY + nmax + GDUMMY
    ADM_gall    = ADM_gall_1d * ADM_gall_1d

    ADM_gall_in = ( nmax+GDUMMY ) * ( nmax+GDUMMY ) !--- inner grid number (e.g., 33x33 for gl05)

    if ( ADM_vlayer == 1 ) then
       ADM_kmin = 1
       ADM_kmax = 1
       ADM_kall = 1
    else
       ADM_kmin = KDUMMY + 1
       ADM_kmax = KDUMMY + ADM_vlayer
       ADM_kall = KDUMMY + ADM_vlayer + KDUMMY
    endif

    ADM_lall = ADM_prc_rnum(ADM_prc_me)

    allocate( ADM_have_sgp(ADM_lall) )
    ADM_have_sgp(:) = .false.

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)
       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          ADM_have_sgp(l) = .true.
       endif
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       ADM_have_pl = .true.
    else
       ADM_have_pl = .false.
    endif

    ! 2010.4.26 M.Satoh; 2010.5.11 M.Satoh
    ! ADM_l_me: this spans from 1 to ADM_lall, if effective.
    ! Otherwise, ADM_l_me = 0 should be set. see mod_history
    ADM_l_me = 0

    !--- make suffix for list-vector loop.
    call ADM_mk_suffix

    call output_info

    return
  end subroutine ADM_setup

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine input_mnginfo
  !>
  subroutine input_mnginfo( fname )
    use mod_misc,  only :&
       MISC_get_available_fid
    implicit none

    character(len=ADM_MAXFNAME), intent(in) :: fname

    integer :: num_of_rgn !--- number of region

    namelist / rgn_info / &
         num_of_rgn

    integer :: rgnid                    !--- region ID
    integer :: sw(ADM_RID:ADM_DIR) = -1 !--- south-west region info
    integer :: nw(ADM_RID:ADM_DIR) = -1 !--- nouth-west region info
    integer :: ne(ADM_RID:ADM_DIR) = -1 !--- nouth-east region info
    integer :: se(ADM_RID:ADM_DIR) = -1 !--- south-east region info

    namelist / rgn_link_info / &
         rgnid, &
         sw,    &
         nw,    &
         ne,    &
         se

    integer :: num_of_proc !--- number of run-processes

    namelist /proc_info/ &
         num_of_proc

    integer :: peid                         !--- process ID
    integer :: num_of_mng                   !--- number of regions be managed
    integer :: mng_rgnid(PRC_RGN_NMAX) = -1 !--- managed region ID

    namelist /rgn_mng_info/ &
         peid,       &
         num_of_mng, &
         mng_rgnid

    integer :: fid, ierr
    integer :: l, m, n
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[mnginfo]/Category[common share]'

    fid = MISC_get_available_fid()
    open( unit   = fid,         &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

    !=> [add] H.Yashiro 20120611
    ! ERROR if filename are not defined
    if ( ierr /= 0 ) then
       write(ADM_LOG_FID,*) 'xxx mnginfo file is not found! STOP. ', trim(fname)
       call ADM_proc_stop
    endif
    !<= [add] H.Yashiro 20120611

    read(fid,nml=rgn_info)
    if ( num_of_rgn /= ADM_rgn_nmax ) then
       write(ADM_LOG_FID,*) 'xxx No match for region number! STOP.'
       write(ADM_LOG_FID,*) 'xxx ADM_rgn_nmax= ',ADM_rgn_nmax,' num_of_rgn=',num_of_rgn
       call ADM_proc_stop
    endif

    allocate( ADM_rgn_etab( ADM_RID:ADM_DIR, &
                            ADM_SW:ADM_SE,   &
                            ADM_rgn_nmax     ) )

    do l = 1, ADM_rgn_nmax
       read(fid,nml=rgn_link_info)

       ADM_rgn_etab(ADM_RID:ADM_DIR,ADM_SW,rgnid) = sw(ADM_RID:ADM_DIR)
       ADM_rgn_etab(ADM_RID:ADM_DIR,ADM_NW,rgnid) = nw(ADM_RID:ADM_DIR)
       ADM_rgn_etab(ADM_RID:ADM_DIR,ADM_NE,rgnid) = ne(ADM_RID:ADM_DIR)
       ADM_rgn_etab(ADM_RID:ADM_DIR,ADM_SE,rgnid) = se(ADM_RID:ADM_DIR)
    enddo

    read(fid,nml=proc_info)
    if ( ADM_prc_all /= num_of_proc ) then
       write(ADM_LOG_FID,*) ' xxx No match for  process number! STOP.'
       write(ADM_LOG_FID,*) ' xxx ADM_prc_all= ',ADM_prc_all,' num_of_proc=',num_of_proc
       call ADM_proc_stop
    endif

    if ( ADM_prc_all /= num_of_proc ) then
       write(ADM_LOG_FID,*) 'Msg : Sub[ADM_input_mngtab]/Mod[admin]'
       write(ADM_LOG_FID,*) ' --- No match for process number!'
       call ADM_proc_stop
    endif

    allocate( ADM_prc_rnum(ADM_prc_all)              )
    allocate( ADM_prc_tab (PRC_RGN_NMAX,ADM_prc_all) )
    allocate( ADM_rgn2prc (ADM_rgn_nmax)             )
    ADM_prc_tab = -1 ! [Fix] 11/06/30  T.Seiki, fill undefined value

    do m = 1, ADM_prc_all
       read(fid,nml=rgn_mng_info)

       ADM_prc_rnum(m)      = num_of_mng
       ADM_prc_tab (:,peid) = mng_rgnid(:)
       do n = 1, num_of_mng
          ADM_rgn2prc(mng_rgnid(n)) = peid
       enddo
    enddo

    call setup_vtab

    close(fid)

    return
  end subroutine input_mnginfo

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine setup_vtab
  !>
  subroutine setup_vtab
    implicit none

    integer :: nrid(ADM_VLINK_NMAX)
    integer :: nvid(ADM_VLINK_NMAX)
    integer :: vnum

    integer :: l, k, ll, v
    !---------------------------------------------------------------------------

    allocate( ADM_rgn_vnum( ADM_W:ADM_S, &
                            ADM_rgn_nmax ) )

    allocate( ADM_rgn_vtab( ADM_RID:ADM_DIR,&
                            ADM_W:ADM_S,    &
                            ADM_rgn_nmax,   &
                            ADM_VLINK_NMAX  ) )

    allocate( ADM_rgn_vtab_pl( ADM_RID:ADM_DIR, &
                               ADM_RGN_NMAX_PL, &
                               ADM_VLINK_NMAX   ) )

    do l = 1, ADM_rgn_nmax
       do k = ADM_W, ADM_S
          call set_vinfo(vnum,nrid,nvid,l,k)

          ADM_rgn_vnum(k,l)           = vnum
          ADM_rgn_vtab(ADM_RID,k,l,:) = nrid(:)
          ADM_rgn_vtab(ADM_DIR,k,l,:) = nvid(:)
       enddo
    enddo

    do l = 1, ADM_rgn_nmax
       if ( ADM_rgn_vnum(ADM_N,l) == ADM_VLINK_NMAX ) then
          ll = l
          exit
       endif
    enddo
    ADM_rgnid_npl_mng = ll

    do v = 1, ADM_VLINK_NMAX
       ADM_rgn_vtab_pl(ADM_RID,ADM_NPL,v) = ADM_rgn_vtab(ADM_RID,ADM_N,ll,v)
       ADM_rgn_vtab_pl(ADM_DIR,ADM_NPL,v) = ADM_rgn_vtab(ADM_DIR,ADM_N,ll,v)
    enddo

    do l = 1, ADM_rgn_nmax
       if ( ADM_rgn_vnum(ADM_S,l) == ADM_VLINK_NMAX ) then
          ll = l
          exit
       endif
    enddo
    ADM_rgnid_spl_mng = ll

    do v = 1, ADM_VLINK_NMAX
       ADM_rgn_vtab_pl(ADM_RID,ADM_SPL,v) = ADM_rgn_vtab(ADM_RID,ADM_S,ll,v)
       ADM_rgn_vtab_pl(ADM_DIR,ADM_SPL,v) = ADM_rgn_vtab(ADM_DIR,ADM_S,ll,v)
    enddo

    return
  end subroutine setup_vtab

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine set_vinfo
  !>
  subroutine set_vinfo( vert_num, nrgnid, nvertid, rgnid, vertid )
    implicit none

    integer,intent(out) :: vert_num
    integer,intent(out) :: nrgnid (:)
    integer,intent(out) :: nvertid(:)
    integer,intent(in)  :: rgnid
    integer,intent(in)  :: vertid

    integer :: eid, rid
    integer :: eid_new, rid_new
    !---------------------------------------------------------------------------

    vert_num = 0

    rid = rgnid
    eid = vertid
    select case(vertid)
    case(ADM_W)
       eid = ADM_SW
    case(ADM_N)
       eid = ADM_NW
    case(ADM_E)
       eid = ADM_NE
    case(ADM_S)
       eid = ADM_SE
    endselect

    nvertid(:) = -1
    nrgnid (:) = -1
    do
       rid_new = ADM_rgn_etab(ADM_RID,eid,rid)
       eid_new = ADM_rgn_etab(ADM_DIR,eid,rid) - 1

       if( eid_new == 0 ) eid_new = 4
       rid = rid_new
       eid = eid_new

       vert_num = vert_num + 1

       nrgnid (vert_num) = rid
       nvertid(vert_num) = eid

       if( rid == rgnid ) exit
    enddo

    return
  end subroutine set_vinfo

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine ADM_mk_suffix
  !>
  subroutine ADM_mk_suffix
    implicit none

    integer :: gall_in
    integer :: i, j, n
    !---------------------------------------------------------------------------

    gall_in = ADM_gmax-ADM_gmin+1

    !--- ADM_IooJoo
    ADM_IooJoo_nmax = ( gall_in ) * ( gall_in )
    allocate( ADM_IooJoo(ADM_IooJoo_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       ADM_IooJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJmo
    ADM_IooJmo_nmax = ( gall_in ) * ( gall_in+1 )
    allocate( ADM_IooJmo(ADM_IooJmo_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin,   ADM_gmax
       ADM_IooJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJop
    ADM_IooJop_nmax = ( gall_in ) * ( gall_in+1 )
    allocate( ADM_IooJop(ADM_IooJop_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin, ADM_gmax+1
    do i = ADM_gmin, ADM_gmax
       ADM_IooJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJmp
    ADM_IooJmp_nmax = ( gall_in ) * ( gall_in+2 )
    allocate( ADM_IooJmp(ADM_IooJmp_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin,   ADM_gmax
       ADM_IooJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJoo
    ADM_ImoJoo_nmax = ( gall_in+1 ) * ( gall_in )
    allocate( ADM_ImoJoo(ADM_ImoJoo_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin,   ADM_gmax
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJmo
    ADM_ImoJmo_nmax = ( gall_in+1 ) * ( gall_in+1 )
    allocate( ADM_ImoJmo(ADM_ImoJmo_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJop
    ADM_ImoJop_nmax = ( gall_in+1 ) * ( gall_in+1 )
    allocate( ADM_ImoJop(ADM_ImoJop_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin,   ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJmp
    ADM_ImoJmp_nmax = ( gall_in+1 ) * ( gall_in+2 )
    allocate( ADM_ImoJmp(ADM_ImoJmp_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJoo
    ADM_IopJoo_nmax = ( gall_in+1 ) * ( gall_in )
    allocate( ADM_IopJoo(ADM_IopJoo_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax+1
       ADM_IopJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJmo
    ADM_IopJmo_nmax = ( gall_in+1 ) * ( gall_in+1 )
    allocate( ADM_IopJmo(ADM_IopJmo_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin,   ADM_gmax+1
       ADM_IopJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJop
    ADM_IopJop_nmax = ( gall_in+1 ) * ( gall_in+1 )
    allocate( ADM_IopJop(ADM_IopJop_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin, ADM_gmax+1
    do i = ADM_gmin, ADM_gmax+1
       ADM_IopJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJmp
    ADM_IopJmp_nmax = ( gall_in+1 ) * ( gall_in+2 )
    allocate( ADM_IopJmp(ADM_IopJmp_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin, ADM_gmax+1
       ADM_IopJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJoo
    ADM_ImpJoo_nmax = ( gall_in+2 ) * ( gall_in )
    allocate( ADM_ImpJoo(ADM_ImpJoo_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin,   ADM_gmax
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJmo
    ADM_ImpJmo_nmax = ( gall_in+2 ) * ( gall_in+1 )
    allocate( ADM_ImpJmo(ADM_ImpJmo_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJop
    ADM_ImpJop_nmax = ( gall_in+2 ) * ( gall_in+1 )
    allocate( ADM_ImpJop(ADM_ImpJop_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin,   ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJmp
    ADM_ImpJmp_nmax = ( gall_in+2 ) * ( gall_in+2 )
    allocate( ADM_ImpJmp(ADM_ImpJmp_nmax,ADM_GIJ_nmax) )
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    return
  contains
    !---------------------------------------------------------------------------
    integer function suf(i,j)
      implicit none

      integer :: i, j
      !-------------------------------------------------------------------------

      suf = ADM_gall_1d * (j-1) + i

    end function suf

  end subroutine ADM_mk_suffix

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine output_info
  !>
  subroutine output_info
    implicit none

    integer :: n, k, m
    integer :: rgnid
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,'(1x,A)'   ) '====== Process management info. ======'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Total number of process           : ', ADM_prc_all
    write(ADM_LOG_FID,'(1x,A,I7)') '--- My Process rank                   : ', ADM_prc_me
    write(ADM_LOG_FID,'(1x,A)'   ) '====== Region/Grid topology info. ======'
    write(ADM_LOG_FID,'(1x,A,A)' ) '--- Grid sysytem                      : ', trim(ADM_HGRID_SYSTEM)
    write(ADM_LOG_FID,'(1x,A,I7)') '--- #  of diamond                     : ', NMAX_DMD
    write(ADM_LOG_FID,'(1x,A)'   ) '====== Region management info. ======'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Region level (RL)                 : ', ADM_rlevel
    write(ADM_LOG_FID,'(1x,A,I7,3(A,I4),A)') '--- Total number of region            : ', ADM_rgn_nmax, &
                                             ' (', 2**ADM_rlevel, ' x', 2**ADM_rlevel, ' x', NMAX_DMD, ' )'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- #  of region per process          : ', ADM_lall
    write(ADM_LOG_FID,'(1x,A)'   ) '--- ID of region in my process        : '
    write(ADM_LOG_FID,*) ADM_prc_tab(1:ADM_lall, ADM_prc_me)

    write(ADM_LOG_FID,'(1x,A,I7)') '--- Region ID, contains north pole    : ', ADM_rgnid_npl_mng
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Region ID, contains south pole    : ', ADM_rgnid_spl_mng
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Process rank, managing north pole : ', ADM_prc_npl
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Process rank, managing south pole : ', ADM_prc_spl
    write(ADM_LOG_FID,'(1x,A)'   ) '====== Grid management info. ======'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Grid level (GL)                   : ', ADM_glevel
    write(ADM_LOG_FID,'(1x,A,I7,2(A,I4),A,I7,A)') '--- Total number of grid (horizontal) : ',  &
                                                  4**(ADM_glevel-ADM_rlevel)*ADM_rgn_nmax, &
                                                  ' (', 2**(ADM_glevel-ADM_rlevel),         &
                                                  ' x', 2**(ADM_glevel-ADM_rlevel),         &
                                                  ' x', ADM_rgn_nmax, ' )'
    write(ADM_LOG_FID,'(1x,A,I7)') '--- Number of vertical layer          : ', ADM_kmax-ADM_kmin+1

    if ( ADM_debug ) then
       do n = 1, ADM_lall
          rgnid = ADM_prc_tab(n, ADM_prc_me)
          write(ADM_LOG_FID,*) ' --- Link information for region', rgnid

          write(ADM_LOG_FID,*) '     < edge link >   --- ( rgnid , edgid )'
          do k = ADM_SW, ADM_SE
             write(ADM_LOG_FID,*) '     (',rgnid,',',k,') -> ',         &
                                  '(', ADM_rgn_etab(ADM_RID,k,rgnid),   &
                                  ',', ADM_rgn_etab(ADM_DIR,k,rgnid), ')'
          enddo

          write(ADM_LOG_FID,*) '     < vertex link > --- ( rgnid , edgid )'
          do k = ADM_W, ADM_S
             write(ADM_LOG_FID,*) '     (',rgnid,',',k,') : ', ADM_rgn_vnum(k,rgnid), 'point link'
             do m = 1, ADM_rgn_vnum(k,rgnid)
                write(ADM_LOG_FID,*) '                -> ',                  &
                                     '(', ADM_rgn_vtab(ADM_RID,k,rgnid,m),   &
                                     ',', ADM_rgn_vtab(ADM_DIR,k,rgnid,m), ')'
             enddo
          enddo

       enddo

       write(ADM_LOG_FID,*) ' --- Table of corresponding between region ID and process ID'
       write(ADM_LOG_FID,*) '    region ID :  process ID'
       do n = 1, ADM_rgn_nmax
          write(ADM_LOG_FID,'(I13,I14)') n, ADM_rgn2prc(n)
       enddo
    endif

    return
  end subroutine output_info

end module mod_adm
!-------------------------------------------------------------------------------
