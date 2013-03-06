!-------------------------------------------------------------------------------
!>
!! Prognostic variable module
!!
!! @par Description
!!         This module contains the prognostic variables
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)   Imported from igdc-4.34
!! @li      2006-05-17 (H.Tomita)   Add optiion 'no_input'.
!! @li      2006-06-02 (H.Tomita)   Avoid the blocking of communication in prgvar_set_in.
!! @li      2007-02-09 (M.Satoh)    Change default: input_direct_access, output_direct_access
!! @li      2007-06-05 (M.Satoh)    add option: TRC_VMAX_INPUT
!! @li      2008-04-12 (T.Mitsui)   add prgvar_get_noq
!! @li      2008-04-12 (T.Mitsui)   add Nc,Nr,Ni,Ns,Ng restart and diagnosis option
!! @li      2008-05-24 (T.Mitsui)   trivial fix (only for log file)
!! @li      2009-04-14 (T.Mitsui)   add option diag_qi (Re-analysis data implicitly contains ice)
!! @li      2009-08-18 (T.Mitusi)   add option opt_qcqi_to_qv (qc and qi derived from Re-analysis may not be reliable)
!! @li      2011-07-22 (T.Ohno)     add subroutines for plane hgrid systems
!! @li      2011-09-03 (H.Yashiro)  New I/O
!! @li      2011-11-30 (S.Iga)      Avoid irregal isend/irecv buffer array (thanks to T.Inoue)
!! @li      2012-07-23 (H.Yashiro)  [add] Water Isotope by K.Yoshimura
!!
!<
module mod_prgvar
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_MAXFNAME, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: prgvar_setup
  public :: prgvar_get
  public :: prgvar_get_noq
  public :: prgvar_get_withdiag
  public :: prgvar_set
  public :: prgvar_get_in
  public :: prgvar_get_in_withdiag
  public :: prgvar_set_in

  public :: prgvar_get_plane_nopl     ! [add] T.Ohno 110722
  public :: prgvar_get_noq_plane_nopl ! [add] T.Ohno 110722
  public :: prgvar_set_plane_nopl     ! [add] T.Ohno 110722

  public :: restart_input
  public :: restart_output

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=ADM_MAXFNAME), public, save :: restart_input_basename  = ''
  character(len=ADM_MAXFNAME), public, save :: restart_output_basename = ''

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: cnvvar_prg2diag
  private :: cnvvar_diag2prg

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: PRG_vmax0 = 6

  integer, private, parameter :: I_RHOG     =  1 ! Density x G^{1/2}
  integer, private, parameter :: I_RHOGVX   =  2 ! Density x G^{1/2} x Horizontal velocity (X-direction)
  integer, private, parameter :: I_RHOGVY   =  3 ! Density x G^{1/2} x Horizontal velocity (Y-direction)
  integer, private, parameter :: I_RHOGVZ   =  4 ! Density x G^{1/2} x Horizontal velocity (Z-direction)
  integer, private, parameter :: I_RHOGW    =  5 ! Density x G^{1/2} x Vertical   velocity
  integer, private, parameter :: I_RHOGE    =  6 ! Density x G^{1/2} x Energy

  integer, private, save      :: I_RHOGQstr =  7 ! tracers
  integer, private, save      :: I_RHOGQend = -1 !

  character(len=16), private, save :: PRG_name(PRG_vmax0)
  data PRG_name / 'rhog', 'rhogvx', 'rhogvy', 'rhogvz', 'rhogw', 'rhoge' /

  integer, private, parameter :: DIAG_vmax0 = 6

  integer, private, parameter :: I_pre  =  1 ! Pressure
  integer, private, parameter :: I_tem  =  2 ! Temperature
  integer, private, parameter :: I_vx   =  3 ! Horizontal velocity (X-direction)
  integer, private, parameter :: I_vy   =  4 ! Horizontal velocity (Y-direction)
  integer, private, parameter :: I_vz   =  5 ! Horizontal velocity (Z-direction)
  integer, private, parameter :: I_w    =  6 ! Vertical   velocity

  integer, private, save      :: I_qstr =  7 ! tracers
  integer, private, save      :: I_qend = -1 !

  character(len=16), private, save :: DIAG_name(DIAG_vmax0)
  data DIAG_name / 'pre', 'tem', 'vx', 'vy', 'vz', 'w' /

  integer, private,              save :: PRG_vmax  ! total number of prognostic variables
  integer, private,              save :: DIAG_vmax

  real(8), private, allocatable, save :: PRG_var   (:,:,:,:) ! container
  real(8), private, allocatable, save :: PRG_var_pl(:,:,:,:)

  real(8), private, allocatable, save :: PRG_var1   (:,:,:,:) ! container
  real(8), private, allocatable, save :: PRG_var1_pl(:,:,:,:)

  real(8), private, allocatable, save :: DIAG_var   (:,:,:,:) ! container
  real(8), private, allocatable, save :: DIAG_var_pl(:,:,:,:)

  integer, private, save :: TRC_vmax_input ! number of input tracer variables

  character(len=ADM_MAXFNAME), private, save :: layername            = ''       ! [add] H.Yashiro 20110826
  character(len=ADM_MAXFNAME), private, save :: input_io_mode        = 'LEGACY' ! [add] H.Yashiro 20110819
  character(len=ADM_MAXFNAME), private, save :: output_io_mode       = 'LEGACY' ! [add] H.Yashiro 20110819
  logical,                     private, save :: input_direct_access  = .true.
  logical,                     private, save :: output_direct_access = .true.
  logical,                     private, save :: allow_missingq       = .false.  ! [add] H.Yashiro 20110906
  logical,                     private, save :: opt_diag_qi          = .false.
  logical,                     private, save :: opt_qcqi_to_qv       = .false.  ! [Add] 09/08/18 T.Mitsui

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine prgvar_setup
  !>
  subroutine prgvar_setup
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_kall,      &
       ADM_lall,      &
       ADM_lall_pl
    use mod_runconf, only: &
       TRC_vmax,             &
       TRC_name,             &
       opt_2moment_water,    &
       flag_diagnose_number, &
       RAIN_TYPE,            &
       I_QC,                 &
       I_QR,                 &
       I_QI,                 &
       I_QS,                 &
       I_QG,                 &
       I_NC,                 &
       I_NR,                 &
       I_NI,                 &
       I_NS,                 &
       I_NG
    implicit none

    character(len=ADM_MAXFNAME) :: input_basename    = ''
    character(len=ADM_MAXFNAME) :: output_basename   = 'restart'
    character(len=ADM_MAXFNAME) :: restart_layername = ''

    namelist / RESTARTPARAM /  &
         TRC_vmax_input,       &
         input_basename,       &
         output_basename,      &
         restart_layername,    & !--- [add] H.Yashiro 20110826
         input_io_mode,        & !--- [add] H.Yashiro 20110819
         output_io_mode,       & !--- [add] H.Yashiro 20110819
         input_direct_access,  &
         output_direct_access, &
         allow_missingq,       & !--- [add] H.Yashiro 20110906
         opt_diag_qi,          & !--- option to diagnose qi by temperature, 09/04/14 T.Mitsui
         opt_qcqi_to_qv          !--- option to convert qc and qi into qv,  09/08/18 T.Mitsui

    logical, allocatable :: flag_exist(:) ! [Add] 08/04/12 T.Mitsui

    integer :: ierr
    integer :: nq
    !---------------------------------------------------------------------------

    TRC_vmax_input = TRC_vmax

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[prgvar]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=RESTARTPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** RESTARTPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist RESTARTPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist RESTARTPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,RESTARTPARAM)

    restart_input_basename  = input_basename
    restart_output_basename = output_basename
    layername               = restart_layername

    PRG_vmax   = PRG_vmax0  + TRC_vmax
    I_RHOGQend = PRG_vmax

    DIAG_vmax  = DIAG_vmax0 + TRC_vmax
    I_qend     = DIAG_vmax

    allocate( flag_exist(TRC_vmax) )
    flag_exist(:) = .true.

    if ( TRC_vmax_input < TRC_vmax ) then
       flag_exist(TRC_vmax_input+1:TRC_vmax) = .false. ! 08/04/12 T.Mitsui

       ! 08/04/12 T.Mitsui, mp_xxxx should set
       ! default value of Nx corresponding to qx
       if ( opt_2moment_water ) then
          if    ( trim(RAIN_TYPE) == "CLOUD_PARAM" ) then

             ! in CLOUD-PARAM framework, we use QC instead of QI
             if( flag_exist(I_QC) .AND. (.NOT. flag_exist(I_NC)) ) flag_diagnose_number(I_NC) = .true.
             if( flag_exist(I_QC) .AND. (.NOT. flag_exist(I_NI)) ) flag_diagnose_number(I_NI) = .true.

          elseif( trim(RAIN_TYPE) == "WARM" ) then

             ! in warm rain framework, we use QC instead of QI
             if( flag_exist(I_QC) .AND. (.NOT. flag_exist(I_NC)) ) flag_diagnose_number(I_NC) = .true.
             if( flag_exist(I_QR) .AND. (.NOT. flag_exist(I_NR)) ) flag_diagnose_number(I_NR) = .true.
             if( flag_exist(I_QC) .AND. (.NOT. flag_exist(I_NI)) ) flag_diagnose_number(I_NI) = .true.

          elseif( trim(RAIN_TYPE) == "COLD" ) then

             if( flag_exist(I_QC) .AND. (.NOT. flag_exist(I_NC)) ) flag_diagnose_number(I_NC) = .true.
             if( flag_exist(I_QR) .AND. (.NOT. flag_exist(I_NR)) ) flag_diagnose_number(I_NR) = .true.
             if( flag_exist(I_QI) .AND. (.NOT. flag_exist(I_NI)) ) flag_diagnose_number(I_NI) = .true.
             if( flag_exist(I_QS) .AND. (.NOT. flag_exist(I_NS)) ) flag_diagnose_number(I_NS) = .true.
             if( flag_exist(I_QG) .AND. (.NOT. flag_exist(I_NG)) ) flag_diagnose_number(I_NG) = .true.

          endif
       endif
    endif

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '--- Restart treatment for tracers'
    if ( allow_missingq ) then
       write(ADM_LOG_FID,*) '*** Allow missing tracer in restart file.'
       write(ADM_LOG_FID,*) '*** Value will set to zero for missing tracer.'
    endif
    write(ADM_LOG_FID,*) '|========================================'
    write(ADM_LOG_FID,*) '    varname         :read?     :diagnose?'
    do nq = 1, TRC_vmax
       write(ADM_LOG_FID,'(1x,A,A16,A,L1,A,L1)') '--- ', TRC_name(nq), ':',&
                          flag_exist(nq),          '        :', &
                          flag_diagnose_number(nq)
    enddo

    ! -> [add] H.Yashiro 20110819
    write(ADM_LOG_FID,*) '*** io_mode for restart, input : ', trim(input_io_mode)
    if    ( input_io_mode == 'ADVANCED' ) then
    elseif( input_io_mode == 'LEGACY'   ) then
    elseif( input_io_mode == 'IDEAL'    ) then
    else
       write(ADM_LOG_FID,*) 'xxx Invalid input_io_mode. STOP.'
       call ADM_proc_stop
    endif

    write(ADM_LOG_FID,*) '*** io_mode for restart, output: ', trim(output_io_mode)
    if    ( output_io_mode == 'ADVANCED' ) then
    elseif( output_io_mode == 'LEGACY'   ) then
    else
       write(ADM_LOG_FID,*) 'xxx Invalid output_io_mode. STOP'
       call ADM_proc_stop
    endif
    ! <- [add] H.Yashiro 20110819

    allocate( PRG_var    (ADM_gall,   ADM_kall,ADM_lall,   PRG_vmax) )
    allocate( PRG_var_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax) )

    allocate( PRG_var1   (ADM_gall,   ADM_kall,ADM_lall,   PRG_vmax) )
    allocate( PRG_var1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax) )

    allocate( DIAG_var   (ADM_gall,   ADM_kall,ADM_lall,   DIAG_vmax) )
    allocate( DIAG_var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,DIAG_vmax) )

    return
  end subroutine prgvar_setup

  !-----------------------------------------------------------------------------
  !>
  !> get prognostic variables from prg[num].
  !>
  subroutine prgvar_get( &
       rhog,   rhog_pl,   &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       rhoge,  rhoge_pl,  &
       rhogq,  rhogq_pl,  &
       num                )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(out) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogq    (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(8), intent(out) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    integer, intent(in)  :: num

    integer :: n, k, l, nq
    !---------------------------------------------------------------------------

    if ( num == 0 ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          rhog  (n,k,l) = PRG_var(n,k,l,I_RHOG  )
          rhogvx(n,k,l) = PRG_var(n,k,l,I_RHOGVX)
          rhogvy(n,k,l) = PRG_var(n,k,l,I_RHOGVY)
          rhogvz(n,k,l) = PRG_var(n,k,l,I_RHOGVZ)
          rhogw (n,k,l) = PRG_var(n,k,l,I_RHOGW )
          rhoge (n,k,l) = PRG_var(n,k,l,I_RHOGE )
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall
          rhogq(n,k,l,nq) = PRG_var(n,k,l,PRG_vmax0+nq)
       enddo
       enddo
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then

          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl
             rhog_pl  (n,k,l) = PRG_var_pl(n,k,l,I_RHOG  )
             rhogvx_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVX)
             rhogvy_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVY)
             rhogvz_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVZ)
             rhogw_pl (n,k,l) = PRG_var_pl(n,k,l,I_RHOGW )
             rhoge_pl (n,k,l) = PRG_var_pl(n,k,l,I_RHOGE )
          enddo
          enddo
          enddo

          do nq = 1, TRC_vmax
          do l  = 1, ADM_lall_pl
          do k  = 1, ADM_kall
          do n  = 1, ADM_gall_pl
             rhogq_pl(n,k,l,nq) = PRG_var_pl(n,k,l,PRG_vmax0+nq)
          enddo
          enddo
          enddo
          enddo

       endif

    elseif( num == 1 ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          rhog  (n,k,l) = PRG_var1(n,k,l,I_RHOG  )
          rhogvx(n,k,l) = PRG_var1(n,k,l,I_RHOGVX)
          rhogvy(n,k,l) = PRG_var1(n,k,l,I_RHOGVY)
          rhogvz(n,k,l) = PRG_var1(n,k,l,I_RHOGVZ)
          rhogw (n,k,l) = PRG_var1(n,k,l,I_RHOGW )
          rhoge (n,k,l) = PRG_var1(n,k,l,I_RHOGE )
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall
          rhogq(n,k,l,nq) = PRG_var1(n,k,l,PRG_vmax0+nq)
       enddo
       enddo
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then

          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl
             rhog_pl  (n,k,l) = PRG_var1_pl(n,k,l,I_RHOG  )
             rhogvx_pl(n,k,l) = PRG_var1_pl(n,k,l,I_RHOGVX)
             rhogvy_pl(n,k,l) = PRG_var1_pl(n,k,l,I_RHOGVY)
             rhogvz_pl(n,k,l) = PRG_var1_pl(n,k,l,I_RHOGVZ)
             rhogw_pl (n,k,l) = PRG_var1_pl(n,k,l,I_RHOGW )
             rhoge_pl (n,k,l) = PRG_var1_pl(n,k,l,I_RHOGE )
          enddo
          enddo
          enddo

          do nq = 1, TRC_vmax
          do l  = 1, ADM_lall_pl
          do k  = 1, ADM_kall
          do n  = 1, ADM_gall_pl
             rhogq_pl(n,k,l,nq) = PRG_var1_pl(n,k,l,PRG_vmax0+nq)
          enddo
          enddo
          enddo
          enddo

       endif

    endif

    return
  end subroutine prgvar_get

  !----------------------------------------------------------------------------
  !>
  !> 08/04/12 [Add] T.Mitsui, only tracer Q is not used here
  !>
  subroutine prgvar_get_noq( &
       rhog,   rhog_pl,   &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       rhoge,  rhoge_pl   )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    implicit none

    real(8), intent(out) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: n, k, l
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_gall
       rhog  (n,k,l) = PRG_var(n,k,l,I_RHOG  )
       rhogvx(n,k,l) = PRG_var(n,k,l,I_RHOGVX)
       rhogvy(n,k,l) = PRG_var(n,k,l,I_RHOGVY)
       rhogvz(n,k,l) = PRG_var(n,k,l,I_RHOGVZ)
       rhogw (n,k,l) = PRG_var(n,k,l,I_RHOGW )
       rhoge (n,k,l) = PRG_var(n,k,l,I_RHOGE )
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do n = 1, ADM_gall_pl
          rhog_pl  (n,k,l) = PRG_var_pl(n,k,l,I_RHOG  )
          rhogvx_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVX)
          rhogvy_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVY)
          rhogvz_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVZ)
          rhogw_pl (n,k,l) = PRG_var_pl(n,k,l,I_RHOGW )
          rhoge_pl (n,k,l) = PRG_var_pl(n,k,l,I_RHOGE )
       enddo
       enddo
       enddo

    endif

    return
  end subroutine prgvar_get_noq

  !-----------------------------------------------------------------------------
  !>
  !> get prognostic variables from prg[num].
  !>
  subroutine prgvar_get_withdiag( &
       rhog,   rhog_pl,   &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       rhoge,  rhoge_pl,  &
       rhogq,  rhogq_pl,  &
       rho,    rho_pl,    &
       pre,    pre_pl,    &
       tem,    tem_pl,    &
       vx,     vx_pl,     &
       vy,     vy_pl,     &
       vz,     vz_pl,     &
       w,      w_pl,      &
       q,      q_pl       )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_vmtr, only: &
       VMTR_GSGAM2,    &
       VMTR_GSGAM2_pl
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(out) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogq    (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(8), intent(out) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(8), intent(out) :: rho      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: pre      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: tem      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: vx       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: vy       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: vz       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: w        (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: q        (ADM_gall,   ADM_kall,ADM_lall,   TRC_vmax)
    real(8), intent(out) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    integer :: n, k, l, nq
    !---------------------------------------------------------------------------

    call cnvvar_prg2diag( PRG_var (:,:,:,:), PRG_var_pl (:,:,:,:), & !--- [IN]
                          DIAG_var(:,:,:,:), DIAG_var_pl(:,:,:,:)  ) !--- [OUT]

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_gall
       rhog  (n,k,l) = PRG_var(n,k,l,I_RHOG  )
       rhogvx(n,k,l) = PRG_var(n,k,l,I_RHOGVX)
       rhogvy(n,k,l) = PRG_var(n,k,l,I_RHOGVY)
       rhogvz(n,k,l) = PRG_var(n,k,l,I_RHOGVZ)
       rhogw (n,k,l) = PRG_var(n,k,l,I_RHOGW )
       rhoge (n,k,l) = PRG_var(n,k,l,I_RHOGE )

       rho   (n,k,l) = PRG_var(n,k,l,I_RHOG) / VMTR_GSGAM2(n,k,l)

       pre   (n,k,l) = DIAG_var(n,k,l,I_pre)
       tem   (n,k,l) = DIAG_var(n,k,l,I_tem)
       vx    (n,k,l) = DIAG_var(n,k,l,I_vx )
       vy    (n,k,l) = DIAG_var(n,k,l,I_vy )
       vz    (n,k,l) = DIAG_var(n,k,l,I_vz )
       w     (n,k,l) = DIAG_var(n,k,l,I_w  )
    enddo
    enddo
    enddo

    do nq = 1, TRC_vmax
    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do n  = 1, ADM_gall
       rhogq(n,k,l,nq) = PRG_var(n,k,l,PRG_vmax0+nq)

       q    (n,k,l,nq) = DIAG_var(n,k,l,DIAG_vmax0+nq)
    enddo
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do n = 1, ADM_gall_pl
          rhog_pl  (n,k,l) = PRG_var_pl(n,k,l,I_RHOG  )
          rhogvx_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVX)
          rhogvy_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVY)
          rhogvz_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVZ)
          rhogw_pl (n,k,l) = PRG_var_pl(n,k,l,I_RHOGW )
          rhoge_pl (n,k,l) = PRG_var_pl(n,k,l,I_RHOGE )

          rho_pl   (n,k,l) = PRG_var_pl(n,k,l,I_RHOG) / VMTR_GSGAM2_pl(n,k,l)

          pre_pl   (n,k,l) = DIAG_var_pl(n,k,l,I_pre)
          tem_pl   (n,k,l) = DIAG_var_pl(n,k,l,I_tem)
          vx_pl    (n,k,l) = DIAG_var_pl(n,k,l,I_vx )
          vy_pl    (n,k,l) = DIAG_var_pl(n,k,l,I_vy )
          vz_pl    (n,k,l) = DIAG_var_pl(n,k,l,I_vz )
          w_pl     (n,k,l) = DIAG_var_pl(n,k,l,I_w  )
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall_pl
          rhogq_pl(n,k,l,nq) = PRG_var_pl(n,k,l,PRG_vmax0+nq)

          q_pl    (n,k,l,nq) = DIAG_var_pl(n,k,l,DIAG_vmax0+nq)
       enddo
       enddo
       enddo
       enddo

    endif

    return
  end subroutine prgvar_get_withdiag

  !-----------------------------------------------------------------------------
  !>
  !> set prognostic variables to prg[num]- and COMMUNICATION.
  !>
  subroutine prgvar_set( &
       rhog,   rhog_pl,   &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       rhoge,  rhoge_pl,  &
       rhogq,  rhogq_pl,  &
       num                )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax
    use mod_comm, only: &
       COMM_data_transfer
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(in)  :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogq    (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(8), intent(in)  :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    integer, intent(in)  :: num

    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    integer :: n, k, l, nq
    !---------------------------------------------------------------------------

    if ( num == 0 ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          PRG_var(n,k,l,I_RHOG  ) = rhog  (n,k,l)
          PRG_var(n,k,l,I_RHOGVX) = rhogvx(n,k,l)
          PRG_var(n,k,l,I_RHOGVY) = rhogvy(n,k,l)
          PRG_var(n,k,l,I_RHOGVZ) = rhogvz(n,k,l)
          PRG_var(n,k,l,I_RHOGW ) = rhogw (n,k,l)
          PRG_var(n,k,l,I_RHOGE ) = rhoge (n,k,l)
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall
          PRG_var(n,k,l,PRG_vmax0+nq) = rhogq(n,k,l,nq)
       enddo
       enddo
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then

          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl
             PRG_var_pl(n,k,l,I_RHOG  ) = rhog_pl  (n,k,l)
             PRG_var_pl(n,k,l,I_RHOGVX) = rhogvx_pl(n,k,l)
             PRG_var_pl(n,k,l,I_RHOGVY) = rhogvy_pl(n,k,l)
             PRG_var_pl(n,k,l,I_RHOGVZ) = rhogvz_pl(n,k,l)
             PRG_var_pl(n,k,l,I_RHOGW ) = rhogw_pl (n,k,l)
             PRG_var_pl(n,k,l,I_RHOGE ) = rhoge_pl (n,k,l)
          enddo
          enddo
          enddo

          do nq = 1, TRC_vmax
          do l  = 1, ADM_lall_pl
          do k  = 1, ADM_kall
          do n  = 1, ADM_gall_pl
             PRG_var_pl(n,k,l,PRG_vmax0+nq) = rhogq_pl(n,k,l,nq)
          enddo
          enddo
          enddo
          enddo

       endif

       ! communication
       call COMM_data_transfer( PRG_var, PRG_var_pl )

       PRG_var(suf(ADM_gall_1d,1),:,:,:) = PRG_var(suf(ADM_gmax+1,ADM_gmin),:,:,:)
       PRG_var(suf(1,ADM_gall_1d),:,:,:) = PRG_var(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    elseif( num == 1 ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          PRG_var1(n,k,l,I_RHOG  ) = rhog  (n,k,l)
          PRG_var1(n,k,l,I_RHOGVX) = rhogvx(n,k,l)
          PRG_var1(n,k,l,I_RHOGVY) = rhogvy(n,k,l)
          PRG_var1(n,k,l,I_RHOGVZ) = rhogvz(n,k,l)
          PRG_var1(n,k,l,I_RHOGW ) = rhogw (n,k,l)
          PRG_var1(n,k,l,I_RHOGE ) = rhoge (n,k,l)
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall
          PRG_var1(n,k,l,PRG_vmax0+nq) = rhogq(n,k,l,nq)
       enddo
       enddo
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then

          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl
             PRG_var1_pl(n,k,l,I_RHOG  ) = rhog_pl  (n,k,l)
             PRG_var1_pl(n,k,l,I_RHOGVX) = rhogvx_pl(n,k,l)
             PRG_var1_pl(n,k,l,I_RHOGVY) = rhogvy_pl(n,k,l)
             PRG_var1_pl(n,k,l,I_RHOGVZ) = rhogvz_pl(n,k,l)
             PRG_var1_pl(n,k,l,I_RHOGW ) = rhogw_pl (n,k,l)
             PRG_var1_pl(n,k,l,I_RHOGE ) = rhoge_pl (n,k,l)
          enddo
          enddo
          enddo

          do nq = 1, TRC_vmax
          do l  = 1, ADM_lall_pl
          do k  = 1, ADM_kall
          do n  = 1, ADM_gall_pl
             PRG_var1_pl(n,k,l,PRG_vmax0+nq) = rhogq_pl(n,k,l,nq)
          enddo
          enddo
          enddo
          enddo

       endif

       ! communication
       call COMM_data_transfer( PRG_var1, PRG_var1_pl )

       PRG_var1(suf(ADM_gall_1d,1),:,:,:) = PRG_var1(suf(ADM_gmax+1,ADM_gmin),:,:,:)
       PRG_var1(suf(1,ADM_gall_1d),:,:,:) = PRG_var1(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    endif

    return
  end subroutine prgvar_set

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine prgvar_get_in
  !>
  subroutine prgvar_get_in( &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogvz, &
       rhogw,  &
       rhoge,  &
       rhogq   )
    use mod_adm, only: &
       ADM_gall,        &
       ADM_kall,        &
       ADM_lall,        &
       ADM_IopJop_nmax, &
       ADM_IopJop,      &
       ADM_GIoJo
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(out) :: rhog  (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvx(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvy(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvz(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogw (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhoge (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogq (ADM_IopJop_nmax,ADM_kall,ADM_lall,TRC_vmax)

    integer :: n, k, l, nq, ij
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_IopJop_nmax
       ij = ADM_IopJop(n,ADM_GIoJo)

       rhog  (n,k,l) = PRG_var(ij,k,l,I_RHOG  )
       rhogvx(n,k,l) = PRG_var(ij,k,l,I_RHOGVX)
       rhogvy(n,k,l) = PRG_var(ij,k,l,I_RHOGVY)
       rhogvz(n,k,l) = PRG_var(ij,k,l,I_RHOGVZ)
       rhogw (n,k,l) = PRG_var(ij,k,l,I_RHOGW )
       rhoge (n,k,l) = PRG_var(ij,k,l,I_RHOGE )
    enddo
    enddo
    enddo

    do nq = 1, TRC_vmax
    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do n  = 1, ADM_IopJop_nmax
       ij = ADM_IopJop(n,ADM_GIoJo)

       rhogq(n,k,l,nq) = PRG_var(ij,k,l,PRG_vmax0+nq)
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine prgvar_get_in

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine prgvar_get_in
  !>
  subroutine prgvar_get_in_withdiag( &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogvz, &
       rhogw,  &
       rhoge,  &
       rhogq,  &
       rho,    &
       pre,    &
       tem,    &
       vx,     &
       vy,     &
       vz,     &
       w,      &
       q       )
    use mod_adm, only: &
       ADM_gall,        &
       ADM_kall,        &
       ADM_lall,        &
       ADM_IopJop_nmax, &
       ADM_IopJop,      &
       ADM_GIoJo
    use mod_vmtr, only: &
       VMTR_GSGAM2
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(out) :: rhog  (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvx(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvy(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvz(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogw (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhoge (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogq (ADM_IopJop_nmax,ADM_kall,ADM_lall,TRC_vmax)
    real(8), intent(out) :: rho   (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: pre   (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: tem   (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: vx    (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: vy    (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: vz    (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: w     (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(out) :: q     (ADM_IopJop_nmax,ADM_kall,ADM_lall,TRC_vmax)

    integer :: n, k, l, nq, ij
    !---------------------------------------------------------------------------

    call cnvvar_prg2diag( PRG_var (:,:,:,:), PRG_var_pl (:,:,:,:), & !--- [IN]
                          DIAG_var(:,:,:,:), DIAG_var_pl(:,:,:,:)  ) !--- [OUT]

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_IopJop_nmax
       ij = ADM_IopJop(n,ADM_GIoJo)

       rhog  (n,k,l) = PRG_var(ij,k,l,I_RHOG  )
       rhogvx(n,k,l) = PRG_var(ij,k,l,I_RHOGVX)
       rhogvy(n,k,l) = PRG_var(ij,k,l,I_RHOGVY)
       rhogvz(n,k,l) = PRG_var(ij,k,l,I_RHOGVZ)
       rhogw (n,k,l) = PRG_var(ij,k,l,I_RHOGW )
       rhoge (n,k,l) = PRG_var(ij,k,l,I_RHOGE )

       rho   (n,k,l) = PRG_var(ij,k,l,I_RHOG) / VMTR_GSGAM2(ij,k,l)

       pre   (n,k,l) = DIAG_var(ij,k,l,I_pre)
       tem   (n,k,l) = DIAG_var(ij,k,l,I_tem)
       vx    (n,k,l) = DIAG_var(ij,k,l,I_vx )
       vy    (n,k,l) = DIAG_var(ij,k,l,I_vy )
       vz    (n,k,l) = DIAG_var(ij,k,l,I_vz )
       w     (n,k,l) = DIAG_var(ij,k,l,I_w  )
    enddo
    enddo
    enddo

    do nq = 1, TRC_vmax
    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do n  = 1, ADM_IopJop_nmax
       ij = ADM_IopJop(n,ADM_GIoJo)

       rhogq(n,k,l,nq) = PRG_var(ij,k,l,PRG_vmax0+nq)

       q    (n,k,l,nq) = DIAG_var(ij,k,l,DIAG_vmax0+nq)
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine prgvar_get_in_withdiag

  !-----------------------------------------------------------------------------
  !>
  !> set prognostic variables to prg[num] and COMMUNICATION.
  !>
  subroutine prgvar_set_in( &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogvz, &
       rhogw,  &
       rhoge,  &
       rhogq   )
    use mod_adm, only : &
       ADM_gall,        &
       ADM_kall,        &
       ADM_lall,        &
       ADM_IopJop_nmax, &
       ADM_IopJop,      &
       ADM_GIoJo,       & 
       ADM_gall_1d,     &
       ADM_gmax,        &
       ADM_gmin
    use mod_comm, only: &
       COMM_data_transfer,       &
       COMM_data_transfer_rgn2pl
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(in) :: rhog  (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvx(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvy(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvz(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogw (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhoge (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogq (ADM_IopJop_nmax,ADM_kall,ADM_lall,TRC_vmax)


    integer :: n, k, l, nq, ij

    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_IopJop_nmax
       ij = ADM_IopJop(n,ADM_GIoJo)

       PRG_var(ij,k,l,I_RHOG  ) = rhog  (n,k,l)
       PRG_var(ij,k,l,I_RHOGVX) = rhogvx(n,k,l)
       PRG_var(ij,k,l,I_RHOGVY) = rhogvy(n,k,l)
       PRG_var(ij,k,l,I_RHOGVZ) = rhogvz(n,k,l)
       PRG_var(ij,k,l,I_RHOGW ) = rhogw (n,k,l)
       PRG_var(ij,k,l,I_RHOGE ) = rhoge (n,k,l)
    enddo
    enddo
    enddo

    do nq = 1, TRC_vmax
    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do n = 1, ADM_IopJop_nmax
       ij = ADM_IopJop(n,ADM_GIoJo)

       PRG_var(ij,k,l,PRG_vmax0+nq) = rhogq(n,k,l,nq)
    enddo
    enddo
    enddo
    enddo

    ! communication
    call COMM_data_transfer_rgn2pl( PRG_var, PRG_var_pl, ADM_kall, PRG_vmax )

    call COMM_data_transfer( PRG_var, PRG_var_pl )

    PRG_var(suf(ADM_gall_1d,1),:,:,:) = PRG_var(suf(ADM_gmax+1,ADM_gmin),:,:,:)
    PRG_var(suf(1,ADM_gall_1d),:,:,:) = PRG_var(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    return
  end subroutine prgvar_set_in

  !-----------------------------------------------------------------------------
  !>
  !> get prognostic variables from prg[num].
  !>
  subroutine prgvar_get_plane_nopl( &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogw,  &
       rhoge,  &
       rhogq,  &
       num     )
    use mod_adm, only: &
       ADM_gall, &
       ADM_kall, &
       ADM_lall
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(out) :: rhog  (ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogw (ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhoge (ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogq (ADM_gall,ADM_kall,ADM_lall,TRC_vmax)
    integer, intent(in)  :: num

    integer :: n, k, l, nq
    !---------------------------------------------------------------------------

    if ( num == 0 ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          rhog  (n,k,l) = PRG_var(n,k,l,I_RHOG  )
          rhogvx(n,k,l) = PRG_var(n,k,l,I_RHOGVX)
          rhogvy(n,k,l) = PRG_var(n,k,l,I_RHOGVY)
          rhogw (n,k,l) = PRG_var(n,k,l,I_RHOGW )
          rhoge (n,k,l) = PRG_var(n,k,l,I_RHOGE )
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall
          rhogq(n,k,l,nq) = PRG_var(n,k,l,PRG_vmax0+nq)
       enddo
       enddo
       enddo
       enddo

    elseif( num == 1 ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          rhog  (n,k,l) = PRG_var1(n,k,l,I_RHOG  )
          rhogvx(n,k,l) = PRG_var1(n,k,l,I_RHOGVX)
          rhogvy(n,k,l) = PRG_var1(n,k,l,I_RHOGVY)
          rhogw (n,k,l) = PRG_var1(n,k,l,I_RHOGW )
          rhoge (n,k,l) = PRG_var1(n,k,l,I_RHOGE )
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall
          rhogq(n,k,l,nq) = PRG_var1(n,k,l,PRG_vmax0+nq)
       enddo
       enddo
       enddo
       enddo

    endif

    return
  end subroutine prgvar_get_plane_nopl

  !-----------------------------------------------------------------------------
  !>
  !> get prognostic variables from prg[num].
  !>
  subroutine prgvar_get_noq_plane_nopl( &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogw,  &
       rhoge   )
    use mod_adm, only: &
       ADM_gall, &
       ADM_kall, &
       ADM_lall
    implicit none

    real(8), intent(out) :: rhog  (ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhogw (ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: rhoge (ADM_gall,ADM_kall,ADM_lall)

    integer :: n, k, l
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_gall
       rhog  (n,k,l) = PRG_var(n,k,l,I_RHOG  )
       rhogvx(n,k,l) = PRG_var(n,k,l,I_RHOGVX)
       rhogvy(n,k,l) = PRG_var(n,k,l,I_RHOGVY)
       rhogw (n,k,l) = PRG_var(n,k,l,I_RHOGW )
       rhoge (n,k,l) = PRG_var(n,k,l,I_RHOGE )
    enddo
    enddo
    enddo

    return
  end subroutine prgvar_get_noq_plane_nopl

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine prgvar_set_plane_nopl
  !>
  subroutine prgvar_set_plane_nopl( &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogw,  &
       rhoge,  &
       rhogq,  &
       num     )
    use mod_adm, only: &
       ADM_gall,    &
       ADM_lall,    &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax
    use mod_comm, only: &
       COMM_data_transfer
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(8), intent(in) :: rhog  (ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogw (ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhoge (ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogq (ADM_gall,ADM_kall,ADM_lall,TRC_vmax)
    integer, intent(in) :: num

    integer :: n, k, l, nq

    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    if ( num == 0 ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          PRG_var(n,k,l,I_RHOG  ) = rhog  (n,k,l)
          PRG_var(n,k,l,I_RHOGVX) = rhogvx(n,k,l)
          PRG_var(n,k,l,I_RHOGVY) = rhogvy(n,k,l)
          PRG_var(n,k,l,I_RHOGW ) = rhogw (n,k,l)
          PRG_var(n,k,l,I_RHOGE ) = rhoge (n,k,l)
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall
          PRG_var(n,k,l,PRG_vmax0+nq) = rhogq(n,k,l,nq)
       enddo
       enddo
       enddo
       enddo

       call COMM_data_transfer( PRG_var, PRG_var_pl )

       PRG_var(suf(ADM_gall_1d,1),:,:,:) = PRG_var(suf(ADM_gmax+1,ADM_gmin),:,:,:)
       PRG_var(suf(1,ADM_gall_1d),:,:,:) = PRG_var(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    elseif( num == 1 ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          PRG_var1(n,k,l,I_RHOG  ) = rhog  (n,k,l)
          PRG_var1(n,k,l,I_RHOGVX) = rhogvx(n,k,l)
          PRG_var1(n,k,l,I_RHOGVY) = rhogvy(n,k,l)
          PRG_var1(n,k,l,I_RHOGW ) = rhogw (n,k,l)
          PRG_var1(n,k,l,I_RHOGE ) = rhoge (n,k,l)
       enddo
       enddo
       enddo

       do nq = 1, TRC_vmax
       do l  = 1, ADM_lall
       do k  = 1, ADM_kall
       do n  = 1, ADM_gall
          PRG_var1(n,k,l,PRG_vmax0+nq) = rhogq(n,k,l,nq)
       enddo
       enddo
       enddo
       enddo

       call COMM_data_transfer( PRG_var1, PRG_var1_pl )

       PRG_var1(suf(ADM_gall_1d,1),:,:,:) = PRG_var1(suf(ADM_gmax+1,ADM_gmin),:,:,:)
       PRG_var1(suf(1,ADM_gall_1d),:,:,:) = PRG_var1(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    endif

    return
  end subroutine prgvar_set_plane_nopl

  !-----------------------------------------------------------------------------
  subroutine restart_input( basename )
    use mod_misc, only: &
       MISC_get_available_fid, &
       MISC_make_idstr
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_lall,    &
       ADM_kall,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_fio, only: &
       FIO_input
    use mod_comm, only: &
       COMM_var
    use mod_gtl, only: &
       GTL_max, &
       GTL_min
    use mod_runconf, only: &
       opt_2moment_water,    &
       flag_diagnose_number, &
       TRC_vmax, &
       TRC_name, &
       I_QC,     &
       I_QI,     &
       I_NI!,     &
!       ISOTOPE,  & ! [add] K.Yoshimura 20120414
!       ISO_MAX,  & ! [add] K.Yoshimura 20120414
!       ISO_STR,  & ! [add] K.Yoshimura 20120414
!       ISO_STR2    ! [add] K.Yoshimura 20120414
    use mod_dycoretest, only :  & ! [add] R.Yoshida 20121019
       dycore_input
    implicit none

    character(len=ADM_MAXFNAME), intent(in) :: basename

    character(len=ADM_MAXFNAME) :: fname

    real(8) :: val_max, val_min
    logical :: nonzero

    integer :: fid
    integer :: l, rgnid, nq
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '*** read restart/initial data'

    ! -> [add]&[Mod] H.Yashiro 20110819
    if ( input_io_mode == 'ADVANCED' ) then

       do nq = 1, DIAG_vmax0
          call FIO_input(DIAG_var(:,:,:,nq),basename,DIAG_name(nq), & 
                         layername,1,ADM_kall,1                     )
       enddo

       do nq = 1, TRC_vmax_input
          call FIO_input(DIAG_var(:,:,:,DIAG_vmax0+nq),basename,TRC_name(nq), &
                         layername,1,ADM_kall,1,                              &
                         allow_missingq=allow_missingq                        )
       enddo

       call comm_var( DIAG_var, DIAG_var_pl, ADM_kall, DIAG_vmax, comm_type=2, NSval_fix=.true. )

    elseif( input_io_mode == 'LEGACY' ) then

       if ( input_direct_access ) then !--- direct access ( defalut )

          do l = 1, ADM_lall
             rgnid = ADM_prc_tab(l,ADM_prc_me)
             call MISC_make_idstr(fname,trim(basename),'rgn',rgnid)
             fid = MISC_get_available_fid()
             open( unit   = fid,                 &
                   file   = trim(fname),         &
                   form   = 'unformatted',       &
                   access = 'direct',            &
                   recl   = ADM_gall*ADM_kall*8, &
                   status = 'old'                )

             do nq = 1, DIAG_vmax0+TRC_vmax_input
                read(fid,rec=nq) DIAG_var(:,:,l,nq)
             enddo

             close(fid)
          enddo

          call comm_var( DIAG_var, DIAG_var_pl, ADM_kall, DIAG_vmax, comm_type=2, NSval_fix=.true. )

       else !--- sequential access

          do l = 1, ADM_lall
             rgnid = ADM_prc_tab(l,ADM_prc_me)
             call MISC_make_idstr(fname,trim(basename),'rgn',rgnid)
             fid = MISC_get_available_fid()
             open( unit   = fid,                 &
                   file   = trim(fname),         &
                   form   = 'unformatted',       &
                   access = 'sequential',        &
                   status = 'old'                )

             do nq = 1, DIAG_vmax0+TRC_vmax_input
                read(fid) DIAG_var(:,:,l,nq)
             enddo

             close(fid)
          enddo

          if ( ADM_prc_me == ADM_prc_pl ) then
             fname = trim(basename)//'.pl'
             fid = MISC_get_available_fid()
             open( unit   = fid,                 &
                   file   = trim(fname),         &
                   form   = 'unformatted',       &
                   access = 'sequential',        &
                   status = 'old'                )

             do nq = 1, DIAG_vmax0+TRC_vmax_input
                read(fid) DIAG_var_pl(:,:,:,nq)
             enddo

             close(fid)
          endif
       endif !--- direct/sequencial
    elseif( input_io_mode == 'IDEAL' ) then

       write(ADM_LOG_FID,*) '*** make ideal initials'

       call dycore_input( DIAG_var(:,:,:,:) ) !--- [OUT]

       call comm_var( DIAG_var, DIAG_var_pl, ADM_kall, DIAG_vmax, comm_type=2, NSval_fix=.true. )

    endif !--- io_mode
    ! <- [add]&[Mod] H.Yashiro 20110819

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '====== data range check : diagnostic variables ======'
    do nq = 1, DIAG_vmax0
       val_max = GTL_max( DIAG_var   (:,:,:,nq),       &
                          DIAG_var_pl(:,:,:,nq),       &
                          ADM_kall, ADM_kmin, ADM_kmax )
       val_min = GTL_min( DIAG_var   (:,:,:,nq),       &
                          DIAG_var_pl(:,:,:,nq),       &
                          ADM_kall, ADM_kmin, ADM_kmax )

       write(ADM_LOG_FID,'(1x,A,A16,2(A,1PE24.17))') '--- ', DIAG_name(nq), ': max=', val_max, ', min=', val_min
    enddo

    do nq = 1, TRC_vmax
       val_max = GTL_max( DIAG_var   (:,:,:,DIAG_vmax0+nq), &
                          DIAG_var_pl(:,:,:,DIAG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax      )

       if ( val_max <= 0.D0 ) then
          nonzero = .false.
       else
          nonzero = .true.
       endif

       val_min = GTL_min( DIAG_var   (:,:,:,DIAG_vmax0+nq), &
                          DIAG_var_pl(:,:,:,DIAG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax, nonzero)

       write(ADM_LOG_FID,'(1x,A,A16,2(A,1PE24.17))') '--- ', TRC_name(nq),  ': max=', val_max, ', min=', val_min
    enddo

!    if ( trim(ISOTOPE) == "ON" ) then ! [add] K.Yoshimura 20110414
!       do nq = 1, ISO_MAX/2
!          DIAG_var(:,:,:,ISO_STR +5+nq) = DIAG_var(:,:,:,6+nq) * 0.99D0
!          DIAG_var(:,:,:,ISO_STR2+5+nq) = DIAG_var(:,:,:,6+nq) * 0.92D0
!       enddo
!    endif

    call cnvvar_diag2prg( PRG_var (:,:,:,:), PRG_var_pl (:,:,:,:), & !--- [OUT]
                          DIAG_var(:,:,:,:), DIAG_var_pl(:,:,:,:)  ) !--- [IN]


    ! [Add] 09/04/14 T.Mitsui
    ! "option is chosen" and "QI is tracer"
    if ( I_QI <= 0 ) then
       if ( opt_diag_qi ) then
          call diag_qi( ADM_gall,              & !--- [IN]
                        ADM_kall,              & !--- [IN]
                        ADM_lall,              & !--- [IN]
                        DIAG_var(:,:,:,I_tem), & !--- [IN]
                        DIAG_var(:,:,:,I_QC),  & !--- [INOUT]
                        DIAG_var(:,:,:,I_QI),  & !--- [INOUT]
                        PRG_var (:,:,:,I_QC),  & !--- [INOUT]
                        PRG_var (:,:,:,I_QI)   ) !--- [INOUT]
          if ( opt_2moment_water ) then
             if ( I_NI <= 0 ) then
                flag_diagnose_number(I_NI) = .true.
             endif
          endif
       endif
    endif

    ![Add] 09/08/18 T.Mitsui
    ! qc and qi are converted into qv and temperature is modified.
    if ( opt_qcqi_to_qv )then
       call convert_qcqi_to_qv( ADM_gall,                             & !--- [IN]
                                ADM_kall,                             & !--- [IN]
                                ADM_lall,                             & !--- [IN]
                                TRC_vmax,                             & !--- [IN]
                                DIAG_var(:,:,:,I_tem),                & !--- [INOUT]
                                DIAG_var(:,:,:,I_qstr:I_qend),        & !--- [INOUT]
                                PRG_var (:,:,:,I_RHOGQstr:I_RHOGQend) ) !--- [INOUT]
    endif

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '====== data range check : prognostic variables ======'
    do nq = 1, DIAG_vmax0
       val_max = GTL_max( PRG_var   (:,:,:,nq),        &
                          PRG_var_pl(:,:,:,nq),        &
                          ADM_kall, ADM_kmin, ADM_kmax )
       val_min = GTL_min( PRG_var   (:,:,:,nq),        &
                          PRG_var_pl(:,:,:,nq),        &
                          ADM_kall, ADM_kmin, ADM_kmax )

       write(ADM_LOG_FID,'(1x,A,A16,2(A,1PE24.17))') '--- ', PRG_name(nq), ': max=', val_max, ', min=', val_min
    enddo

    do nq = 1, TRC_vmax
       val_max = GTL_max( PRG_var   (:,:,:,PRG_vmax0+nq), &
                          PRG_var_pl(:,:,:,PRG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax    )

       if ( val_max <= 0.D0 ) then
          nonzero = .false.
       else
          nonzero = .true.
       endif

       val_min = GTL_min( PRG_var   (:,:,:,PRG_vmax0+nq),      &
                          PRG_var_pl(:,:,:,PRG_vmax0+nq),      &
                          ADM_kall, ADM_kmin, ADM_kmax, nonzero)

       write(ADM_LOG_FID,'(1x,A,A16,2(A,1PE24.17))') '--- rhog * ', TRC_name(nq),  ': max=', val_max, ', min=', val_min
    enddo

    return
  end subroutine restart_input

  !-----------------------------------------------------------------------------
  subroutine restart_output( basename )
    use mod_misc, only: &
       MISC_get_available_fid, &
       MISC_make_idstr
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_lall,    &
       ADM_kall,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_fio, only : & ! [add] H.Yashiro 20110819
       FIO_output, &
       FIO_HSHORT, &
       FIO_HMID,   &
       FIO_REAL8
    use mod_comm, only: &
       COMM_var
    use mod_time, only : &
       TIME_CTIME
    use mod_gtl, only: &
       GTL_max, &
       GTL_min
    use mod_runconf, only : &
       TRC_vmax, &
       TRC_name, &
       WLABEL
    implicit none

    character(len=ADM_MAXFNAME), intent(in) :: basename

    character(len=FIO_HMID)   :: desc = 'INITIAL/RESTART_data_of_prognostic_variables' ! [add] H.Yashiro 20110819

    character(len=FIO_HSHORT) :: DLABEL(DIAG_vmax0)
    data DLABEL / 'Pressure ',        &
                  'Temperature ',     &
                  'H-Velocity(XDIR)', &
                  'H-Velocity(YDIR)', &
                  'H-Velocity(ZDIR)', &
                  'V-Velocity '       /

    character(len=FIO_HSHORT) :: DUNIT(DIAG_vmax0)
    data DUNIT /  'Pa',  &
                  'K',   &
                  'm/s', &
                  'm/s', &
                  'm/s', &
                  'm/s'  /

    character(len=FIO_HSHORT) :: WUNIT = 'kg/kg'

    character(len=ADM_MAXFNAME) :: fname

    real(8) :: val_max, val_min
    logical :: nonzero

    integer :: fid
    integer :: l, rgnid, nq
    !---------------------------------------------------------------------------

    call cnvvar_prg2diag( PRG_var (:,:,:,:), PRG_var_pl (:,:,:,:), & !--- [IN]
                          DIAG_var(:,:,:,:), DIAG_var_pl(:,:,:,:)  ) !--- [OUT]

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '====== data range check : prognostic variables ======'
    do nq = 1, DIAG_vmax0
       val_max = GTL_max( DIAG_var   (:,:,:,nq),       &
                          DIAG_var_pl(:,:,:,nq),       &
                          ADM_kall, ADM_kmin, ADM_kmax )
       val_min = GTL_min( DIAG_var   (:,:,:,nq),       &
                          DIAG_var_pl(:,:,:,nq),       &
                          ADM_kall, ADM_kmin, ADM_kmax )

       write(ADM_LOG_FID,'(1x,A,A16,2(A,1PE24.17))') '--- ', DIAG_name(nq), ': max=', val_max, ', min=', val_min
    enddo

    do nq = 1, TRC_vmax
       val_max = GTL_max( DIAG_var   (:,:,:,DIAG_vmax0+nq), &
                          DIAG_var_pl(:,:,:,DIAG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax      )

       if ( val_max <= 0.D0 ) then
          nonzero = .false.
       else
          nonzero = .true.
       endif

       val_min = GTL_min( DIAG_var   (:,:,:,DIAG_vmax0+nq), &
                          DIAG_var_pl(:,:,:,DIAG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax, nonzero)

       write(ADM_LOG_FID,'(1x,A,A16,2(A,1PE24.17))') '--- ', TRC_name(nq),  ': max=', val_max, ', min=', val_min
    enddo

    ! -> [add] H.Yashiro 20110819
    if ( output_io_mode == 'ADVANCED' ) then

       do nq = 1, DIAG_vmax0
          call FIO_output( DIAG_var(:,:,:,nq), basename, desc, '', DIAG_name(nq), DLABEL(nq), '', DUNIT(nq), &
                           FIO_REAL8, layername, 1, ADM_kall, 1, TIME_CTIME, TIME_CTIME                      )
       enddo

       do nq = 1, TRC_vmax
          call FIO_output( DIAG_var(:,:,:,DIAG_vmax0+nq), basename, desc, '', TRC_name(nq), WLABEL(nq), '', WUNIT, &
                           FIO_REAL8, layername, 1, ADM_kall, 1, TIME_CTIME, TIME_CTIME                            )
       enddo

    elseif( output_io_mode == 'LEGACY' ) then

       if ( output_direct_access ) then !--- direct access ( defalut )

          do l = 1, ADM_lall
             rgnid = ADM_prc_tab(l,ADM_prc_me)
             call MISC_make_idstr(fname,trim(basename),'rgn',rgnid)
             fid = MISC_get_available_fid()
             open( unit   = fid,                 &
                   file   = trim(fname),         &
                   form   = 'unformatted',       &
                   access = 'direct',            &
                   recl   = ADM_gall*ADM_kall*8, &
                   status = 'unknown'            )

             do nq = 1, DIAG_vmax
                write(fid,rec=nq) DIAG_var(:,:,l,nq)
             enddo

             close(fid)
          enddo

       else !--- sequential access

          do l = 1, ADM_lall
             rgnid = ADM_prc_tab(l,ADM_prc_me)
             call MISC_make_idstr(fname,trim(basename),'rgn',rgnid)
             fid = MISC_get_available_fid()
             open( unit   = fid,                 &
                   file   = trim(fname),         &
                   form   = 'unformatted',       &
                   access = 'sequential',        &
                   status = 'unknown'            )

             do nq = 1, DIAG_vmax
                write(fid) DIAG_var(:,:,l,nq)
             enddo

             close(fid)
          enddo

          if ( ADM_prc_me == ADM_prc_pl ) then
             fname = trim(basename)//'.pl'
             fid = misc_get_available_fid()
             open( unit   = fid,                 &
                   file   = trim(fname),         &
                   form   = 'unformatted',       &
                   access = 'sequential',        &
                   status = 'unknown'            )

             do nq = 1, DIAG_vmax
                write(fid) DIAG_var_pl(:,:,:,nq)
             enddo

             close(fid)
          endif

       endif !--- direct/sequencial

    endif  !--- io_mode
    ! <- [add] H.Yashiro 20110819

    return
  end subroutine restart_output

  !-----------------------------------------------------------------------------
  subroutine cnvvar_diag2prg( &
       prg,  prg_pl, & !--- [OUT]
       diag, diag_pl ) !--- [IN]
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_grd, only:  &
       GRD_afac, &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_GSGAM2,    &
       VMTR_GSGAM2_pl, &
       VMTR_GSGAM2H,   &
       VMTR_GSGAM2H_pl
    use mod_runconf, only: &
       TRC_vmax, &
       NQW_STR,  &
       NQW_END
    use mod_thrmdyn, only: &
       THRMDYN_qd_ijkl,  &
       THRMDYN_rho_ijkl, &
       THRMDYN_ein_ijkl
    implicit none

    real(8), intent(out) :: prg    (ADM_gall,   ADM_kall,ADM_lall,   PRG_vmax )
    real(8), intent(out) :: prg_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax )
    real(8), intent(in)  :: diag   (ADM_gall,   ADM_kall,ADM_lall,   DIAG_vmax)
    real(8), intent(in)  :: diag_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,DIAG_vmax)

    real(8) :: qd       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: qd_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rho      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ein      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: ein_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: rhog_h   (ADM_gall,   ADM_kall)
    real(8) :: rhog_h_pl(ADM_gall_pl,ADM_kall)

    integer :: ij, k, l, iv
    !---------------------------------------------------------------------------

    !--- calculation of dry mass concentration
    call THRMDYN_qd_ijkl ( ADM_gall, ADM_kall, ADM_lall, & !--- [IN]
                           TRC_vmax, NQW_STR, NQW_END,   & !--- [IN]
                           qd  (:,:,:),                  & !--- [OUT]
                           diag(:,:,:,I_qstr:I_qend)     ) !--- [IN]

    !--- calculation  of density 
    call THRMDYN_rho_ijkl( ADM_gall, ADM_kall, ADM_lall, & !--- [IN]
                           rho(:,:,:),                   & !--- [OUT]
                           diag(:,:,:,I_pre),            & !--- [IN]
                           diag(:,:,:,I_tem),            & !--- [IN]
                           qd  (:,:,:),                  & !--- [IN]
                           diag(:,:,:,I_qstr)            ) !--- [IN]

    !--- calculation of internal energy
    call THRMDYN_ein_ijkl( ADM_gall, ADM_kall, ADM_lall, & !--- [IN]
                           TRC_vmax, NQW_STR, NQW_END,   & !--- [IN]
                           ein(:,:,:),                   & !--- [OUT]
                           diag(:,:,:,I_tem),            & !--- [IN]
                           qd  (:,:,:),                  & !--- [IN]
                           diag(:,:,:,I_qstr:I_qend)     ) !--- [IN]

    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do ij = 1, ADM_gall
       prg(ij,k,l,I_RHOG  ) = rho(ij,k,l) * VMTR_GSGAM2(ij,k,l)
       prg(ij,k,l,I_RHOGVX) = prg(ij,k,l,I_RHOG) * diag(ij,k,l,I_vx)
       prg(ij,k,l,I_RHOGVY) = prg(ij,k,l,I_RHOG) * diag(ij,k,l,I_vy)
       prg(ij,k,l,I_RHOGVZ) = prg(ij,k,l,I_RHOG) * diag(ij,k,l,I_vz)
       prg(ij,k,l,I_RHOGE ) = prg(ij,k,l,I_RHOG) * ein (ij,k,l)
    enddo
    enddo
    enddo

    do iv = 1, TRC_vmax
    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do ij = 1, ADM_gall
       prg(ij,k,l,PRG_vmax0+iv) = prg(ij,k,l,I_RHOG) * diag(ij,k,l,DIAG_vmax0+iv)
    enddo
    enddo
    enddo
    enddo

    do l  = 1, ADM_lall
       !------ interpolation of rhog_h
       do k  = 2, ADM_kall
       do ij = 1, ADM_gall
          rhog_h(ij,k) = VMTR_GSGAM2H(ij,k,l) * 0.5D0 * ( GRD_afac(k) * rho(ij,k  ,l) &
                                                        + GRD_bfac(k) * rho(ij,k-1,l) )
       enddo
       enddo
       do ij = 1, ADM_gall
          rhog_h(ij,1) = rhog_h(ij,2)
       enddo

       do k  = 1, ADM_kall
       do ij = 1, ADM_gall
          prg(ij,k,l,I_RHOGW) = rhog_h(ij,k) * diag(ij,k,l,I_w)
       enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then

       !--- calculation of dry mass concentration
       call THRMDYN_qd_ijkl ( ADM_gall_pl, ADM_kall, ADM_lall_pl, & !--- [IN]
                              TRC_vmax, NQW_STR, NQW_END,         & !--- [IN]
                              qd_pl  (:,:,:),                     & !--- [OUT]
                              diag_pl(:,:,:,I_qstr:I_qend)        ) !--- [IN]

       !--- calculation  of density 
       call THRMDYN_rho_ijkl( ADM_gall_pl, ADM_kall, ADM_lall_pl, & !--- [IN]
                              rho_pl(:,:,:),                      & !--- [OUT]
                              diag_pl(:,:,:,I_pre),               & !--- [IN]
                              diag_pl(:,:,:,I_tem),               & !--- [IN]
                              qd_pl  (:,:,:),                     & !--- [IN]
                              diag_pl(:,:,:,I_qstr)               ) !--- [IN]

       !--- calculation of internal energy
       call THRMDYN_ein_ijkl( ADM_gall_pl, ADM_kall, ADM_lall_pl, & !--- [IN]
                              TRC_vmax, NQW_STR, NQW_END,         & !--- [IN]
                              ein_pl(:,:,:),                      & !--- [OUT]
                              diag_pl(:,:,:,I_tem),               & !--- [IN]
                              qd_pl  (:,:,:),                     & !--- [IN]
                              diag_pl(:,:,:,I_qstr:I_qend)        ) !--- [IN]

       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do ij = 1, ADM_gall_pl
          prg_pl(ij,k,l,I_RHOG  ) = rho_pl(ij,k,l) * VMTR_GSGAM2_pl(ij,k,l)
          prg_pl(ij,k,l,I_RHOGVX) = prg_pl(ij,k,l,I_RHOG) * diag_pl(ij,k,l,I_vx)
          prg_pl(ij,k,l,I_RHOGVY) = prg_pl(ij,k,l,I_RHOG) * diag_pl(ij,k,l,I_vy)
          prg_pl(ij,k,l,I_RHOGVZ) = prg_pl(ij,k,l,I_RHOG) * diag_pl(ij,k,l,I_vz)
          prg_pl(ij,k,l,I_RHOGE ) = prg_pl(ij,k,l,I_RHOG) * ein_pl (ij,k,l)
       enddo
       enddo
       enddo

       do iv = 1, TRC_vmax
       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do ij = 1, ADM_gall_pl
          prg_pl(ij,k,l,PRG_vmax0+iv) = prg_pl(ij,k,l,I_RHOG) * diag_pl(ij,k,l,DIAG_vmax0+iv)
       enddo
       enddo
       enddo
       enddo

       do l  = 1, ADM_lall_pl
          !------ interpolation of rhog_h
          do k  = 2, ADM_kall
          do ij = 1, ADM_gall_pl
             rhog_h_pl(ij,k) = VMTR_GSGAM2H_pl(ij,k,l) * 0.5D0 * ( GRD_afac(k) * rho_pl(ij,k  ,l) &
                                                                 + GRD_bfac(k) * rho_pl(ij,k-1,l) )
          enddo
          enddo
          do ij = 1, ADM_gall_pl
             rhog_h_pl(ij,1) = rhog_h_pl(ij,2)
          enddo

          do k  = 1, ADM_kall
          do ij = 1, ADM_gall_pl
             prg_pl(ij,k,l,I_RHOGW) = rhog_h_pl(ij,k) * diag_pl(ij,k,l,I_w)
          enddo
          enddo
       enddo

    endif

    return
  end subroutine cnvvar_diag2prg

  !-----------------------------------------------------------------------------
  subroutine cnvvar_prg2diag(&
       prg,  prg_pl, & !--- [OUT]
       diag, diag_pl ) !--- [IN]
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_grd, only:  &
       GRD_afac, &
       GRD_bfac
    use mod_vmtr, only : &
       VMTR_GSGAM2,    &
       VMTR_GSGAM2_pl, &
       VMTR_GSGAM2H,   &
       VMTR_GSGAM2H_pl
    use mod_runconf, only: &
       TRC_vmax, &
       NQW_STR,  &
       NQW_END
    use mod_thrmdyn, only: &
       THRMDYN_qd_ijkl,    &
       THRMDYN_tempre_ijkl
    implicit none

    real(8), intent(in)  :: prg    (ADM_gall,   ADM_kall,ADM_lall,   PRG_vmax )
    real(8), intent(in)  :: prg_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax )
    real(8), intent(out) :: diag   (ADM_gall,   ADM_kall,ADM_lall,   DIAG_vmax)
    real(8), intent(out) :: diag_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,DIAG_vmax)

    real(8) :: qd       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: qd_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rho      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ein      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: ein_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: rhog_h   (ADM_gall,   ADM_kall)
    real(8) :: rhog_h_pl(ADM_gall_pl,ADM_kall)

    integer :: ij, k, l, iv
    !---------------------------------------------------------------------------

    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do ij = 1, ADM_gall
       rho (ij,k,l)      = prg(ij,k,l,I_RHOG  ) / VMTR_GSGAM2(ij,k,l)
       diag(ij,k,l,I_vx) = prg(ij,k,l,I_RHOGVX) / prg(ij,k,l,I_RHOG)
       diag(ij,k,l,I_vy) = prg(ij,k,l,I_RHOGVY) / prg(ij,k,l,I_RHOG)
       diag(ij,k,l,I_vz) = prg(ij,k,l,I_RHOGVZ) / prg(ij,k,l,I_RHOG)
       ein (ij,k,l)      = prg(ij,k,l,I_RHOGE ) / prg(ij,k,l,I_RHOG)
    enddo
    enddo
    enddo

    do iv = 1, TRC_vmax
    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do ij = 1, ADM_gall
       diag(ij,k,l,DIAG_vmax0+iv) = prg(ij,k,l,PRG_vmax0+iv) / prg(ij,k,l,I_RHOG)
    enddo
    enddo
    enddo
    enddo

    !--- calculation of dry mass concentration
    call THRMDYN_qd_ijkl ( ADM_gall, ADM_kall, ADM_lall, & !--- [IN]
                           TRC_vmax, NQW_STR, NQW_END,   & !--- [IN]
                           qd  (:,:,:),                  & !--- [OUT]
                           diag(:,:,:,I_qstr:I_qend)     ) !--- [IN]

    !--- calculation of tem, pre
    call THRMDYN_tempre_ijkl( ADM_gall, ADM_kall, ADM_lall, & !--- [IN]
                              TRC_vmax, NQW_STR, NQW_END,   & !--- [IN]
                              diag(:,:,:,I_tem),            & !--- [OUT]
                              diag(:,:,:,I_pre),            & !--- [OUT]
                              ein(:,:,:),                   & !--- [IN]
                              rho(:,:,:),                   & !--- [IN]
                              qd  (:,:,:),                  & !--- [IN]
                              diag(:,:,:,I_qstr:I_qend)     ) !--- [IN]

    do l  = 1, ADM_lall
       !------ interpolation of rhog_h
       do k  = 2, ADM_kall
       do ij = 1, ADM_gall
          rhog_h(ij,k) = VMTR_GSGAM2H(ij,k,l) * 0.5D0 * ( GRD_afac(k) * rho(ij,k  ,l) &
                                                        + GRD_bfac(k) * rho(ij,k-1,l) )
       enddo
       enddo
       do ij = 1, ADM_gall
          rhog_h(ij,1) = rhog_h(ij,2)
       enddo

       do k  = 1, ADM_kall
       do ij = 1, ADM_gall
          diag(ij,k,l,I_w) = prg(ij,k,l,I_RHOGW) / rhog_h(ij,k)
       enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then

       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do ij = 1, ADM_gall_pl
          rho_pl (ij,k,l)      = prg_pl(ij,k,l,I_RHOG  ) / VMTR_GSGAM2_pl(ij,k,l)
          diag_pl(ij,k,l,I_vx) = prg_pl(ij,k,l,I_RHOGVX) / prg_pl(ij,k,l,I_RHOG)
          diag_pl(ij,k,l,I_vy) = prg_pl(ij,k,l,I_RHOGVY) / prg_pl(ij,k,l,I_RHOG)
          diag_pl(ij,k,l,I_vz) = prg_pl(ij,k,l,I_RHOGVZ) / prg_pl(ij,k,l,I_RHOG)
          ein_pl (ij,k,l)      = prg_pl(ij,k,l,I_RHOGE ) / prg_pl(ij,k,l,I_RHOG)
       enddo
       enddo
       enddo

       do iv = 1, TRC_vmax
       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do ij = 1, ADM_gall_pl
          diag_pl(ij,k,l,DIAG_vmax0+iv) = prg_pl(ij,k,l,PRG_vmax0+iv) / prg_pl(ij,k,l,I_RHOG)
       enddo
       enddo
       enddo
       enddo

       !--- calculation of dry mass concentration
       call THRMDYN_qd_ijkl ( ADM_gall_pl, ADM_kall, ADM_lall_pl, & !--- [IN]
                              TRC_vmax, NQW_STR, NQW_END,         & !--- [IN]
                              qd_pl  (:,:,:),                     & !--- [OUT]
                              diag_pl(:,:,:,I_qstr:I_qend)        ) !--- [IN]

       !--- calculation of tem, pre
       call THRMDYN_tempre_ijkl( ADM_gall_pl, ADM_kall, ADM_lall_pl, & !--- [IN]
                                 TRC_vmax, NQW_STR, NQW_END,         & !--- [IN]
                                 diag_pl(:,:,:,I_tem),               & !--- [OUT]
                                 diag_pl(:,:,:,I_pre),               & !--- [OUT]
                                 ein_pl(:,:,:),                      & !--- [IN]
                                 rho_pl(:,:,:),                      & !--- [IN]
                                 qd_pl  (:,:,:),                     & !--- [IN]
                                 diag_pl(:,:,:,I_qstr:I_qend)        ) !--- [IN]

       do l  = 1, ADM_lall_pl
          !------ interpolation of rhog_h
          do k  = 2, ADM_kall
          do ij = 1, ADM_gall_pl
             rhog_h_pl(ij,k) = VMTR_GSGAM2H_pl(ij,k,l) * 0.5D0 * ( GRD_afac(k) * rho_pl(ij,k  ,l) &
                                                                 + GRD_bfac(k) * rho_pl(ij,k-1,l) )
          enddo
          enddo
          do ij = 1, ADM_gall_pl
             rhog_h_pl(ij,1) = rhog_h_pl(ij,2)
          enddo

          do k  = 1, ADM_kall
          do ij = 1, ADM_gall_pl
             diag_pl(ij,k,l,I_w) = prg_pl(ij,k,l,I_RHOGW) / rhog_h_pl(ij,k)
          enddo
          enddo
       enddo

    endif

    return
  end subroutine cnvvar_prg2diag

  !-----------------------------------------------------------------------------
  ! 09/04/14 [Add] T.Mitsui
  ! GCMs, including Reanalysis data, diagnose qi with qc and temperature implicitly
  ! and qi is never provided by them.
  ! Here we diagnose qi without release of latent heat.
  ! This assumption is equivalent with treatment of GCMs who never predict qi.
  !
  ! Therefore we should diagnose in this module before calculating rhoge
  ! because we should keep coservation laws for mass and energy.
  subroutine diag_qi( &
       gall, kall, lall, & ! in
       tem,                   & ! in
       qc, qi, rhogqc, rhogqi ) ! inout
    use mod_adm, only: &
       ADM_CTL_FID,  &
       ADM_LOG_FID
    implicit none

    integer, intent(in)    :: gall
    integer, intent(in)    :: kall
    integer, intent(in)    :: lall
    real(8), intent(in)    :: tem(gall,kall,lall)
    real(8), intent(inout) :: qc(gall,kall,lall)
    real(8), intent(inout) :: qi(gall,kall,lall)
    real(8), intent(inout) :: rhogqc(gall,kall,lall)
    real(8), intent(inout) :: rhogqi(gall,kall,lall)
    ! Reference: MIROC4.1, G98, Rogers and Yau(1989)(book) and so on.
    ! -15deg. is the most effective temperature of Bergeron process
    real(8), save :: tem_low = 258.15d0 
    !   0deg. is begining of Bergeron process
    real(8), save :: tem_up  = 273.15d0 
    !
    real(8), parameter :: dtem_min=0.1d0
    real(8), parameter :: tem_low_min = 173.15d0
    !
    namelist /nm_restart_diag_qi/ &
         tem_low, tem_up
    !
    real(8) :: r_dtem
    real(8) :: liquid_ratio
    !
    real(8) :: qc_pre(gall,kall,lall)
    real(8) :: rhogqc_pre(gall,kall,lall)
    real(8) :: dqc(gall,kall,lall)
    real(8) :: drhogqc(gall,kall,lall)
    !
    integer :: ij,k,l
    !---------------------------------------------------------------------------

    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=nm_restart_diag_qi,end=100)
100 continue
    ! 09/04/14 [Add] T.Mitsui, filter
    tem_low =  max(tem_low_min, min(tem_low, tem_up-dtem_min))
    tem_up  =  max(tem_up, tem_low+dtem_min)
    write(ADM_LOG_FID,nml=nm_restart_diag_qi)
    !
    r_dtem = 1.d0/(tem_up - tem_low)
    qc_pre(:,:,:)     = qc(:,:,:)
    rhogqc_pre(:,:,:) = rhogqc(:,:,:)
    !
    do l=1, lall
       do k=1, kall
          do ij=1, gall
             ! ratio=1 => all water,  ratio=0 => all ice.
             liquid_ratio  = min(1.d0, max(0.d0, (tem(ij,k,l)-tem_low)*r_dtem))
             dqc(ij,k,l)     = (liquid_ratio-1.d0)*qc(ij,k,l)
             drhogqc(ij,k,l) = (liquid_ratio-1.d0)*rhogqc(ij,k,l)
             !
             qc(ij,k,l)      = qc(ij,k,l)     + dqc(ij,k,l)
             qi(ij,k,l)      = qi(ij,k,l)     - dqc(ij,k,l)
             rhogqc(ij,k,l)  = rhogqc(ij,k,l) + drhogqc(ij,k,l)
             rhogqi(ij,k,l)  = rhogqi(ij,k,l) - drhogqc(ij,k,l)
          end do
       end do
    end do
    !
    write(ADM_LOG_FID,'(a)') "*** Diagnosis of QI, check max and min value "
    do k=1, kall
       write(ADM_LOG_FID,'(a,i5)') "Layer number  = ", k
       write(ADM_LOG_FID,'(a,2f16.6)') "tem(max,min)  = ",maxval(tem(:,k,:)), minval(tem(:,k,:))
       write(ADM_LOG_FID,'(a,e16.6,a,e16.6)') "qc(pre) max   =",maxval(qc_pre(:,k,:)),&
                                              " => qc(post) max   =" , maxval(qc(:,k,:))
       write(ADM_LOG_FID,'(a,e16.6,a,e16.6)') "qc(pre) min   =",minval(qc_pre(:,k,:)),&
                                              " => qc(post) min   =" , minval(qc(:,k,:))
       write(ADM_LOG_FID,'(a,e16.6,a,e16.6)') "dqc     min   =",minval(dqc(:,k,:)),&
                                              ",   qi(post) max   =" , maxval(qi(:,k,:))
    enddo

    return
  end subroutine diag_qi

  !-----------------------------------------------------------------------------
  ! [Add] 09/08/18 T.MItsui. convert qc and qi into qv(immediate evaporation).
  subroutine convert_qcqi_to_qv( &
       ijdim, kdim, lall, nqmax, &
       tem, q, rhogq )
    use mod_runconf, only: &
         I_QV, I_QC, I_QI, &
         I_NC, I_NI,       &
         NQW_STR, NQW_END, &
         opt_2moment_water
    use mod_cnst, only: &
         CNST_LH00,     &
         CNST_LHF00,    &
         CNST_CVV,      &
         CNST_CL,       &
         CNST_CI
    use mod_thrmdyn, only: &
         thrmdyn_cv
    implicit none
    !
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: lall
    integer, intent(in) :: nqmax
    real(8), intent(inout) :: tem(ijdim,kdim,lall)
    real(8), intent(inout) :: q(ijdim,kdim,lall,nqmax)
    real(8), intent(inout) :: rhogq(ijdim,kdim,lall,nqmax)
    !
    real(8) :: qd(ijdim,kdim)
    real(8) :: cva(ijdim,kdim)
    real(8) :: dqv(ijdim,kdim), dqi(ijdim,kdim)
    real(8) :: drhogqv(ijdim,kdim)
    real(8) :: drhoe(ijdim,kdim)
    !
    integer :: l_region,nq
    !---------------------------------------------------------------------------

    do l_region=1,lall
       qd(:,:) = 1.d0
       do nq=NQW_STR, NQW_END
          qd(:,:) = qd(:,:) - q(:,:,l_region,nq)
       end do
       call thrmdyn_cv( &
            ijdim,      & !--- in
            cva,        & !--- out
            q,qd        ) !--- in
       !
       if( NQW_END >= I_QI )then
          dqv(:,:)     =   q(:,:,l_region,I_QC) + q(:,:,l_region,I_QI)
          dqi(:,:)     = - q(:,:,l_region,I_QI)
          drhogqv(:,:) = rhogq(:,:,l_region,I_QC) + rhogq(:,:,l_region,I_QI)
          drhoe(:,:)   = -CNST_LH00*dqv(:,:) + CNST_LHF00*dqi(:,:) &
               - (CNST_CVV*tem(:,:,l_region)-CNST_CL*tem(:,:,l_region))*dqv(:,:) &
               - (CNST_CI *tem(:,:,l_region)-CNST_CL*tem(:,:,l_region))*dqi(:,:)
          !
          tem(:,:,l_region)    = tem(:,:,l_region) + drhoe(:,:)/cva(:,:)
          q(:,:,l_region,I_QV) = q(:,:,l_region,I_QV) + dqv(:,:)
          q(:,:,l_region,I_QC) = 0.d0
          q(:,:,l_region,I_QI) = 0.d0
          rhogq(:,:,l_region,I_QV) = rhogq(:,:,l_region,I_QV) + drhogqv(:,:)
          rhogq(:,:,l_region,I_QC) = 0.d0
          rhogq(:,:,l_region,I_QI) = 0.d0          
          if( opt_2moment_water )then
             q(:,:,l_region,I_NC)     = 0.d0
             q(:,:,l_region,I_NI)     = 0.d0
             rhogq(:,:,l_region,I_NC) = 0.d0
             rhogq(:,:,l_region,I_NI) = 0.d0          
          end if
       else
          dqv(:,:)     = q(:,:,l_region,I_QC) 
          drhogqv(:,:) = rhogq(:,:,l_region,I_QC) 
          !
          drhoe(:,:)   = -CNST_LH00*dqv(:,:) &
               - (CNST_CVV*tem(:,:,l_region)-CNST_CL*tem(:,:,l_region))*dqv(:,:) 
          !
          tem(:,:,l_region)    = tem(:,:,l_region) + drhoe(:,:)/cva(:,:)
          q(:,:,l_region,I_QV) = q(:,:,l_region,I_QV) + dqv(:,:)
          q(:,:,l_region,I_QC) = 0.d0
          rhogq(:,:,l_region,I_QV) = rhogq(:,:,l_region,I_QV) + drhogqv(:,:)
          rhogq(:,:,l_region,I_QC) = 0.d0
          if( opt_2moment_water )then
             q(:,:,l_region,I_NC) = 0.d0
             rhogq(:,:,l_region,I_NC) = 0.d0
          end if
       end if
    end do

    return
  end subroutine convert_qcqi_to_qv    

end module mod_prgvar
!-------------------------------------------------------------------------------
