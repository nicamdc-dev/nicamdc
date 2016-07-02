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
  use mod_precision
  use mod_stdio
  use mod_debug
  use mod_runconf, only: &
     PRG_vmax0,  &
     I_RHOG,     &
     I_RHOGVX,   &
     I_RHOGVY,   &
     I_RHOGVZ,   &
     I_RHOGW,    &
     I_RHOGE,    &
     I_RHOGQstr, &
     I_RHOGQend, &
     DIAG_vmax0, &
     I_pre,      &
     I_tem,      &
     I_vx,       &
     I_vy,       &
     I_vz,       &
     I_w,        &
     I_qstr,     &
     I_qend
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

  public :: restart_input
  public :: restart_output

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: PRG_var (:,:,:,:)
  real(RP), public, allocatable :: PRG_var1(:,:,:,:)
  real(RP), public, allocatable :: DIAG_var(:,:,:,:)

  character(len=H_LONG), public :: restart_input_basename  = ''
  character(len=H_LONG), public :: restart_output_basename = ''

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: PRG_var_pl (:,:,:,:)
  real(RP), private, allocatable :: PRG_var1_pl(:,:,:,:)
  real(RP), private, allocatable :: DIAG_var_pl(:,:,:,:)

  integer, private :: TRC_vmax_input ! number of input tracer variables

  character(len=H_LONG), private :: layername      = ''
  character(len=H_LONG), private :: input_io_mode  = 'ADVANCED'
  character(len=H_LONG), private :: output_io_mode = 'ADVANCED'
  logical,               private :: allow_missingq = .false.

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine prgvar_setup
  !>
  subroutine prgvar_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_kall,      &
       ADM_lall,      &
       ADM_lall_pl
    use mod_runconf, only: &
       PRG_vmax,  &
       DIAG_vmax, &
       TRC_vmax
    implicit none

    character(len=H_LONG) :: input_basename    = ''
    character(len=H_LONG) :: output_basename   = 'restart'
    character(len=H_LONG) :: restart_layername = ''

    namelist / RESTARTPARAM / &
       TRC_vmax_input,    &
       input_basename,    &
       output_basename,   &
       restart_layername, &
       input_io_mode,     &
       output_io_mode,    &
       allow_missingq

    integer :: ierr
    !---------------------------------------------------------------------------

    TRC_vmax_input = TRC_vmax

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[prgvar]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=RESTARTPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** RESTARTPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist RESTARTPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist RESTARTPARAM. STOP.'
       call PRC_MPIstop
    endif
    write(IO_FID_LOG,nml=RESTARTPARAM)

    restart_input_basename  = input_basename
    restart_output_basename = output_basename
    layername               = restart_layername

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** io_mode for restart, input : ', trim(input_io_mode)
    if    ( input_io_mode == 'ADVANCED'     ) then
    elseif( input_io_mode == 'LEGACY'       ) then
    elseif( input_io_mode == 'IDEAL'        ) then
    elseif( input_io_mode == 'IDEAL_TRACER' ) then
    else
       write(IO_FID_LOG,*) 'xxx Invalid input_io_mode. STOP.'
       call PRC_MPIstop
    endif

    write(IO_FID_LOG,*) '*** io_mode for restart, output: ', trim(output_io_mode)
    if    ( output_io_mode == 'ADVANCED' ) then
    elseif( output_io_mode == 'LEGACY'   ) then
    else
       write(IO_FID_LOG,*) 'xxx Invalid output_io_mode. STOP'
       call PRC_MPIstop
    endif

    if ( allow_missingq ) then
       write(IO_FID_LOG,*) '*** Allow missing tracer in restart file.'
       write(IO_FID_LOG,*) '*** Value will set to zero for missing tracer.'
    endif

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
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(RP), intent(out) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogq    (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP), intent(out) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
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

       if ( ADM_have_pl ) then

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

       if ( ADM_have_pl ) then

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
  subroutine prgvar_get_noq( &
       rhog,   rhog_pl,   &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       rhoge,  rhoge_pl   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    implicit none

    real(RP), intent(out) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

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

    if ( ADM_have_pl ) then

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
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_vmtr, only: &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl
    use mod_runconf, only: &
       TRC_vmax
    use mod_cnvvar, only: &
       cnvvar_prg2diag
    implicit none

    real(RP), intent(out) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogq    (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP), intent(out) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(RP), intent(out) :: rho      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: pre      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: tem      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vx       (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vy       (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vz       (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: w        (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: q        (ADM_gall,   ADM_kall,ADM_lall,   TRC_vmax)
    real(RP), intent(out) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

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

       rho   (n,k,l) = PRG_var(n,k,l,I_RHOG) * VMTR_RGSGAM2(n,k,l)

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

    if ( ADM_have_pl ) then

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do n = 1, ADM_gall_pl
          rhog_pl  (n,k,l) = PRG_var_pl(n,k,l,I_RHOG  )
          rhogvx_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVX)
          rhogvy_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVY)
          rhogvz_pl(n,k,l) = PRG_var_pl(n,k,l,I_RHOGVZ)
          rhogw_pl (n,k,l) = PRG_var_pl(n,k,l,I_RHOGW )
          rhoge_pl (n,k,l) = PRG_var_pl(n,k,l,I_RHOGE )

          rho_pl   (n,k,l) = PRG_var_pl(n,k,l,I_RHOG) * VMTR_RGSGAM2_pl(n,k,l)

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
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(RP), intent(in)  :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogq    (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP), intent(in)  :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    integer, intent(in)  :: num

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

       if ( ADM_have_pl ) then

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

       if ( ADM_have_pl ) then

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
       ADM_kall,        &
       ADM_lall,        &
       ADM_IopJop_nmax, &
       ADM_IopJop,      &
       ADM_GIoJo
    use mod_runconf, only: &
       TRC_vmax
    implicit none

    real(RP), intent(out) :: rhog  (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvx(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvy(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvz(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogw (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhoge (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogq (ADM_IopJop_nmax,ADM_kall,ADM_lall,TRC_vmax)

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
       ADM_kall,        &
       ADM_lall,        &
       ADM_IopJop_nmax, &
       ADM_IopJop,      &
       ADM_GIoJo
    use mod_vmtr, only: &
       VMTR_RGSGAM2
    use mod_runconf, only: &
       TRC_vmax
    use mod_cnvvar, only: &
       cnvvar_prg2diag
    implicit none

    real(RP), intent(out) :: rhog  (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvx(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvy(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvz(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogw (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhoge (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogq (ADM_IopJop_nmax,ADM_kall,ADM_lall,TRC_vmax)
    real(RP), intent(out) :: rho   (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: pre   (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: tem   (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: vx    (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: vy    (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: vz    (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: w     (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(out) :: q     (ADM_IopJop_nmax,ADM_kall,ADM_lall,TRC_vmax)

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

       rho   (n,k,l) = PRG_var(ij,k,l,I_RHOG) * VMTR_RGSGAM2(ij,k,l)

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
       ADM_kall,        &
       ADM_lall,        &
       ADM_IopJop_nmax, &
       ADM_IopJop,      &
       ADM_GIoJo
    use mod_comm, only: &
       COMM_var
    use mod_runconf, only: &
       PRG_vmax, &
       TRC_vmax
    implicit none

    real(RP), intent(in) :: rhog  (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogvx(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogvy(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogvz(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogw (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhoge (ADM_IopJop_nmax,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogq (ADM_IopJop_nmax,ADM_kall,ADM_lall,TRC_vmax)

    integer :: n, k, l, nq, ij
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
    call COMM_var( PRG_var, PRG_var_pl, ADM_kall, PRG_vmax )

    return
  end subroutine prgvar_set_in

  !-----------------------------------------------------------------------------
  subroutine restart_input( basename )
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
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
       PRG_name,  &
       DIAG_vmax, &
       DIAG_name, &
       TRC_vmax,  &
       TRC_name
    use mod_cnvvar, only: &
       cnvvar_diag2prg
    use mod_ideal_init, only :  & ! [add] R.Yoshida 20121019
       dycore_input, &
       tracer_input
    implicit none

    character(len=H_LONG), intent(in) :: basename

    character(len=H_LONG) :: fname

    real(RP) :: val_max, val_min
    logical :: nonzero

    integer :: fid
    integer :: l, rgnid, nq
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** read restart/initial data'

    if ( input_io_mode == 'ADVANCED' ) then

       do nq = 1, DIAG_vmax0
          call FIO_input( DIAG_var(:,:,:,nq),basename,DIAG_name(nq), &
                          layername,1,ADM_kall,1                     )
       enddo

       do nq = 1, TRC_vmax_input
          call FIO_input( DIAG_var(:,:,:,DIAG_vmax0+nq),basename,TRC_name(nq), &
                          layername,1,ADM_kall,1,                              &
                          allow_missingq=allow_missingq                        )
       enddo

    elseif( input_io_mode == 'LEGACY' ) then

       do l = 1, ADM_lall
          rgnid = ADM_prc_tab(l,ADM_prc_me)
          call IO_make_idstr(fname,trim(basename),'rgn',rgnid)
          fid = IO_get_available_fid()
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

    elseif( input_io_mode == 'IDEAL' ) then

       write(IO_FID_LOG,*) '*** make ideal initials'

       call dycore_input( DIAG_var(:,:,:,:) ) !--- [OUT]

    elseif( input_io_mode == 'IDEAL_TRACER' ) then

       do nq = 1, DIAG_vmax0
          call FIO_input( DIAG_var(:,:,:,nq),basename,DIAG_name(nq), &
                          layername,1,ADM_kall,1                     )
       enddo

       write(IO_FID_LOG,*) '*** make ideal initials for tracer'

       call tracer_input( DIAG_var(:,:,:,DIAG_vmax0+1:DIAG_vmax0+TRC_vmax) ) !--- [OUT]

    endif !--- io_mode

    call COMM_var( DIAG_var, DIAG_var_pl, ADM_kall, DIAG_vmax )

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '====== data range check : diagnostic variables ======'
    do nq = 1, DIAG_vmax0
       val_max = GTL_max( DIAG_var(:,:,:,nq), DIAG_var_pl(:,:,:,nq), ADM_kall, ADM_kmin, ADM_kmax )
       val_min = GTL_min( DIAG_var(:,:,:,nq), DIAG_var_pl(:,:,:,nq), ADM_kall, ADM_kmin, ADM_kmax )
       write(IO_FID_LOG,'(1x,A,A16,2(A,1PE24.17))') '--- ', DIAG_name(nq), ': max=', val_max, ', min=', val_min
    enddo

    do nq = 1, TRC_vmax
       val_max = GTL_max( DIAG_var   (:,:,:,DIAG_vmax0+nq), &
                          DIAG_var_pl(:,:,:,DIAG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax      )

       if ( val_max <= 0.0_RP ) then
          nonzero = .false.
       else
          nonzero = .true.
       endif

       val_min = GTL_min( DIAG_var   (:,:,:,DIAG_vmax0+nq), &
                          DIAG_var_pl(:,:,:,DIAG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax, nonzero)

       write(IO_FID_LOG,'(1x,A,A16,2(A,1PE24.17))') '--- ', TRC_name(nq),  ': max=', val_max, ', min=', val_min
    enddo

    call cnvvar_diag2prg( PRG_var (:,:,:,:), PRG_var_pl (:,:,:,:), & !--- [OUT]
                          DIAG_var(:,:,:,:), DIAG_var_pl(:,:,:,:)  ) !--- [IN]

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '====== data range check : prognostic variables ======'
    do nq = 1, PRG_vmax0
       val_max = GTL_max( PRG_var(:,:,:,nq), PRG_var_pl(:,:,:,nq), ADM_kall, ADM_kmin, ADM_kmax )
       val_min = GTL_min( PRG_var(:,:,:,nq), PRG_var_pl(:,:,:,nq), ADM_kall, ADM_kmin, ADM_kmax )
       write(IO_FID_LOG,'(1x,A,A16,2(A,1PE24.17))') '--- ', PRG_name(nq), ': max=', val_max, ', min=', val_min
    enddo

    do nq = 1, TRC_vmax
       val_max = GTL_max( PRG_var   (:,:,:,PRG_vmax0+nq), &
                          PRG_var_pl(:,:,:,PRG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax    )

       if ( val_max <= 0.0_RP ) then
          nonzero = .false.
       else
          nonzero = .true.
       endif

       val_min = GTL_min( PRG_var   (:,:,:,PRG_vmax0+nq),      &
                          PRG_var_pl(:,:,:,PRG_vmax0+nq),      &
                          ADM_kall, ADM_kmin, ADM_kmax, nonzero)

       write(IO_FID_LOG,'(1x,A,A16,2(A,1PE24.17))') '--- rhog * ', TRC_name(nq),  ': max=', val_max, ', min=', val_min
    enddo

    return
  end subroutine restart_input

  !-----------------------------------------------------------------------------
  subroutine restart_output( basename )
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_gall,    &
       ADM_lall,    &
       ADM_kall,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_io_param, only: &
       IO_REAL8
    use mod_fio, only: &
       FIO_output
    use mod_time, only : &
       TIME_CTIME
    use mod_gtl, only: &
       GTL_max, &
       GTL_min
    use mod_runconf, only: &
       DIAG_vmax, &
       DIAG_name, &
       TRC_vmax,  &
       TRC_name,  &
       WLABEL
    use mod_cnvvar, only: &
       cnvvar_prg2diag
    implicit none

    character(len=H_LONG), intent(in) :: basename

    character(len=H_MID)   :: desc = 'INITIAL/RESTART_data_of_prognostic_variables' ! [add] H.Yashiro 20110819

    character(len=H_SHORT) :: DLABEL(DIAG_vmax0)
    data DLABEL / 'Pressure ',        &
                  'Temperature ',     &
                  'H-Velocity(XDIR)', &
                  'H-Velocity(YDIR)', &
                  'H-Velocity(ZDIR)', &
                  'V-Velocity '       /

    character(len=H_SHORT) :: DUNIT(DIAG_vmax0)
    data DUNIT /  'Pa',  &
                  'K',   &
                  'm/s', &
                  'm/s', &
                  'm/s', &
                  'm/s'  /

    character(len=H_SHORT) :: WUNIT = 'kg/kg'

    character(len=H_LONG) :: fname

    real(RP) :: val_max, val_min
    logical  :: nonzero

    integer  :: fid
    integer  :: l, rgnid, nq
    !---------------------------------------------------------------------------

    call cnvvar_prg2diag( PRG_var (:,:,:,:), PRG_var_pl (:,:,:,:), & !--- [IN]
                          DIAG_var(:,:,:,:), DIAG_var_pl(:,:,:,:)  ) !--- [OUT]

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '====== data range check : prognostic variables ======'
    do nq = 1, DIAG_vmax0
       val_max = GTL_max( DIAG_var   (:,:,:,nq),       &
                          DIAG_var_pl(:,:,:,nq),       &
                          ADM_kall, ADM_kmin, ADM_kmax )
       val_min = GTL_min( DIAG_var   (:,:,:,nq),       &
                          DIAG_var_pl(:,:,:,nq),       &
                          ADM_kall, ADM_kmin, ADM_kmax )

       write(IO_FID_LOG,'(1x,A,A16,2(A,1PE24.17))') '--- ', DIAG_name(nq), ': max=', val_max, ', min=', val_min
    enddo

    do nq = 1, TRC_vmax
       val_max = GTL_max( DIAG_var   (:,:,:,DIAG_vmax0+nq), &
                          DIAG_var_pl(:,:,:,DIAG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax      )

       if ( val_max <= 0.0_RP ) then
          nonzero = .false.
       else
          nonzero = .true.
       endif

       val_min = GTL_min( DIAG_var   (:,:,:,DIAG_vmax0+nq), &
                          DIAG_var_pl(:,:,:,DIAG_vmax0+nq), &
                          ADM_kall, ADM_kmin, ADM_kmax, nonzero)

       write(IO_FID_LOG,'(1x,A,A16,2(A,1PE24.17))') '--- ', TRC_name(nq),  ': max=', val_max, ', min=', val_min
    enddo

    if ( output_io_mode == 'ADVANCED' ) then

       do nq = 1, DIAG_vmax0
          call FIO_output( DIAG_var(:,:,:,nq), basename, desc, '', DIAG_name(nq), DLABEL(nq), '', DUNIT(nq), &
                           IO_REAL8, layername, 1, ADM_kall, 1, TIME_CTIME, TIME_CTIME                       )
       enddo

       do nq = 1, TRC_vmax
          call FIO_output( DIAG_var(:,:,:,DIAG_vmax0+nq), basename, desc, '', TRC_name(nq), WLABEL(nq), '', WUNIT, &
                           IO_REAL8, layername, 1, ADM_kall, 1, TIME_CTIME, TIME_CTIME                             )
       enddo

    elseif( output_io_mode == 'LEGACY' ) then

       do l = 1, ADM_lall
          rgnid = ADM_prc_tab(l,ADM_prc_me)
          call IO_make_idstr(fname,trim(basename),'rgn',rgnid)
          fid = IO_get_available_fid()
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

    endif  !--- io_mode

    return
  end subroutine restart_output

end module mod_prgvar
!-------------------------------------------------------------------------------
