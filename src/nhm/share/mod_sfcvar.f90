!-------------------------------------------------------------------------------
!
!+  Diagnostic variable module
!
!-------------------------------------------------------------------------------
module mod_sfcvar
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module contains the surface diagnostic variables and
  !       control subroutines for the non-hydrostatic model.
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Add this module
  !      0.05      xx-xx-xx   Change module name.
  !                05-10-21   M.Satoh set default of I_CUMFRC = 0
  !                06-04-18   T.Mitsui add precip_energy_cp,mp and sfcvar_set2
  !                07-05-08   H.Tomita : Delete I_TKE_SFC
  !                07-06-27   T.Mitsui : add sfcvar_set1
  !                07-07-23   K.Suzuki : Add I_CEV for use in SPRINTARS
  !                09-07-10   H.Tomita: Add several variables for energy budget.
  !                10-06-19   A.T.Noda:
  !                                Allow to use a convection parameterization
  !                                with an advanced microphysics schemes,
  !                                such as G98, NSW?,
  !                11-11-28   Y.Yamada:
  !                                Merge Terai-san timer code
  !                                                     into the original code.
  !                12-03-28   T.Seiki: add friction velocity used in SPRINTARS
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: sfcvar_setup
  public :: sfcvar_get
  public :: sfcvar_get1
  public :: sfcvar_get2
  public :: sfcvar_set
  public :: sfcvar_set1     ! add T.Mitsui 07.06.27
  public :: sfcvar_get_in
  public :: sfcvar_set_in
  public :: sfcvar_set1_in
  public :: sfcvar_get1_in
  public :: sfcvar_set2_in
  public :: sfcvar_get2_in
  public :: sfcvar_set2     ! add T.Mitsui 06.04.18
  public :: sfcvar_comm
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  integer, public, parameter :: INDEX_LAND    =  2
  integer, public, parameter :: INDEX_LANDICE =  1
  integer, public, parameter :: INDEX_SEA     =  0
  integer, public, parameter :: INDEX_SEAICE  = -1

  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public & Private parameters
  !
  integer, private, parameter :: DIAG_MAX = 1

  integer, public :: I_PRECIP_TOT
  integer, public :: I_PRECIP_CP
  integer, public :: I_PRECIP_MP
  integer, public :: I_OLR
  integer, public :: I_EVAP_SFC
  integer, public :: I_PRE_SFC
  integer, public :: I_TEM_SFC
  integer, public :: I_TH_SFC
  integer, public :: I_RHO_SFC
  integer, public :: I_RCOSZ
  integer, public :: I_RFLUXS_SU
  integer, public :: I_RFLUXS_SD
  integer, public :: I_RFLUXS_LU
  integer, public :: I_RFLUXS_LD
  integer, public :: I_RFLUX_TOA_SU
  integer, public :: I_RFLUX_TOA_SD
  integer, public :: I_RFLUX_TOA_LU
  integer, public :: I_RFLUX_TOA_LD
  integer, public :: I_RFLUX_TOA_SU_C
  integer, public :: I_RFLUX_TOA_SD_C
  integer, public :: I_RFLUX_TOA_LU_C
  integer, public :: I_RFLUX_TOA_LD_C
  integer, public :: I_QV_SFC
  integer, public :: I_QC_SFC
  integer, public :: I_QR_SFC
  integer, public :: I_SH_FLUX_SFC
  integer, public :: I_LH_FLUX_SFC
  integer, public :: I_TAUX_SFC
  integer, public :: I_TAUY_SFC
  integer, public :: I_TAUZ_SFC
  integer, public :: I_VX10
  integer, public :: I_VY10
  integer, public :: I_VZ10
  integer, public :: I_T2
  integer, public :: I_Q2
  integer, public :: I_CEV
  integer, public :: I_CUMFRC
  integer, public :: I_ALBEDO_SFC
  integer, public :: I_PRECIP
  integer, public :: I_PRECIP_ENERGY
  integer, public :: I_PRECIP_ENERGY_CP
  integer, public :: I_PRECIP_ENERGY_MP
  integer, public :: I_PRECIP_EIN
  integer, public :: I_PRECIP_EIN_CP
  integer, public :: I_PRECIP_EIN_MP
  integer, public :: I_PRECIP_LH
  integer, public :: I_PRECIP_LH_CP
  integer, public :: I_PRECIP_LH_MP
  integer, public :: I_PRECIP_PHI
  integer, public :: I_PRECIP_PHI_CP
  integer, public :: I_PRECIP_PHI_MP
  integer, public :: I_PRECIP_KIN
  integer, public :: I_PRECIP_KIN_CP
  integer, public :: I_PRECIP_KIN_MP
  integer, public :: I_EVAP_ENERGY
  integer, public :: I_EVAP_EIN
  integer, public :: I_EVAP_LH
  integer, public :: I_EVAP_PHI
  integer, public :: I_EVAP_KIN
  integer, public :: I_SFCRAD_ENERGY
  integer, public :: I_TOARAD_ENERGY
  integer, public :: I_VFRICTION

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  integer, allocatable, private :: NVAR(:)

  integer, private :: I_ALL

  REAL(RP), public , allocatable :: sfcvar   (:,:,:,:)
  REAL(RP), private, allocatable :: sfcvar_pl(:,:,:,:)

  integer, private, allocatable :: KMAX(:)
  integer, public , allocatable :: KSTR(:)
  integer, private, allocatable :: KEND(:)
  integer, private              :: KSUM

  REAL(RP), private :: VMISS = -999.D30

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine sfcvar_setup
    !
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_GALL_PL,       &
         ADM_LALL_PL,       &
         ADM_KNONE,         &
         ADM_gall,          &
         ADM_gmin,          &
         ADM_gmax,          &
         ADM_lall,          &
         ADM_gall_1d
    use mod_runconf, only : &
         RAIN_TYPE,         &
         CP_TYPE,           & ! 10/06/15 A.T.Noda
         AE_TYPE,           & ! 12/03/28 T.Seiki
         NRBND,             &
         NRDIR
    !
    implicit none

    integer :: i
    !
    I_PRECIP_TOT         = 1   !--- total precip
    I_PRECIP_CP          = 2   !--- convective precip ( rain, snow )
    I_PRECIP_MP          = 3   !--- dynamical precip  ( rain, snow )
    I_PRECIP             = 4   !--- total precip      ( rain, snow )
    ! 06.04.18 add T.Mitsui
    I_PRECIP_ENERGY_CP   = 5
    I_PRECIP_ENERGY_MP   = 6
    I_PRECIP_ENERGY      = 7
    ! 09.07.03 add H.Tomita
    I_PRECIP_EIN         = 8
    I_PRECIP_EIN_CP      = 9
    I_PRECIP_EIN_MP      = 10
    I_PRECIP_LH          = 11
    I_PRECIP_LH_CP       = 12
    I_PRECIP_LH_MP       = 13
    I_PRECIP_PHI         = 14
    I_PRECIP_KIN         = 15
    I_PRECIP_PHI_CP      = 16
    I_PRECIP_PHI_MP      = 17
    I_PRECIP_KIN_CP      = 18
    I_PRECIP_KIN_MP      = 19
    !
    I_EVAP_ENERGY        = 20
    I_EVAP_EIN           = 21
    I_EVAP_LH            = 22
    I_EVAP_PHI           = 23
    I_EVAP_KIN           = 24
    !
    I_SFCRAD_ENERGY      = 25
    I_TOARAD_ENERGY      = 26
    !
    I_OLR                = 27
    I_EVAP_SFC           = 28
    !
    I_PRE_SFC            = 29
    I_TEM_SFC            = 30
    I_TH_SFC             = 31
    I_RHO_SFC            = 32
    !
    I_RCOSZ              = 33
    I_RFLUXS_SU          = 34
    I_RFLUXS_SD          = 35
    I_RFLUXS_LU          = 36
    I_RFLUXS_LD          = 37
    !
    I_RFLUX_TOA_SU       = 38
    I_RFLUX_TOA_SD       = 39
    I_RFLUX_TOA_LU       = 40
    I_RFLUX_TOA_LD       = 41
    !
    I_RFLUX_TOA_SU_C     = 42
    I_RFLUX_TOA_SD_C     = 43
    I_RFLUX_TOA_LU_C     = 44
    I_RFLUX_TOA_LD_C     = 45
    !
    I_SH_FLUX_SFC        = 46
    I_LH_FLUX_SFC        = 47
    I_TAUX_SFC           = 48
    I_TAUY_SFC           = 49
    I_TAUZ_SFC           = 50
    !
    I_VX10               = 51
    I_VY10               = 52
    I_VZ10               = 53
    I_T2                 = 54
    I_Q2                 = 55
    !
    I_CEV                = 56  !  07/07/23 K.Suzuki added for SPRINTARS
    !
!   if(RAIN_TYPE=='CLOUD_PARAM') then
    if(CP_TYPE/='NONE') then  ! 10/06/15 A.T.Noda
       I_CUMFRC          = 57
       I_ALL             = 57
    else
       I_ALL             = 56
    end if
    !  [Add] 12/03/28 T.Seiki
    if(AE_TYPE /= 'NONE' )then
       I_ALL       = I_ALL+1
       I_VFRICTION = I_ALL
    end if
    !
    ! 04/11/30 M.Satoh
    I_ALL = I_ALL + 1
    I_ALBEDO_SFC = I_ALL
    !
    allocate(KMAX(I_ALL))
    allocate(KSTR(I_ALL))
    allocate(KEND(I_ALL))
    !
    KMAX(:) = ADM_KNONE
    KMAX(I_PRECIP) = 2
    KMAX(I_PRECIP_CP) = 2
    KMAX(I_PRECIP_MP) = 2
    KMAX(I_ALBEDO_SFC) = NRDIR * NRBND
    !
    KSUM = 0
    do i = 1, I_ALL
       KSTR(i) = KSUM + 1
       KSUM = KSUM + KMAX(i)
       KEND(i) = KSUM
    end do
    !
    allocate(sfcvar(    &
         ADM_gall,      &
         KSUM,          &
         ADM_lall,      &
         DIAG_MAX))
    allocate(sfcvar_pl( &
         ADM_GALL_PL,   &
         KSUM,          &
         ADM_lall_pl,   &
         DIAG_MAX))
    !
    sfcvar(:,:,:,:)=0.0D0
    sfcvar_pl(:,:,:,:)=0.0D0
!    ! 2010.5.17. M.Satoh
!    sfcvar(:,:,:,:)=VMISS
!    sfcvar_pl(:,:,:,:)=VMISS
    !
    return
    !
  end subroutine sfcvar_setup
  !-----------------------------------------------------------------------------
  subroutine sfcvar_get( &
       sv, sv_pl,        &  !--- OUT : surface variable
       vid               &  !--- IN  : variable ID
       )
    !------
    !------ get prognostic variables from diag[num].
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_GALL_PL,       &
         ADM_LALL_PL,       &
         ADM_prc_me,        &
         ADM_prc_pl,        &
         ADM_KNONE,         &
         ADM_gall,          &
         ADM_lall,          &
         ADM_gall_1d
    !
    implicit none
    REAL(RP), intent(out) :: sv(ADM_gall,ADM_KNONE,ADM_lall)
    REAL(RP), intent(out) :: sv_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL)
    integer, intent(in)  :: vid
    !
    sv(1:ADM_gall,ADM_KNONE,1:ADM_lall) &
         = sfcvar(1:ADM_gall,KSTR(vid),1:ADM_lall,1)

    if(ADM_prc_me==ADM_prc_pl) then
       sv_pl(1:ADM_gall_pl,ADM_KNONE,1:ADM_lall_pl) &
            = sfcvar_pl(1:ADM_gall_pl,KSTR(vid),1:ADM_lall_pl,1)
    end if
    !
    return
    !
  end subroutine sfcvar_get
  !-----------------------------------------------------------------------------
  subroutine sfcvar_get1( &
       sv, sv_pl,        &  !--- OUT : surface variable
       vid,              &  !--- IN  : variable ID
       mdim              &  !--- IN  : dimension
       )
    !------
    !------ get prognostic variables from diag[num].
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_GALL_PL,       &
         ADM_LALL_PL,       &
         ADM_prc_me,        &
         ADM_prc_pl,        &
         ADM_KNONE,         &
         ADM_gall,          &
         ADM_lall,          &
         ADM_gall_1d
    use mod_comm, only :    &
         COMM_data_transfer
    !
    implicit none
    integer, intent(in) :: mdim
    REAL(RP), intent(out) :: sv(ADM_gall,ADM_KNONE,ADM_lall,mdim)
    REAL(RP), intent(out) :: sv_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL,mdim)
    integer, intent(in)  :: vid
    !
    integer :: m,l,n
    !
    do m=1, mdim
       do l=1,ADM_lall
          do n=1,ADM_gall
             sv(n,ADM_KNONE,l,m) = sfcvar(n,KSTR(vid)+m-1,l,1)
          end do
       end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
       do m=1, mdim
          do l=1,ADM_lall_pl
             do n=1,ADM_gall_pl
                sv_pl(n,ADM_KNONE,l,m) = sfcvar_pl(n,KSTR(vid)+m-1,l,1)
             end do
          end do
       end do
    end if
    !
    return
    !
  end subroutine sfcvar_get1
  !-----------------------------------------------------------------------------
  subroutine sfcvar_get2(  &
       sv, sv_pl,          &  !--- OUT : surface variable
       vid,                &  !--- IN : variable ID
       mdim1, mdim2        &  !--- IN
       )
    !------
    !------ get diagnostic variables
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_lall,          &
         ADM_gall,          &
         ADM_gall_pl,       &
         ADM_lall_pl,       &
         ADM_prc_pl,        &
         ADM_prc_me

    !
    implicit none
    integer, intent(in) :: mdim1, mdim2
    REAL(RP), intent(out) :: sv(ADM_gall,ADM_KNONE,ADM_lall,mdim1,mdim2)
    REAL(RP), intent(out) :: sv_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL,mdim1,mdim2)
    integer, intent(in) :: vid
    !
    integer :: l,n
    integer :: m1,m2, k
    !
    ! assume: vid == I_ALBEDO_SFC
    ! assume: NVAR(I_ALBEDO_SFC) = NRDIR * NRBND
    k = 0
    ! modify T.Mitsui 06.04.18
!!$    do m1 = 1, mdim1
!!$       do m2 = 1, mdim2
    do m2 = 1, mdim2
       do m1 = 1, mdim1
          k = k + 1
          do l=1, ADM_lall
             do n=1, ADM_gall
                sv(n,ADM_KNONE,l,m1,m2) = sfcvar(n,KSTR(vid)+k-1,l,1)
             end do
          end do
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       k = 0
       ! modify T.Mitsui 06.04.18
!!$       do m1 = 1, mdim1
!!$          do m2 = 1, mdim2
       do m2 = 1, mdim2
          do m1 = 1, mdim1
             k = k + 1
             do l=1, ADM_lall_pl
                do n=1, ADM_gall_pl
                   sv_pl(n,ADM_KNONE,l,m1,m2) = sfcvar_pl(n,KSTR(vid)+k-1,l,1)
                end do
             end do
          end do
       end do
    end if
    !
    return
    !
  end subroutine sfcvar_get2
  !-----------------------------------------------------------------------------
  subroutine sfcvar_set(   &
       sv, sv_pl,          &  !--- IN : surface variable
       vid                 &  !--- IN : variable ID
       )
    !------
    !------ set diagnostic variables
    !------ and COMMUNICATION.
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_GALL_PL,       &
         ADM_LALL_PL,       &
         ADM_prc_me,        &
         ADM_prc_pl,        &
         ADM_KNONE,         &
         ADM_gall,          &
         ADM_gmin,          &
         ADM_gmax,          &
         ADM_lall,          &
         ADM_gall_1d
    !
    implicit none
    REAL(RP), intent(in) :: sv(ADM_gall,ADM_KNONE,ADM_lall)
    REAL(RP), intent(in) :: sv_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL)
    integer, intent(in) :: vid
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    sfcvar(1:ADM_gall,KSTR(vid),1:ADM_lall,1) &
         = sv(1:ADM_gall,ADM_KNONE,1:ADM_lall)
    if(ADM_prc_me==ADM_prc_pl) then
       sfcvar_pl(1:ADM_gall_pl,KSTR(vid),1:ADM_lall_pl,1) &
            = sv_pl(1:ADM_gall_pl,ADM_KNONE,1:ADM_lall_pl)
    end if
    !
    sfcvar(suf(ADM_gmax+1,ADM_gmin-1),KSTR(vid),:,1) &
         = sfcvar(suf(ADM_gmax+1,ADM_gmin),KSTR(vid),:,1)
    sfcvar(suf(ADM_gmin-1,ADM_gmax+1),KSTR(vid),:,1) &
         = sfcvar(suf(ADM_gmin,ADM_gmax+1),KSTR(vid),:,1)
    !
    return
    !
  end subroutine sfcvar_set
  !-----------------------------------------------------------------------------
  subroutine sfcvar_get_in(&
       sv,                 &  !--- OUT : surface variable
       vid                 &  !--- IN : variable ID
       )
    !------
    !------ get diagnostic variables
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_lall,          &
         ADM_IopJop_nmax,   &
         ADM_IopJop,        &
         ADM_GIoJo
    !
    implicit none
    REAL(RP), intent(out) :: sv(ADM_IopJop_nmax,ADM_KNONE,ADM_lall)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    !
    do l=1, ADM_lall
       do n=1, ADM_IopJop_nmax
          nn = ADM_IopJop(n,ADM_GIoJo)
          sv(n,ADM_KNONE,l)  = sfcvar(nn,KSTR(vid),l,1)
       end do
    end do
    !
    return
    !
  end subroutine sfcvar_get_in
  !-----------------------------------------------------------------------------
  subroutine sfcvar_set_in( &
       sv,                  &  !--- IN : surface variable
       vid                  &  !--- IN : variable ID
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_gmin,          &
         ADM_gmax,          &
         ADM_lall,          &
         ADM_gall_1d,       &
         ADM_IopJop_nmax,   &
         ADM_IopJop,        &
         ADM_GIoJo
    !
    implicit none
    REAL(RP), intent(in) :: sv(ADM_IopJop_nmax,ADM_KNONE,ADM_lall)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    do l=1, ADM_lall
       do n=1, ADM_IopJop_nmax
          nn = ADM_IopJop(n,ADM_GIoJo)
          sfcvar(nn,KSTR(vid),l,1) = sv(n,ADM_KNONE,l)
       end do
    end do
    !
    sfcvar(suf(ADM_gmax+1,ADM_gmin-1),KSTR(vid),:,1) &
         = sfcvar(suf(ADM_gmax+1,ADM_gmin),KSTR(vid),:,1)
    sfcvar(suf(ADM_gmin-1,ADM_gmax+1),KSTR(vid),:,1) &
         = sfcvar(suf(ADM_gmin,ADM_gmax+1),KSTR(vid),:,1)
    !
  end subroutine sfcvar_set_in
  !-----------------------------------------------------------------------------
  subroutine sfcvar_get1_in(&
       sv,                 &  !--- OUT : surface variable
       vid,                &  !--- IN : variable ID
       mdim                &  !--- IN : dimension
       )
    !------
    !------ get diagnostic variables
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_lall,          &
         ADM_IopJop_nmax,   &
         ADM_IopJop,        &
         ADM_GIoJo
    !
    implicit none
    integer, intent(in) :: mdim
    REAL(RP), intent(out) :: sv(ADM_IopJop_nmax,ADM_KNONE,ADM_lall,mdim)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    integer :: m
    !
    do m=1, mdim
       do l=1, ADM_lall
          do n=1, ADM_IopJop_nmax
             nn = ADM_IopJop(n,ADM_GIoJo)
             sv(n,ADM_KNONE,l,m)  = sfcvar(nn,KSTR(vid)+m-1,l,1)
          end do
       end do
    end do
    !
    return
    !
  end subroutine sfcvar_get1_in
  !
  subroutine sfcvar_set1(   &
       sv, sv_pl,           &  !--- IN : surface variable
       vid,                 &  !--- IN : variable ID
       mdim,                &  !--- IN : dimension
       comm_flag            &  !--- IN, optional
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_gmin,          &
         ADM_gmax,          &
         ADM_gall,          &
         ADM_lall,          &
         ADM_gall_pl,       &
         ADM_lall_pl,       &
         ADM_gall_1d,       &
         adm_prc_me,        &
         adm_prc_pl
    use mod_comm, only :    &
         COMM_data_transfer
    !
    implicit none
    integer, intent(in) :: mdim
    REAL(RP), intent(in) :: sv(ADM_gall,ADM_KNONE,ADM_lall,mdim)
    REAL(RP), intent(in) :: sv_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,mdim)
    integer, intent(in) :: vid
    integer, optional, intent(in) :: comm_flag
    !
    integer :: l,n,nn
    integer :: m
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    integer :: ierr

    do m=1, mdim
       do l=1, ADM_lall
          do n=1, ADM_gall
             sfcvar(n,KSTR(vid)+m-1,l,1) = sv(n,ADM_KNONE,l,m)
          end do
       end do
    end do
    if(ADM_prc_me==ADM_prc_pl) then
       do m=1, mdim
          do l=1,ADM_lall_pl
             do n=1,ADM_gall_pl
                sfcvar_pl(n,KSTR(vid)+m-1,l,1) = sv_pl(n,ADM_KNONE,l,m)
             end do
          end do
       end do
    end if

    if ( present(comm_flag) ) then
       if(comm_flag==1) then
          call COMM_data_transfer(sfcvar(:,:,:,:),sfcvar_pl(:,:,:,:))
       endif

       sfcvar(suf(ADM_gmax+1,ADM_gmin-1),:,:,:) = sfcvar(suf(ADM_gmax+1,ADM_gmin),:,:,:)
       sfcvar(suf(ADM_gmin-1,ADM_gmax+1),:,:,:) = sfcvar(suf(ADM_gmin,ADM_gmax+1),:,:,:)
    endif

    return
  end subroutine sfcvar_set1
  !-----------------------------------------------------------------------------
  subroutine sfcvar_set1_in( &
       sv,                  &  !--- IN : surface variable
       vid,                 &  !--- IN : variable ID
       mdim                 &  !--- IN : dimension
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_gmin,          &
         ADM_gmax,          &
         ADM_lall,          &
         ADM_gall_1d,       &
         ADM_IopJop_nmax,   &
         ADM_IopJop,        &
         ADM_GIoJo
    !
    implicit none
    integer, intent(in) :: mdim
    REAL(RP), intent(in) :: sv(ADM_IopJop_nmax,ADM_KNONE,ADM_lall,mdim)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    integer :: m
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    do m=1, mdim
       do l=1, ADM_lall
          do n=1, ADM_IopJop_nmax
             nn = ADM_IopJop(n,ADM_GIoJo)
             sfcvar(nn,KSTR(vid)+m-1,l,1) = sv(n,ADM_KNONE,l,m)
          end do
       end do
    end do
    sfcvar(suf(ADM_gmax+1,ADM_gmin-1),KSTR(vid):KEND(vid),:,1) &
         = sfcvar(suf(ADM_gmax+1,ADM_gmin),KSTR(vid):KEND(vid),:,1)
    sfcvar(suf(ADM_gmin-1,ADM_gmax+1),KSTR(vid):KEND(vid),:,1) &
         = sfcvar(suf(ADM_gmin,ADM_gmax+1),KSTR(vid):KEND(vid),:,1)
    !
  end subroutine sfcvar_set1_in
  !-----------------------------------------------------------------------------
  subroutine sfcvar_get2_in(&
       sv,                 &  !--- OUT : surface variable
       vid,                &  !--- IN : variable ID
       mdim1, mdim2        &  !--- IN
       )
    !------
    !------ get diagnostic variables
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_lall,          &
         ADM_IopJop_nmax,   &
         ADM_IopJop,        &
         ADM_GIoJo
    !
    implicit none
    integer, intent(in) :: mdim1, mdim2
    REAL(RP), intent(out) :: sv(ADM_IopJop_nmax,ADM_KNONE,ADM_lall,mdim1,mdim2)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    integer :: m1,m2, k
    !
    ! assume: vid == I_ALBEDO_SFC
    ! assume: NVAR(I_ALBEDO_SFC) = NRDIR * NRBND
    k = 0
    ! modify T.Mitsui 06.04.18
!!$    do m1 = 1, mdim1
!!$       do m2 = 1, mdim2
    do m2 = 1, mdim2
       do m1 = 1, mdim1
          k = k + 1
          do l=1, ADM_lall
             do n=1, ADM_IopJop_nmax
                nn = ADM_IopJop(n,ADM_GIoJo)
                sv(n,ADM_KNONE,l,m1,m2) = sfcvar(nn,KSTR(vid)+k-1,l,1)
             end do
          end do
       end do
    end do
    !
    return
    !
  end subroutine sfcvar_get2_in
  !-----------------------------------------------------------------------------
  subroutine sfcvar_set2_in( &
       sv,                  &  !--- IN : surface variable
       vid,                 &  !--- IN : variable ID
       mdim1, mdim2        &  !--- IN
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_gmin,          &
         ADM_gmax,          &
         ADM_lall,          &
         ADM_gall_1d,       &
         ADM_IopJop_nmax,   &
         ADM_IopJop,        &
         ADM_GIoJo
    !
    implicit none
    integer, intent(in) :: mdim1, mdim2
    REAL(RP), intent(in) :: sv(ADM_IopJop_nmax,ADM_KNONE,ADM_lall,mdim1,mdim2)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    integer :: m1, m2, k
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    ! assume: vid == I_ALBEDO_SFC
    ! assume: NVAR(I_ALBEDO_SFC) = NRDIR * NRBND
    !
    k = 0
    ! modify T.Mitsui 06.04.18
!!$    do m1 = 1, mdim1
!!$       do m2 = 1, mdim2
    do m2 = 1, mdim2
       do m1 = 1, mdim1
          k = k + 1
          do l=1, ADM_lall
             do n=1, ADM_IopJop_nmax
                nn = ADM_IopJop(n,ADM_GIoJo)
                sfcvar(nn,KSTR(vid)+k-1,l,1) = sv(n,ADM_KNONE,l,m1,m2)
             end do
          end do
       end do
    end do
    sfcvar(suf(ADM_gmax+1,ADM_gmin-1),KSTR(vid):KEND(vid),:,1) &
         = sfcvar(suf(ADM_gmax+1,ADM_gmin),KSTR(vid):KEND(vid),:,1)
    sfcvar(suf(ADM_gmin-1,ADM_gmax+1),KSTR(vid):KEND(vid),:,1) &
         = sfcvar(suf(ADM_gmin,ADM_gmax+1),KSTR(vid):KEND(vid),:,1)
    !
  end subroutine sfcvar_set2_in
  !-----------------------------------------------------------------------------
  ! ADD T.Mitsui 06.04.18
  subroutine sfcvar_set2( &
       sv,sv_pl,          &  !--- IN : surface variable
       vid,               &  !--- IN : variable ID
       mdim1, mdim2       &  !--- IN
       )
    use mod_adm, only :     &
         ADM_KNONE,         &
         ADM_gall,          &
         ADM_gall_pl,       &
         ADM_lall,          &
         ADM_lall_pl,       &
         adm_prc_me,        &
         adm_prc_pl
    !
    implicit none
    integer, intent(in) :: mdim1, mdim2
    REAL(RP), intent(in) :: sv   (ADM_gall   ,ADM_KNONE,ADM_lall   ,mdim1,mdim2)
    REAL(RP), intent(in) :: sv_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,mdim1,mdim2)
    integer, intent(in) :: vid
    !
    integer :: l,ij
    integer :: m1, m2, k
    !
    ! assume: vid == I_ALBEDO_SFC
    ! assume: NVAR(I_ALBEDO_SFC) = NRDIR * NRBND
    !
    k = 0
    do m2 = 1, mdim2
       do m1 = 1, mdim1
          k = k + 1
          do l=1, ADM_lall
             do ij=1, adm_gall
                sfcvar(ij,KSTR(vid)+k-1,l,1) = sv(ij,ADM_KNONE,l,m1,m2)
             end do
          end do
       end do
    end do
    !
    if(adm_prc_me == adm_prc_pl) then
       k=0
       do m2=1,mdim2
          do m1=1,mdim1
             k = k + 1
             do l=1, adm_lall_pl
                do ij=1,ADM_gall_pl
                   sfcvar_pl(ij,KSTR(vid)+k-1,l,1)=sv_pl(ij,ADM_KNONE,l,m1,m2)
                end do
             end do
          end do
       end do
    end if
    !
    return
  end subroutine sfcvar_set2

  !-----------------------------------------------------------------------------
  subroutine sfcvar_comm
    use mod_comm, only: &
       COMM_var
    implicit none
    !---------------------------------------------------------------------------

    call COMM_var( sfcvar, sfcvar_pl, KSUM, DIAG_MAX )

    return
  end subroutine sfcvar_comm

end module mod_sfcvar
!-------------------------------------------------------------------------------
