!-------------------------------------------------------------------------------
!>
!! Vertical metrics module
!!
!! @par Description
!!         In this module, the vertical metrics is calculated for the
!!         non-hydrostatic icoshaedral model.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2006-08-11 (        )  Trivial bug fix for VMTR_VOLUME in using shallow water model
!!
!<
module mod_vmtr
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: VMTR_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  ! index for VMTR_C2Wfact
  integer, public, parameter :: I_a_GZXH = 1
  integer, public, parameter :: I_b_GZXH = 2
  integer, public, parameter :: I_a_GZYH = 3
  integer, public, parameter :: I_b_GZYH = 4
  integer, public, parameter :: I_a_GZZH = 5
  integer, public, parameter :: I_b_GZZH = 6

  !--- gsqrt X gamma^2 at the interger level
  real(8), public, allocatable, save :: VMTR_GSGAM2   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GSGAM2_pl(:,:,:)

  !--- 1/(gsqrt X gamma^2) at the interger level
  real(8), public, allocatable, save :: VMTR_RGSGAM2   (:,:,:)
  real(8), public, allocatable, save :: VMTR_RGSGAM2_pl(:,:,:)

  !--- gsqrt X gamma^2 at the half-integer level
  real(8), public, allocatable, save :: VMTR_GSGAM2H   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GSGAM2H_pl(:,:,:)

  !--- 1/(gsqrt X gamma^2) at the half-integer level
  real(8), public, allocatable, save :: VMTR_RGSGAM2H   (:,:,:)
  real(8), public, allocatable, save :: VMTR_RGSGAM2H_pl(:,:,:)

  !--- vector G^z at the half-integer level
  real(8), public, allocatable, save :: VMTR_GZXH   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GZXH_pl(:,:,:)
  real(8), public, allocatable, save :: VMTR_GZYH   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GZYH_pl(:,:,:)
  real(8), public, allocatable, save :: VMTR_GZZH   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GZZH_pl(:,:,:)

  !--- vector G^z at the integer level
  real(8), public, allocatable, save :: VMTR_GZX   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GZX_pl(:,:,:)
  real(8), public, allocatable, save :: VMTR_GZY   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GZY_pl(:,:,:)
  real(8), public, allocatable, save :: VMTR_GZZ   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GZZ_pl(:,:,:)

  !--- gsqrt X gamma at the half-integer level
  real(8), public, allocatable, save :: VMTR_GSGAMH   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GSGAMH_pl(:,:,:)

  !--- 1/gamma at the integer level
  real(8), public, allocatable, save :: VMTR_RGAM   (:,:,:)
  real(8), public, allocatable, save :: VMTR_RGAM_pl(:,:,:)

  !--- 1/gamma at the half-integer level
  real(8), public, allocatable, save :: VMTR_RGAMH   (:,:,:)
  real(8), public, allocatable, save :: VMTR_RGAMH_pl(:,:,:)

  !--- 1/gamma^2 at the integer level
  real(8), public, allocatable, save :: VMTR_RGAM2   (:,:,:)
  real(8), public, allocatable, save :: VMTR_RGAM2_pl(:,:,:)

  !--- gamma^2 at the half-integer level
  real(8), public, allocatable, save :: VMTR_GAM2H   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GAM2H_pl(:,:,:)

  !--- gamma^2 at the integer level
  real(8), public, allocatable, save :: VMTR_GAM2   (:,:,:)
  real(8), public, allocatable, save :: VMTR_GAM2_pl(:,:,:)

  !--- 1/gamma^2 at the half-integer level
  real(8), public, allocatable, save :: VMTR_RGAM2H   (:,:,:)
  real(8), public, allocatable, save :: VMTR_RGAM2H_pl(:,:,:)

  !--- 1/gsqrt at the half-integer level
  real(8), public, allocatable, save :: VMTR_RGSH   (:,:,:)
  real(8), public, allocatable, save :: VMTR_RGSH_pl(:,:,:)

  !--- volume at the integer level
  real(8), public, allocatable, save :: VMTR_VOLUME   (:,:,:)
  real(8), public, allocatable, save :: VMTR_VOLUME_pl(:,:,:)

  ! [Add] 20120717 H.Yashiro
  !--- geopotential at the integer level
  real(8), public, allocatable, save :: VMTR_PHI   (:,:,:)
  real(8), public, allocatable, save :: VMTR_PHI_pl(:,:,:)

  !--- factor for integer to half integer level
  real(8), public, allocatable, save :: VMTR_C2Wfact   (:,:,:,:)
  real(8), public, allocatable, save :: VMTR_C2Wfact_pl(:,:,:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
!  logical, private, save :: deep = .true. [mod] 20120704 H.Yashiro
  logical, private, save :: deep = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !>
  !> Setup the vertical metrics
  !>
  subroutine VMTR_setup
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_LOG_FID,   &
       ADM_CTL_FID,   &
       ADM_prc_me,    &
       ADM_prc_pl,    &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_kall,      &
       ADM_kmin,      &
       ADM_kmax,      &
       ADM_KNONE,     &
       ADM_gmin,      &
       ADM_gall_1d
    use mod_cnst, only: &
       CNST_ERADIUS, &
       CNST_EGRAV
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_Z,     &
       GRD_ZH,    &
       GRD_vz,    &
       GRD_vz_pl, &
       GRD_dgz,   &
       GRD_dgzh,  &
       GRD_afac,  &
       GRD_bfac,  &
       GRD_grid_type
    use mod_gmtr, only: &
       GMTR_P_AREA,  &
       GMTR_P_var,   &
       GMTR_P_var_pl
!    use mod_oprt_plane, only: &
!       OPRT_PLANE_gradient
    use mod_oprt, only: &
       OPRT_gradient,         &
       OPRT_horizontalize_vec
    implicit none

    integer, parameter :: var_max = 6

    integer, parameter :: JXH     = 1
    integer, parameter :: JYH     = 2
    integer, parameter :: JZH     = 3
    integer, parameter :: JX      = 4
    integer, parameter :: JY      = 5
    integer, parameter :: JZ      = 6

    real(8) :: var   (ADM_gall,   ADM_kall,ADM_lall,   var_max)
    real(8) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,var_max)

    real(8) :: GSQRT    (ADM_gall,   ADM_kall)
    real(8) :: GSQRT_pl (ADM_gall_pl,ADM_kall)
    real(8) :: GSQRTH   (ADM_gall,   ADM_kall)
    real(8) :: GSQRTH_pl(ADM_gall_pl,ADM_kall)

    real(8) :: GAM, GAMH

    namelist / VMTRPARAM / &
        deep

    integer :: ierr
    integer :: n, k, l
    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i) 
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[vmtr]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=VMTRPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** VMTRPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist VMTRPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist VMTRPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,VMTRPARAM)

    !--- initialization
    allocate( VMTR_GSGAM2     (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGSGAM2    (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GSGAM2H    (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2H_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGSGAM2H   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2H_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GSGAMH     (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAMH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GZXH       (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GZXH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GZYH       (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GZYH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GZZH       (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GZZH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GZX        (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GZX_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GZY        (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GZY_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GZZ        (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GZZ_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGAM       (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAM_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGAMH      (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAMH_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GAM2       (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GAM2_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGAM2      (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAM2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_GAM2H      (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GAM2H_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGAM2H     (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAM2H_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGSH       (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_VOLUME     (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_VOLUME_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_PHI        (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_PHI_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_C2Wfact   (6,ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_C2Wfact_pl(6,ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    !--- if 1 layer model( shallow water model ),
    if ( ADM_kall == ADM_KNONE ) then

       VMTR_VOLUME   (:,:,:) = GMTR_P_var   (:,:,:,GMTR_P_AREA)
       VMTR_VOLUME_pl(:,:,:) = GMTR_P_var_pl(:,:,:,GMTR_P_AREA)

       return
    endif

    var   (:,:,:,:) = 0.D0
    var_pl(:,:,:,:) = 0.D0

    if ( trim(GRD_grid_type) == 'ON_PLANE' ) then
!       --- calculation of Jxh, Jyh, and Jzh
!       call OPRT_PLANE_gradient( var   (:,:,:,JXH),   & !--- [OUT]
!                                 var   (:,:,:,JYH),   & !--- [OUT]
!                                 GRD_vz(:,:,:,GRD_ZH) ) !--- [IN]
!       --- calculation of Jx, Jy, and Jz
!       call OPRT_PLANE_gradient( var   (:,:,:,JX),   & !--- [OUT]
!                                 var   (:,:,:,JY),   & !--- [OUT]
!                                 GRD_vz(:,:,:,GRD_Z) ) !--- [IN]
    else
       !--- calculation of Jxh, Jyh, and Jzh
       call OPRT_gradient( var(:,:,:,JXH),       var_pl(:,:,:,JXH),      & !--- [OUT]
                           var(:,:,:,JYH),       var_pl(:,:,:,JYH),      & !--- [OUT]
                           var(:,:,:,JZH),       var_pl(:,:,:,JZH),      & !--- [OUT]
                           GRD_vz(:,:,:,GRD_ZH), GRD_vz_pl(:,:,:,GRD_ZH) ) !--- [IN]

       call OPRT_horizontalize_vec( var(:,:,:,JXH), var_pl(:,:,:,JXH), & !--- [INOUT]
                                    var(:,:,:,JYH), var_pl(:,:,:,JYH), & !--- [INOUT]
                                    var(:,:,:,JZH), var_pl(:,:,:,JZH)  ) !--- [INOUT]
       !--- calculation of Jx, Jy, and Jz
       call OPRT_gradient( var(:,:,:,JX),       var_pl(:,:,:,JX),      & !--- [OUT]
                           var(:,:,:,JY),       var_pl(:,:,:,JY),      & !--- [OUT]
                           var(:,:,:,JZ),       var_pl(:,:,:,JZ),      & !--- [OUT]
                           GRD_vz(:,:,:,GRD_Z), GRD_vz_pl(:,:,:,GRD_Z) ) !--- [IN]

       call OPRT_horizontalize_vec( var(:,:,:,JX), var_pl(:,:,:,JX), & !--- [INOUT]
                                    var(:,:,:,JY), var_pl(:,:,:,JY), & !--- [INOUT]
                                    var(:,:,:,JZ), var_pl(:,:,:,JZ)  ) !--- [INOUT]
    endif

    !--- fill HALO
    call COMM_data_transfer( var, var_pl )
    var(suf(1,ADM_gall_1d),:,:,:) = var(suf(ADM_gmin,ADM_gmin),:,:,:)
    var(suf(ADM_gall_1d,1),:,:,:) = var(suf(ADM_gmin,ADM_gmin),:,:,:)

    do l = 1, ADM_lall

       !--- calculation of G^{1/2} at integer points
       !---   G^{1/2} = dz/dgz
       do k = ADM_kmin, ADM_kmax
       do n = 1, ADM_gall
          GSQRT(n,k) = ( GRD_vz(n,k+1,l,GRD_ZH) - GRD_vz(n,k,l,GRD_ZH) ) / GRD_dgz(k)
       enddo
       enddo
       do n = 1, ADM_gall
          GSQRT(n,ADM_kmin-1) = GSQRT(n,ADM_kmin)
          GSQRT(n,ADM_kmax+1) = GSQRT(n,ADM_kmax)
       enddo

       !--- calculation of G^{1/2} at half points
       do k = ADM_kmin, ADM_kmax+1
       do n = 1, ADM_gall
          GSQRTH(n,k) = ( GRD_vz(n,k,l,GRD_Z) - GRD_vz(n,k-1,l,GRD_Z) ) / GRD_dgzh(k)
       enddo
       enddo
       do n = 1, ADM_gall
          GSQRTH(n,ADM_kmin-1) = GSQRTH(n,ADM_kmin)
       enddo

       if ( deep ) then
          do k = 1, ADM_kall
          do n = 1, ADM_gall
             !--- calculation gamma at integer/half points
             GAM  = 1.D0 + GRD_vz(n,k,l,GRD_Z)  / CNST_ERADIUS
             GAMH = 1.D0 + GRD_vz(n,k,l,GRD_ZH) / CNST_ERADIUS

             VMTR_GSGAMH  (n,k,l) = GSQRTH(n,k) * GAMH

             VMTR_GAM2    (n,k,l) = GAM  * GAM
             VMTR_GAM2H   (n,k,l) = GAMH * GAMH

             VMTR_GSGAM2  (n,k,l) = GSQRT (n,k) * GAM  * GAM
             VMTR_GSGAM2H (n,k,l) = GSQRTH(n,k) * GAMH * GAMH

             VMTR_RGSH    (n,k,l) = 1.D0 / GSQRTH(n,k)
             VMTR_RGAM    (n,k,l) = 1.D0 / GAM
             VMTR_RGAMH   (n,k,l) = 1.D0 / GAMH
             VMTR_RGAM2   (n,k,l) = 1.D0 / VMTR_GAM2   (n,k,l)
             VMTR_RGAM2H  (n,k,l) = 1.D0 / VMTR_GAM2H  (n,k,l)
             VMTR_RGSGAM2 (n,k,l) = 1.D0 / VMTR_GSGAM2 (n,k,l)
             VMTR_RGSGAM2H(n,k,l) = 1.D0 / VMTR_GSGAM2H(n,k,l)
          enddo
          enddo
       else
          do k = 1, ADM_kall
          do n = 1, ADM_gall
             ! GAM  = 1.D0
             ! GAMH = 1.D0

             VMTR_GSGAMH  (n,k,l) = GSQRTH(n,k)

             VMTR_GAM2    (n,k,l) = 1.D0
             VMTR_GAM2H   (n,k,l) = 1.D0

             VMTR_GSGAM2  (n,k,l) = GSQRT (n,k)
             VMTR_GSGAM2H (n,k,l) = GSQRTH(n,k)

             VMTR_RGSH    (n,k,l) = 1.D0 / GSQRTH(n,k)
             VMTR_RGAM    (n,k,l) = 1.D0
             VMTR_RGAMH   (n,k,l) = 1.D0
             VMTR_RGAM2   (n,k,l) = 1.D0
             VMTR_RGAM2H  (n,k,l) = 1.D0
             VMTR_RGSGAM2 (n,k,l) = 1.D0 / GSQRT (n,k)
             VMTR_RGSGAM2H(n,k,l) = 1.D0 / GSQRTH(n,k)
          enddo
          enddo
       endif

       do k = 1, ADM_kall
       do n = 1, ADM_gall
          !--- calculation of GZ at half points
          !---    GZX = - JX / G^{1/2}
          !---    GZY = - JY / G^{1/2}
          !---    GZZ = - JZ / G^{1/2}
          VMTR_GZXH(n,k,l) = -var(n,k,l,JXH) / GSQRTH(n,k)
          VMTR_GZYH(n,k,l) = -var(n,k,l,JYH) / GSQRTH(n,k)
          VMTR_GZZH(n,k,l) = -var(n,k,l,JZH) / GSQRTH(n,k)
          !--- calculation of GZ at full points
          VMTR_GZX (n,k,l) = -var(n,k,l,JX)  / GSQRT (n,k)
          VMTR_GZY (n,k,l) = -var(n,k,l,JY)  / GSQRT (n,k)
          VMTR_GZZ (n,k,l) = -var(n,k,l,JZ)  / GSQRT (n,k)

          !--- calculation of volume
          VMTR_VOLUME(n,k,l) = GMTR_P_var(n,ADM_KNONE,l,GMTR_P_AREA) &
                             * VMTR_GSGAM2(n,k,l)                    &
                             * GRD_dgz(k)

          !--- calculation of geopotential
          VMTR_PHI(n,k,l) = GRD_vz(n,k,l,GRD_Z) * CNST_EGRAV
       enddo
       enddo

       do k = ADM_kmin, ADM_kall
       do n = 1, ADM_gall
          !--- calculation of factor for integer to half integer level with Gz
          VMTR_C2Wfact(I_a_GZXH,n,k,l) = 0.5D0 * GRD_afac(k) / VMTR_GSGAM2(n,k  ,l) * VMTR_GSGAM2H(n,k,l) &
                                       * VMTR_GZXH(n,k,l)
          VMTR_C2Wfact(I_b_GZXH,n,k,l) = 0.5D0 * GRD_bfac(k) / VMTR_GSGAM2(n,k-1,l) * VMTR_GSGAM2H(n,k,l) &
                                       * VMTR_GZXH(n,k,l)
          VMTR_C2Wfact(I_a_GZYH,n,k,l) = 0.5D0 * GRD_afac(k) / VMTR_GSGAM2(n,k  ,l) * VMTR_GSGAM2H(n,k,l) &
                                       * VMTR_GZYH(n,k,l)
          VMTR_C2Wfact(I_b_GZYH,n,k,l) = 0.5D0 * GRD_bfac(k) / VMTR_GSGAM2(n,k-1,l) * VMTR_GSGAM2H(n,k,l) &
                                       * VMTR_GZYH(n,k,l)
          VMTR_C2Wfact(I_a_GZZH,n,k,l) = 0.5D0 * GRD_afac(k) / VMTR_GSGAM2(n,k  ,l) * VMTR_GSGAM2H(n,k,l) &
                                       * VMTR_GZZH(n,k,l)
          VMTR_C2Wfact(I_b_GZZH,n,k,l) = 0.5D0 * GRD_bfac(k) / VMTR_GSGAM2(n,k-1,l) * VMTR_GSGAM2H(n,k,l) &
                                       * VMTR_GZZH(n,k,l)
       enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

       !--- calculation of G^{1/2} at integer points
       !---   G^{1/2} = dz/dgz
       do k = ADM_kmin, ADM_kmax
       do n = 1, ADM_gall_pl
          GSQRT_pl(n,k) = ( GRD_vz_pl(n,k+1,l,GRD_ZH) - GRD_vz_pl(n,k,l,GRD_ZH) ) / GRD_dgz(k)
       enddo
       enddo
       do n = 1, ADM_gall_pl
          GSQRT_pl(n,ADM_kmin-1) = GSQRT_pl(n,ADM_kmin)
          GSQRT_pl(n,ADM_kmax+1) = GSQRT_pl(n,ADM_kmax)
       enddo

       !--- calculation of G^{1/2} at half points
       do k = ADM_kmin, ADM_kmax+1
       do n = 1, ADM_gall_pl
          GSQRTH_pl(n,k) = ( GRD_vz_pl(n,k,l,GRD_Z) - GRD_vz_pl(n,k-1,l,GRD_Z) ) / GRD_dgzh(k)
       enddo
       enddo
       do n = 1, ADM_gall_pl
          GSQRTH_pl(n,ADM_kmin-1) = GSQRTH_pl(n,ADM_kmin)
       enddo

       if ( deep ) then
          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl
             !--- calculation gamma at integer/half points
             GAM  = 1.D0 + GRD_vz_pl(n,k,l,GRD_Z)  / CNST_ERADIUS
             GAMH = 1.D0 + GRD_vz_pl(n,k,l,GRD_ZH) / CNST_ERADIUS

             VMTR_GSGAMH_pl  (n,k,l) = GSQRTH(n,k) * GAMH

             VMTR_GAM2_pl    (n,k,l) = GAM  * GAM
             VMTR_GAM2H_pl   (n,k,l) = GAMH * GAMH

             VMTR_GSGAM2_pl  (n,k,l) = GSQRT_pl (n,k) * GAM  * GAM
             VMTR_GSGAM2H_pl (n,k,l) = GSQRTH_pl(n,k) * GAMH * GAMH

             VMTR_RGSH_pl    (n,k,l) = 1.D0 / GSQRTH_pl(n,k)
             VMTR_RGAM_pl    (n,k,l) = 1.D0 / GAM
             VMTR_RGAMH_pl   (n,k,l) = 1.D0 / GAMH
             VMTR_RGAM2_pl   (n,k,l) = 1.D0 / VMTR_GAM2_pl (n,k,l)
             VMTR_RGAM2H_pl  (n,k,l) = 1.D0 / VMTR_GAM2H_pl(n,k,l)
             VMTR_RGSGAM2_pl (n,k,l) = 1.D0 / VMTR_GSGAM2_pl (n,k,l)
             VMTR_RGSGAM2H_pl(n,k,l) = 1.D0 / VMTR_GSGAM2H_pl(n,k,l)
          enddo
          enddo
       else
          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl
             ! GAM  = 1.D0
             ! GAMH = 1.D0

             VMTR_GSGAMH_pl  (n,k,l) = GSQRTH_pl(n,k)

             VMTR_GAM2_pl    (n,k,l) = 1.D0
             VMTR_GAM2H_pl   (n,k,l) = 1.D0

             VMTR_GSGAM2_pl  (n,k,l) = GSQRT_pl (n,k)
             VMTR_GSGAM2H_pl (n,k,l) = GSQRTH_pl(n,k)

             VMTR_RGSH_pl    (n,k,l) = 1.D0 / GSQRTH_pl(n,k)
             VMTR_RGAM_pl    (n,k,l) = 1.D0
             VMTR_RGAMH_pl   (n,k,l) = 1.D0
             VMTR_RGAM2_pl   (n,k,l) = 1.D0
             VMTR_RGAM2H_pl  (n,k,l) = 1.D0
             VMTR_RGSGAM2_pl (n,k,l) = 1.D0 / GSQRT_pl (n,k)
             VMTR_RGSGAM2H_pl(n,k,l) = 1.D0 / GSQRTH_pl(n,k)
          enddo
          enddo
       endif

       do k = 1, ADM_kall
       do n = 1, ADM_gall_pl
          !--- calculation of GZ at half points
          !---    GZX = - JX / G^{1/2}
          !---    GZY = - JY / G^{1/2}
          !---    GZZ = - JZ / G^{1/2}
          VMTR_GZXH_pl(n,k,l) = -var_pl(n,k,l,JXH) / GSQRTH_pl(n,k)
          VMTR_GZYH_pl(n,k,l) = -var_pl(n,k,l,JYH) / GSQRTH_pl(n,k)
          VMTR_GZZH_pl(n,k,l) = -var_pl(n,k,l,JZH) / GSQRTH_pl(n,k)
          !--- calculation of GZ at full points
          VMTR_GZX_pl (n,k,l) = -var_pl(n,k,l,JX)  / GSQRT_pl (n,k)
          VMTR_GZY_pl (n,k,l) = -var_pl(n,k,l,JY)  / GSQRT_pl (n,k)
          VMTR_GZZ_pl (n,k,l) = -var_pl(n,k,l,JZ)  / GSQRT_pl (n,k)

          !--- calculation of volume
          VMTR_VOLUME_pl(n,k,l) = GMTR_P_var_pl(n,ADM_KNONE,l,GMTR_P_AREA) &
                                * VMTR_GSGAM2_pl(n,k,l)                    &
                                * GRD_dgz(k)

          !--- calculation of geopotential
          VMTR_PHI_pl(n,k,l) = GRD_vz_pl(n,k,l,GRD_Z) * CNST_EGRAV
       enddo
       enddo

       do k = ADM_kmin, ADM_kall
       do n = 1, ADM_gall_pl
          VMTR_C2Wfact_pl(I_a_GZXH,n,k,l) = 0.5D0 * GRD_afac(k) &
                                          * VMTR_GSGAM2H_pl(n,k,l) / VMTR_GSGAM2_pl(n,k  ,l) * VMTR_GZXH_pl(n,k,l)
          VMTR_C2Wfact_pl(I_b_GZXH,n,k,l) = 0.5D0 * GRD_bfac(k) &
                                          * VMTR_GSGAM2H_pl(n,k,l) / VMTR_GSGAM2_pl(n,k-1,l) * VMTR_GZXH_pl(n,k,l)
          VMTR_C2Wfact_pl(I_a_GZYH,n,k,l) = 0.5D0 * GRD_afac(k) &
                                          * VMTR_GSGAM2H_pl(n,k,l) / VMTR_GSGAM2_pl(n,k  ,l) * VMTR_GZYH_pl(n,k,l)
          VMTR_C2Wfact_pl(I_b_GZYH,n,k,l) = 0.5D0 * GRD_bfac(k) &
                                          * VMTR_GSGAM2H_pl(n,k,l) / VMTR_GSGAM2_pl(n,k-1,l) * VMTR_GZYH_pl(n,k,l)
          VMTR_C2Wfact_pl(I_a_GZZH,n,k,l) = 0.5D0 * GRD_afac(k) &
                                          * VMTR_GSGAM2H_pl(n,k,l) / VMTR_GSGAM2_pl(n,k  ,l) * VMTR_GZZH_pl(n,k,l)
          VMTR_C2Wfact_pl(I_b_GZZH,n,k,l) = 0.5D0 * GRD_bfac(k) &
                                          * VMTR_GSGAM2H_pl(n,k,l) / VMTR_GSGAM2_pl(n,k-1,l) * VMTR_GZZH_pl(n,k,l)
       enddo
       enddo

       enddo
    endif

    return
  end subroutine VMTR_setup

end module mod_vmtr
!-------------------------------------------------------------------------------
