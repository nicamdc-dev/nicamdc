!-------------------------------------------------------------------------------
!> Module vertical metrics
!!
!! @par Description
!!         Vertical metrics on the icosahedral grid system
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_vmtr
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof

  use mod_adm, only: &
     ADM_lall,    &
     ADM_lall_pl, &
     ADM_gall,    &
     ADM_gall_pl, &
     ADM_kall
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: VMTR_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: I_a = 1      ! index for W2Cfact
  integer, public, parameter :: I_b = 2

  integer, public, parameter :: I_c = 1      ! index for C2Wfact
  integer, public, parameter :: I_d = 2

  integer, public, parameter :: I_a_GZXH = 1 ! index for C2WfactGz
  integer, public, parameter :: I_b_GZXH = 2
  integer, public, parameter :: I_a_GZYH = 3
  integer, public, parameter :: I_b_GZYH = 4
  integer, public, parameter :: I_a_GZZH = 5
  integer, public, parameter :: I_b_GZZH = 6

#ifdef _FIXEDINDEX_
  real(RP), public              :: VMTR_GAM2H       (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_GAM2H_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_GSGAM2      (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_GSGAM2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_GSGAM2H     (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_GSGAM2H_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

  real(RP), public              :: VMTR_RGSQRTH     (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGSQRTH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGAM        (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGAM_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGAMH       (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGAMH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGSGAM2     (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGSGAM2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGSGAM2H    (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGSGAM2H_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

  real(RP), public              :: VMTR_W2Cfact     (ADM_gall   ,ADM_kall,2,ADM_lall   )
  real(RP), public              :: VMTR_W2Cfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl)
  real(RP), public              :: VMTR_C2Wfact     (ADM_gall   ,ADM_kall,2,ADM_lall   )
  real(RP), public              :: VMTR_C2Wfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl)
  real(RP), public              :: VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   )
  real(RP), public              :: VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl)

  real(RP), public              :: VMTR_VOLUME      (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_VOLUME_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

  real(RP), public              :: VMTR_PHI         (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_PHI_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
#else
  real(RP), public, allocatable :: VMTR_GAM2H       (:,:,:)   ! Gamma^2 at the half level
  real(RP), public, allocatable :: VMTR_GAM2H_pl    (:,:,:)
  real(RP), public, allocatable :: VMTR_GSGAM2      (:,:,:)   ! G^1/2 X Gamma^2 at the full level
  real(RP), public, allocatable :: VMTR_GSGAM2_pl   (:,:,:)
  real(RP), public, allocatable :: VMTR_GSGAM2H     (:,:,:)   ! G^1/2 X Gamma^2 at the half level
  real(RP), public, allocatable :: VMTR_GSGAM2H_pl  (:,:,:)

  real(RP), public, allocatable :: VMTR_RGSQRTH     (:,:,:)   ! 1 / G^1/2 at the half level
  real(RP), public, allocatable :: VMTR_RGSQRTH_pl  (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAM        (:,:,:)   ! 1 / Gamma at the integer level
  real(RP), public, allocatable :: VMTR_RGAM_pl     (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAMH       (:,:,:)   ! 1 / Gamma at the half level
  real(RP), public, allocatable :: VMTR_RGAMH_pl    (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2     (:,:,:)   ! 1 / (G^1/2 X Gamma^2) at the full level
  real(RP), public, allocatable :: VMTR_RGSGAM2_pl  (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2H    (:,:,:)   ! 1 / (G^1/2 X Gamma^2) at the half level
  real(RP), public, allocatable :: VMTR_RGSGAM2H_pl (:,:,:)

  real(RP), public, allocatable :: VMTR_W2Cfact     (:,:,:,:) ! factor for half to full level
  real(RP), public, allocatable :: VMTR_W2Cfact_pl  (:,:,:,:)
  real(RP), public, allocatable :: VMTR_C2Wfact     (:,:,:,:) ! factor for full to half level
  real(RP), public, allocatable :: VMTR_C2Wfact_pl  (:,:,:,:)
  real(RP), public, allocatable :: VMTR_C2WfactGz   (:,:,:,:) ! factor for full to half level with Gz
  real(RP), public, allocatable :: VMTR_C2WfactGz_pl(:,:,:,:)

  real(RP), public, allocatable :: VMTR_VOLUME      (:,:,:)   ! volume at the full level
  real(RP), public, allocatable :: VMTR_VOLUME_pl   (:,:,:)

  real(RP), public, allocatable :: VMTR_PHI         (:,:,:)   ! geopotential at the full level
  real(RP), public, allocatable :: VMTR_PHI_pl      (:,:,:)
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: VMTR_deep_atmos = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine VMTR_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_kmin,    &
       ADM_kmax
    use mod_const, only: &
       GRAV => CONST_GRAV
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_Z,        &
       GRD_ZH,       &
       GRD_dgz,      &
       GRD_dgzh,     &
       GRD_vz,       &
       GRD_vz_pl,    &
       GRD_afact,    &
       GRD_bfact,    &
       GRD_cfact,    &
       GRD_dfact,    &
       GRD_rscale
    use mod_gmtr, only: &
       GMTR_area,   &
       GMTR_area_pl
    use mod_oprt, only: &
       OPRT_gradient,          &
       OPRT_horizontalize_vec, &
       OPRT_coef_grad,         &
       OPRT_coef_grad_pl
    implicit none

    integer, parameter :: var_max = 6

    integer, parameter :: JXH     = 1
    integer, parameter :: JYH     = 2
    integer, parameter :: JZH     = 3
    integer, parameter :: JX      = 4
    integer, parameter :: JY      = 5
    integer, parameter :: JZ      = 6

    real(RP) :: var   (ADM_gall   ,ADM_kall,ADM_lall   ,var_max)
    real(RP) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,var_max)

    !--- G^1/2
    real(RP) :: GSQRT    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GSQRT_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GSQRTH   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GSQRTH_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- Gamma factor
    real(RP) :: GAM      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GAM_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GAMH     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GAMH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- vector G^z at the full level
    real(RP) :: GZX      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZX_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GZY      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZY_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GZZ      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZZ_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- vector G^z at the half level
    real(RP) :: GZXH     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZXH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GZYH     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZYH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GZZH     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZZH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    namelist / VMTRPARAM / &
       VMTR_deep_atmos

    integer  :: ierr
    integer  :: g, k, l
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[vmtr]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=VMTRPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** VMTRPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist VMTRPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=VMTRPARAM)

#ifndef _FIXEDINDEX_
    allocate( VMTR_GAM2H       (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_GAM2H_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_GSGAM2      (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_GSGAM2H     (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2H_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGSQRTH     (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSQRTH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGAM        (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAM_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGAMH       (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAMH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGSGAM2     (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGSGAM2H    (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2H_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_W2Cfact     (ADM_gall   ,ADM_kall,2,ADM_lall   ) )
    allocate( VMTR_W2Cfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl) )
    allocate( VMTR_C2Wfact     (ADM_gall   ,ADM_kall,2,ADM_lall   ) )
    allocate( VMTR_C2Wfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl) )
    allocate( VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   ) )
    allocate( VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl) )

    allocate( VMTR_VOLUME      (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_VOLUME_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_PHI         (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_PHI_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
#endif

    if( IO_L ) write(IO_FID_LOG,*) '*** setup metrics for 3-D control volume'

    !--- if 1 layer model( shallow water model ),
    if ( ADM_kall == ADM_KNONE ) then

       do k = 1, ADM_kall
          VMTR_VOLUME   (:,k,:) = GMTR_area   (:,:)
          VMTR_VOLUME_pl(:,k,:) = GMTR_area_pl(:,:)
       enddo

       return
    endif

    var   (:,:,:,:) = 0.0_RP
    var_pl(:,:,:,:) = 0.0_RP

    !--- calculation of Jxh, Jyh, and Jzh
    call OPRT_gradient( var           (:,:,:,JXH:JZH), var_pl           (:,:,:,JXH:JZH), & ! [OUT]
                        GRD_vz        (:,:,:,GRD_ZH),  GRD_vz_pl        (:,:,:,GRD_ZH),  & ! [IN]
                        OPRT_coef_grad(:,:,:,:),       OPRT_coef_grad_pl(:,:,:)          ) ! [IN]

    call OPRT_horizontalize_vec( var(:,:,:,JXH), var_pl(:,:,:,JXH), & ! [INOUT]
                                 var(:,:,:,JYH), var_pl(:,:,:,JYH), & ! [INOUT]
                                 var(:,:,:,JZH), var_pl(:,:,:,JZH)  ) ! [INOUT]

    !--- calculation of Jx, Jy, and Jz
    call OPRT_gradient( var           (:,:,:,JX:JZ), var_pl           (:,:,:,JX:JZ), & ! [OUT]
                        GRD_vz        (:,:,:,GRD_Z), GRD_vz_pl        (:,:,:,GRD_Z), & ! [IN]
                        OPRT_coef_grad(:,:,:,:),     OPRT_coef_grad_pl(:,:,:)        ) ! [IN]

    call OPRT_horizontalize_vec( var(:,:,:,JX), var_pl(:,:,:,JX), & ! [INOUT]
                                 var(:,:,:,JY), var_pl(:,:,:,JY), & ! [INOUT]
                                 var(:,:,:,JZ), var_pl(:,:,:,JZ)  ) ! [INOUT]

    !--- fill HALO
    call COMM_data_transfer( var, var_pl )



    !--- G^1/2 = dz/dgz
    do l = 1, ADM_lall
       !--- calculation of G^1/2 at full level
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          GSQRT(g,k,l) = ( GRD_vz(g,k+1,l,GRD_ZH) - GRD_vz(g,k,l,GRD_ZH) ) / GRD_dgz(k)
       enddo
       enddo
       do g = 1, ADM_gall
          GSQRT(g,ADM_kmin-1,l) = GSQRT(g,ADM_kmin,l)
          GSQRT(g,ADM_kmax+1,l) = GSQRT(g,ADM_kmax,l)
       enddo

       !--- calculation of G^1/2 at half level
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          GSQRTH(g,k,l) = ( GRD_vz(g,k,l,GRD_Z) - GRD_vz(g,k-1,l,GRD_Z) ) / GRD_dgzh(k)
       enddo
       enddo
       do g = 1, ADM_gall
          GSQRTH(g,ADM_kmin-1,l) = GSQRTH(g,ADM_kmin,l)
       enddo
    enddo

    !--- Gamma = (a+z) / a
    if ( VMTR_deep_atmos ) then
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          GAM (g,k,l) = 1.0_RP + GRD_vz(g,k,l,GRD_Z)  / GRD_rscale
          GAMH(g,k,l) = 1.0_RP + GRD_vz(g,k,l,GRD_ZH) / GRD_rscale
       enddo
       enddo
       enddo
    else
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          GAM (g,k,l) = 1.0_RP
          GAMH(g,k,l) = 1.0_RP
       enddo
       enddo
       enddo
    endif

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       VMTR_GAM2H   (g,k,l) = GAMH(g,k,l) * GAMH(g,k,l)
       VMTR_GSGAM2  (g,k,l) = GAM (g,k,l) * GAM (g,k,l) * GSQRT (g,k,l)
       VMTR_GSGAM2H (g,k,l) = GAMH(g,k,l) * GAMH(g,k,l) * GSQRTH(g,k,l)

       VMTR_RGSQRTH (g,k,l) = 1.0_RP / GSQRTH(g,k,l)
       VMTR_RGAM    (g,k,l) = 1.0_RP / GAM (g,k,l)
       VMTR_RGAMH   (g,k,l) = 1.0_RP / GAMH(g,k,l)
       VMTR_RGSGAM2 (g,k,l) = 1.0_RP / VMTR_GSGAM2 (g,k,l)
       VMTR_RGSGAM2H(g,k,l) = 1.0_RP / VMTR_GSGAM2H(g,k,l)
    enddo
    enddo
    enddo

    ! full level <-> half level interpolation factor
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          VMTR_C2Wfact(g,k,I_a,l) = GRD_afact(k) * VMTR_RGSGAM2(g,k  ,l) * VMTR_GSGAM2H(g,k,l)
          VMTR_C2Wfact(g,k,I_b,l) = GRD_bfact(k) * VMTR_RGSGAM2(g,k-1,l) * VMTR_GSGAM2H(g,k,l)
       enddo
       enddo
       do g = 1, ADM_gall
          VMTR_C2Wfact(g,ADM_kmin-1,I_a,l) = 0.0_RP
          VMTR_C2Wfact(g,ADM_kmin-1,I_b,l) = 0.0_RP
       enddo

       do k = ADM_kmin-1, ADM_kmax
       do g = 1, ADM_gall
          VMTR_W2Cfact(g,k,I_c,l) = GRD_cfact(k) * VMTR_GSGAM2(g,k,l) * VMTR_RGSGAM2H(g,k+1,l)
          VMTR_W2Cfact(g,k,I_d,l) = GRD_dfact(k) * VMTR_GSGAM2(g,k,l) * VMTR_RGSGAM2H(g,k  ,l)
       enddo
       enddo
       do g = 1, ADM_gall
          VMTR_W2Cfact(g,ADM_kmax+1,I_c,l) = 0.0_RP
          VMTR_W2Cfact(g,ADM_kmax+1,I_d,l) = 0.0_RP
       enddo
    enddo

    ! full level <-> half level interpolation factor with Gz

    !--- Gz(X) = - JX / G^1/2
    !--- Gz(Y) = - JY / G^1/2
    !--- Gz(Z) = - JZ / G^1/2
    do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          GZX (g,k,l) = -var(g,k,l,JX)  / GSQRT (g,k,l)
          GZY (g,k,l) = -var(g,k,l,JY)  / GSQRT (g,k,l)
          GZZ (g,k,l) = -var(g,k,l,JZ)  / GSQRT (g,k,l)
          GZXH(g,k,l) = -var(g,k,l,JXH) / GSQRTH(g,k,l)
          GZYH(g,k,l) = -var(g,k,l,JYH) / GSQRTH(g,k,l)
          GZZH(g,k,l) = -var(g,k,l,JZH) / GSQRTH(g,k,l)
       enddo
       enddo

       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          VMTR_C2WfactGz(g,k,I_a_GZXH,l) = GRD_afact(k) * VMTR_RGSGAM2(g,k  ,l) * VMTR_GSGAM2H(g,k,l) * GZXH(g,k,l)
          VMTR_C2WfactGz(g,k,I_b_GZXH,l) = GRD_bfact(k) * VMTR_RGSGAM2(g,k-1,l) * VMTR_GSGAM2H(g,k,l) * GZXH(g,k,l)
          VMTR_C2WfactGz(g,k,I_a_GZYH,l) = GRD_afact(k) * VMTR_RGSGAM2(g,k  ,l) * VMTR_GSGAM2H(g,k,l) * GZYH(g,k,l)
          VMTR_C2WfactGz(g,k,I_b_GZYH,l) = GRD_bfact(k) * VMTR_RGSGAM2(g,k-1,l) * VMTR_GSGAM2H(g,k,l) * GZYH(g,k,l)
          VMTR_C2WfactGz(g,k,I_a_GZZH,l) = GRD_afact(k) * VMTR_RGSGAM2(g,k  ,l) * VMTR_GSGAM2H(g,k,l) * GZZH(g,k,l)
          VMTR_C2WfactGz(g,k,I_b_GZZH,l) = GRD_bfact(k) * VMTR_RGSGAM2(g,k-1,l) * VMTR_GSGAM2H(g,k,l) * GZZH(g,k,l)
       enddo
       enddo

       do g = 1, ADM_gall
          VMTR_C2WfactGz(g,ADM_kmin-1,I_a_GZXH,l) = 0.0_RP
          VMTR_C2WfactGz(g,ADM_kmin-1,I_b_GZXH,l) = 0.0_RP
          VMTR_C2WfactGz(g,ADM_kmin-1,I_a_GZYH,l) = 0.0_RP
          VMTR_C2WfactGz(g,ADM_kmin-1,I_b_GZYH,l) = 0.0_RP
          VMTR_C2WfactGz(g,ADM_kmin-1,I_a_GZZH,l) = 0.0_RP
          VMTR_C2WfactGz(g,ADM_kmin-1,I_b_GZZH,l) = 0.0_RP
       enddo
    enddo

    !--- calculation of volume, geopotential
    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       VMTR_VOLUME(g,k,l) = GMTR_area(g,l) * VMTR_GSGAM2(g,k,l) * GRD_dgz(k)

       VMTR_PHI   (g,k,l) = GRD_vz(g,k,l,GRD_Z) * GRAV
    enddo
    enddo
    enddo



    if ( ADM_have_pl ) then

       !---   G^1/2 = dz/dgz
       do l = 1, ADM_lall_pl
          !--- calculation of G^1/2 at full level
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             GSQRT_pl(g,k,l) = ( GRD_vz_pl(g,k+1,l,GRD_ZH) - GRD_vz_pl(g,k,l,GRD_ZH) ) / GRD_dgz(k)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             GSQRT_pl(g,ADM_kmin-1,l) = GSQRT_pl(g,ADM_kmin,l)
             GSQRT_pl(g,ADM_kmax+1,l) = GSQRT_pl(g,ADM_kmax,l)
          enddo
          !--- calculation of G^1/2 at half level
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             GSQRTH_pl(g,k,l) = ( GRD_vz_pl(g,k,l,GRD_Z) - GRD_vz_pl(g,k-1,l,GRD_Z) ) / GRD_dgzh(k)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             GSQRTH_pl(g,ADM_kmin-1,l) = GSQRTH_pl(g,ADM_kmin,l)
          enddo
       enddo

       !--- Gamma = (a+z) / a
       if ( VMTR_deep_atmos ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             GAM_pl (g,k,l) = 1.0_RP + GRD_vz_pl(g,k,l,GRD_Z)  / GRD_rscale
             GAMH_pl(g,k,l) = 1.0_RP + GRD_vz_pl(g,k,l,GRD_ZH) / GRD_rscale
          enddo
          enddo
          enddo
       else
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             GAM_pl (g,k,l) = 1.0_RP
             GAMH_pl(g,k,l) = 1.0_RP
          enddo
          enddo
          enddo
       endif

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          VMTR_GAM2H_pl   (g,k,l) = GAMH_pl(g,k,l) * GAMH_pl(g,k,l)
          VMTR_GSGAM2_pl  (g,k,l) = GAM_pl (g,k,l) * GAM_pl (g,k,l) * GSQRT_pl (g,k,l)
          VMTR_GSGAM2H_pl (g,k,l) = GAMH_pl(g,k,l) * GAMH_pl(g,k,l) * GSQRTH_pl(g,k,l)

          VMTR_RGSQRTH_pl (g,k,l) = 1.0_RP / GSQRTH_pl(g,k,l)
          VMTR_RGAM_pl    (g,k,l) = 1.0_RP / GAM_pl (g,k,l)
          VMTR_RGAMH_pl   (g,k,l) = 1.0_RP / GAMH_pl(g,k,l)
          VMTR_RGSGAM2_pl (g,k,l) = 1.0_RP / VMTR_GSGAM2_pl (g,k,l)
          VMTR_RGSGAM2H_pl(g,k,l) = 1.0_RP / VMTR_GSGAM2H_pl(g,k,l)
       enddo
       enddo
       enddo

       ! full level <-> half level interpolation factor
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             VMTR_C2Wfact_pl(g,k,I_a,l) = GRD_afact(k) * VMTR_RGSGAM2_pl(g,k  ,l) * VMTR_GSGAM2H_pl(g,k,l)
             VMTR_C2Wfact_pl(g,k,I_b,l) = GRD_bfact(k) * VMTR_RGSGAM2_pl(g,k-1,l) * VMTR_GSGAM2H_pl(g,k,l)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             VMTR_C2Wfact_pl(g,ADM_kmin-1,I_a,l) = 0.0_RP
             VMTR_C2Wfact_pl(g,ADM_kmin-1,I_b,l) = 0.0_RP
          enddo

          do k = ADM_kmin-1, ADM_kmax
          do g = 1, ADM_gall_pl
             VMTR_W2Cfact_pl(g,k,I_c,l) = GRD_cfact(k) * VMTR_GSGAM2_pl(g,k,l) * VMTR_RGSGAM2H_pl(g,k+1,l)
             VMTR_W2Cfact_pl(g,k,I_d,l) = GRD_dfact(k) * VMTR_GSGAM2_pl(g,k,l) * VMTR_RGSGAM2H_pl(g,k  ,l)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             VMTR_W2Cfact_pl(g,ADM_kmax+1,I_c,l) = 0.0_RP
             VMTR_W2Cfact_pl(g,ADM_kmax+1,I_d,l) = 0.0_RP
          enddo
       enddo

       ! full level <-> half level interpolation factor with Gz

       !--- Gz(X) = - JX / G^1/2
       !--- Gz(Y) = - JY / G^1/2
       !--- Gz(Z) = - JZ / G^1/2
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             GZX_pl (g,k,l) = -var_pl(g,k,l,JX)  / GSQRT_pl (g,k,l)
             GZY_pl (g,k,l) = -var_pl(g,k,l,JY)  / GSQRT_pl (g,k,l)
             GZZ_pl (g,k,l) = -var_pl(g,k,l,JZ)  / GSQRT_pl (g,k,l)
             GZXH_pl(g,k,l) = -var_pl(g,k,l,JXH) / GSQRTH_pl(g,k,l)
             GZYH_pl(g,k,l) = -var_pl(g,k,l,JYH) / GSQRTH_pl(g,k,l)
             GZZH_pl(g,k,l) = -var_pl(g,k,l,JZH) / GSQRTH_pl(g,k,l)
          enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             VMTR_C2WfactGz_pl(g,k,I_a_GZXH,l) = GRD_afact(k) * VMTR_RGSGAM2_pl(g,k  ,l) * VMTR_GSGAM2H_pl(g,k,l) * GZXH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_b_GZXH,l) = GRD_bfact(k) * VMTR_RGSGAM2_pl(g,k-1,l) * VMTR_GSGAM2H_pl(g,k,l) * GZXH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_a_GZYH,l) = GRD_afact(k) * VMTR_RGSGAM2_pl(g,k  ,l) * VMTR_GSGAM2H_pl(g,k,l) * GZYH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_b_GZYH,l) = GRD_bfact(k) * VMTR_RGSGAM2_pl(g,k-1,l) * VMTR_GSGAM2H_pl(g,k,l) * GZYH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_a_GZZH,l) = GRD_afact(k) * VMTR_RGSGAM2_pl(g,k  ,l) * VMTR_GSGAM2H_pl(g,k,l) * GZZH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_b_GZZH,l) = GRD_bfact(k) * VMTR_RGSGAM2_pl(g,k-1,l) * VMTR_GSGAM2H_pl(g,k,l) * GZZH_pl(g,k,l)
          enddo
          enddo

          do g = 1, ADM_gall_pl
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_a_GZXH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_b_GZXH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_a_GZYH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_b_GZYH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_a_GZZH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_b_GZZH,l) = 0.0_RP
          enddo
       enddo

       !--- calculation of volume, geopotential
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          VMTR_VOLUME_pl(g,k,l) = GMTR_area_pl(g,l) * VMTR_GSGAM2_pl(g,k,l) * GRD_dgz(k)

          VMTR_PHI_pl   (g,k,l) = GRD_vz_pl(g,k,l,GRD_Z) * GRAV
       enddo
       enddo
       enddo
    else
       VMTR_GAM2H_pl    (:,:,:)   = 0.0_RP
       VMTR_GSGAM2_pl   (:,:,:)   = 0.0_RP
       VMTR_GSGAM2H_pl  (:,:,:)   = 0.0_RP
       VMTR_RGSQRTH_pl  (:,:,:)   = 0.0_RP
       VMTR_RGAM_pl     (:,:,:)   = 0.0_RP
       VMTR_RGAMH_pl    (:,:,:)   = 0.0_RP
       VMTR_RGSGAM2_pl  (:,:,:)   = 0.0_RP
       VMTR_RGSGAM2H_pl (:,:,:)   = 0.0_RP
       VMTR_W2Cfact_pl  (:,:,:,:) = 0.0_RP
       VMTR_C2Wfact_pl  (:,:,:,:) = 0.0_RP
       VMTR_C2WfactGz_pl(:,:,:,:) = 0.0_RP
       VMTR_VOLUME_pl   (:,:,:)   = 0.0_RP
       VMTR_PHI_pl      (:,:,:)   = 0.0_RP
    endif

    return
  end subroutine VMTR_setup

end module mod_vmtr
