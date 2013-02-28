! To do list
!
! - change Gamma (parameter) --> (variable)
! - change Pr(parameter) --> (variable)
! - Effect of stratos
! - Energy_tot ok?
! - potem_l    ok?
! 
!-------------------------------------------------------------------------------
!
!+  Smagorinsky turbulent diffusion module
!
!-------------------------------------------------------------------------------
module mod_tb_smg
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines for Smagorinsky-type turbulent diffusion
  !       
  ! 
  !++ Current Corresponding Authors : S.Iga and A.Noda
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      11-xx-xx   generated
  !                10-12-03   Iga modified
  !                10-12-27   Iga modified
  !                11-11-28   Y.Yamada: Merge Terai-san code into
  !                                                         the original code.
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only :  &
       ADM_log_fid, &
       ADM_prc_me, &
       ADM_prc_pl, &
       ADM_gall, &
       ADM_kall, &
       ADM_lall, &
       ADM_gall_pl, &
       ADM_lall_pl, &
       ADM_knone,        &
       ADM_kmin,         &
       ADM_kmax, &
       ADM_NSYS
    use mod_vmtr,only:&
!        VMTR_GSCAM,&
         VMTR_GSGAM2,&
         VMTR_RGAM2,&
         VMTR_RGAM2H,&
         VMTR_RGSGAM2,&
         VMTR_RGSGAM2H,&
         VMTR_GAM2,&
         VMTR_GAM2H,&
         VMTR_GSGAM2H,   &
         VMTR_GZX,&
         VMTR_GZXH,&
         VMTR_GZY,&
         VMTR_GZZ,&
         VMTR_GZYH,&
         VMTR_GZZH,&
         VMTR_RGSH,&
         VMTR_RGAM,&
         VMTR_GAM2_PL,&
         VMTR_GAM2H_PL,&
!        VMTR_GSCAM_PL,&
         VMTR_GSGAM2H_pl,&
         VMTR_GSGAM2_PL,&
         VMTR_RGAM2_PL,&
         VMTR_RGAM2H_PL,&
         VMTR_RGSGAM2_PL,&
         VMTR_RGSGAM2H_PL,&
         VMTR_GZX_PL,&
         VMTR_GZY_PL,&
         VMTR_GZXH_PL,&
         VMTR_GZYH_PL,&
         VMTR_GZZ_PL,&
         VMTR_GZZH_PL,&
         VMTR_RGSH_PL,&
         VMTR_RGAM_PL
    use mod_history
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: tb_smg_setup
  public :: tb_smg_driver
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !------ Dependency of horizontal grid size 
  logical, private, save :: dep_hgrid = .false.
  !------ coefficient for diffusion (horizontal)
  integer, private, save :: lap_order_hdiff = 2
  !
  real(8), private, save :: K_coef_minlim = 0!3.0D+9
  real(8), private, save :: K_coef_maxlim = 1.0D+99
  real(8), private, save :: LENGTH_maxlim = 1.0D+99
  real(8), private, save :: SMG_CS = 0.2d0
  real(8), private, save :: SMALL=1.0d-10
  real(8), private, save :: GAMMA=1.0d0  ! (kind of aspect ratio?)
  real(8), private, save :: Pr=1.0d0
  logical, private, save :: stratos_effect = .false.
  real(8), private, save :: beta_theta!=1.0d0
  real(8), private, save :: beta_q!=1.0d0


!  logical, private, save :: DEEP_EFFECT = .true. ! --> why ?  meaningless here

  real(8), private, save :: horiz_dx2
  !
  logical, private, save :: first = .true.
  !
  ! for smg_oprt

  real(8), allocatable, private, save:: smg_oprt_cxh(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cxh_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cyh(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cyh_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_czh(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_czh_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cx(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cx_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cy(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cy_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cz(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_cz_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_rgamH(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_rgamH_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_rgam(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_rgam_pl(:,:,:)
!  real(8), allocatable, private, save:: smg_oprt_GzGz(:,:,:)
!  real(8), allocatable, private, save:: smg_oprt_GzGz_pl(:,:,:)
!  real(8), allocatable, private, save:: smg_oprt_GzGzh(:,:,:)
!  real(8), allocatable, private, save:: smg_oprt_GzGzh_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_gsqrt(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_gsqrt_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_gsqrtH(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_gsqrtH_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_GAM(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_GAM_pl(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_GAMH(:,:,:)
  real(8), allocatable, private, save:: smg_oprt_GAMH_pl(:,:,:)

  real(8), allocatable, private, save:: var(:,:,:,:,:)
  real(8), allocatable, private, save:: var_pl(:,:,:,:,:)
  real(8), allocatable, private, save:: varh(:,:,:,:,:)
  real(8), allocatable, private, save:: varh_pl(:,:,:,:,:)

   integer, private, parameter ::  IVX=1
   integer, private, parameter ::  IVY=2
   integer, private, parameter ::  IVZ=3
   integer, private, parameter ::  ipotem=4
   integer, private, parameter ::  IMAX=ipotem


  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine tb_smg_setup
    !
    use mod_adm, only :  &
       ADM_ctl_fid
    !
    implicit none
    !
    integer :: ierr
                                  
    namelist / SMGPARAM / &
         dep_hgrid,             & !--- depnedency of horzontal grid size
         K_coef_minlim,        & !--- K_coef (minimum limit)
         K_coef_maxlim,        & !--- K_coef (maximum limit)
         LENGTH_maxlim,       & !--- 
         SMG_CS,                & !--- Smagorinsky constant
         stratos_effect                 !--- stratos effect (default: false)
!         DEEP_EFFECT            
    !
!    integer :: k,l,n
    !
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=SMGPARAM,iostat = ierr)
    if(ierr<0) then
       write(ADM_LOG_FID,*) &
            'Msg : Sub[smg_setup]/Mod[smg]'
       write(ADM_LOG_FID,*) &
            ' *** Not found namelist.'
       write(ADM_LOG_FID,*) &
            ' *** Use default values.'
    else if(ierr>0) then
       write(*,*) &
            'Msg : Sub[smg_setup]/Mod[smg]'
       write(*,*) &
            ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
    end if
    !
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=SMGPARAM,iostat = ierr)
    write(ADM_LOG_FID,nml=SMGPARAM)

!    if ((trim(divdamp_type)/="DIRECT").or.(trim(hdiff_type)/="DIRECT")) then
!       area(:,ADM_KNONE,:) = GMTR_area(:,:)
!       area_pl(:,ADM_KNONE,:) = GMTR_area_pl(:,:)
!       horiz_dx2=GTL_global_mean_without_aw(area,area_pl)
!    end if
    !

    call tb_smg_oprt_init()

  end subroutine tb_smg_setup
  !-----------------------------------------------------------------------------
  subroutine tb_smg_driver(&
       nl,&
       rho,   rho_pl,                  &  !--- IN : density
       rhog,  rhog_pl,                 &  !--- IN : density
       rhogq , rhogq_pl,                       &  !--- IN : q
       vx, vx_pl,                      &  !--- IN : Vx
       vy, vy_pl,                      &  !--- IN : Vy
       vz, vz_pl,                      &  !--- IN : Vz
       w , w_pl,                       &  !--- IN : w
!       temd, temd_pl,                  &  !--- IN : temperature
       tem, tem_pl,                  &  !--- IN : temperature
       q , q_pl,                       &  !--- IN : q
       potem, potem_pl,                &  !--- IN : potential temperature (please prepare)
       frhog, frhog_pl,                &  !--- INOUT : tend. of rhog
       frhogvx, frhogvx_pl,            &  !--- INOUT : tend. of rhogvx
       frhogvy, frhogvy_pl,            &  !--- INOUT : tend. of rhogvy
       frhogvz, frhogvz_pl,            &  !--- INOUT : tend. of rhogvz
       frhogw,  frhogw_pl,             &  !--- INOUT : tend. of rhogw
       frhoge,  frhoge_pl,             &  !--- INOUT : tend. of rhoge
       frhogetot,  frhogetot_pl,       &  !--- INOUT : tend. of rhoge
       frhogq,  frhogq_pl              &  !--- INOUT : tend. of rhogq
       )
    !------ 
    !------ Numerical diffusion
    !------     1. Calculation region of frhogvx,...., frhoge
    !------                    : (:,:,:)
    !------ 
    !
    use mod_adm, only : &
         ADM_comm_run_world
    use mod_gmtr, only :  &
         GMTR_area,&
         GMTR_area_pl
    use mod_grd, only :  &
         GRD_rscale,&
         GRD_xdir,&
         GRD_ydir,&
         GRD_zdir,&
         GRD_x,&
         GRD_x_pl,&
!         GRD_dgz,        &!101201
         GRD_vz,       &!101201
         GRD_Z,        &!101201
         GRD_ZH,        &
         GRD_vz_pl,    &
         GRD_afac,       &
         GRD_bfac,       &
         GRD_gz,         &
         GRD_gzh
    use mod_cnst, only : &
         CNST_CV,        &
         CNST_CP,        &
         CNST_PI
    use mod_runconf, only : &
         I_QV
    use mod_oprt, only :        &
         OPRT_horizontalize_vec,&
         OPRT_laplacian,        &
         OPRT_diffusion
    use mod_comm, only :        &
         COMM_data_transfer
    use mod_runconf, only :         &
         TRC_VMAX
    use mod_bsstate, only :         &
         rho_bs,rho_bs_pl
    use mod_time, only : &
         TIME_DTL
    !
    implicit none
    !
    integer, intent(in) :: nl
    real(8), intent(in) :: rho(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rho_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhog(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhog_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogq(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8), intent(in) :: rhogq_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,TRC_VMAX)
    !
    real(8), intent(in) :: vx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: w(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: w_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
!    real(8), intent(inout) :: w(ADM_gall,ADM_kall,ADM_lall)! inout for debug
!    real(8), intent(inout) :: w_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)! inout for debug
!    real(8), intent(in) :: temd(ADM_gall,ADM_kall,ADM_lall)
!    real(8), intent(in) :: temd_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: tem(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: tem_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: q(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8), intent(in) :: q_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,TRC_VMAX)
    !
    real(8), intent(inout) :: frhog(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhog_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    real(8), intent(inout) :: frhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogw_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhoge(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhoge_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogetot(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogetot_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogq(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8), intent(inout) :: frhogq_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,TRC_VMAX)

    !full level growth late
    real(8) :: grhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhogvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: grhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhogvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: grhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhogvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: grhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhogw_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: grhoge(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhoge_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
!    real(8) :: grhogetot(ADM_gall,ADM_kall,ADM_lall)
!    real(8) :: grhogetot_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: grhogq(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: grhogq_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,TRC_VMAX)
    !half level growth late
    real(8) :: grhogvxh(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhogvxh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: grhogvyh(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhogvyh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: grhogvzh(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhogvzh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8), intent(in) :: potem(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: potem_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: pi(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: pi_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    integer :: p, id
    integer:: i,j,k,i2,j2,l,n,  idir,ivar
    integer :: nq
    !
    real(8) :: del_xyz2, SMG_CS2, LENGTH_maxlim2
    real(8)::sij(ADM_gall,ADM_kall,ADM_lall,3,3)           ! full level
    real(8)::sij_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,3,3)  ! full level
    real(8)::sijsij(ADM_gall,ADM_kall,ADM_lall)           ! full level
    real(8)::sijsij_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! full level
    real(8)::sijsijh(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8)::sijsijh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    real(8)::sijh(ADM_gall,ADM_kall,ADM_lall,3,3)
    real(8)::sijh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,3,3)  ! partial ui / partial xj
    !
    !------ diffusion coefficient
    real(8)::K_coef(ADM_gall,ADM_kall,ADM_lall)            ! full level
    real(8)::K_coef_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)            ! full level
    real(8)::K_coefh(ADM_gall,ADM_kall,ADM_lall)            ! half level
    real(8)::K_coefh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)            ! half level
    !
    real(8)::abs_vxh(ADM_gall,ADM_kall,ADM_lall)            ! half level
    real(8)::abs_vxh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)   ! half level
    real(8)::abs_vyh(ADM_gall,ADM_kall,ADM_lall)            ! half level
    real(8)::abs_vyh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)   ! half level
    real(8)::abs_vzh(ADM_gall,ADM_kall,ADM_lall)            ! half level
    real(8)::abs_vzh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)   ! half level
    !
    real(8)::abs_vx(ADM_gall,ADM_kall,ADM_lall)            ! full level
    real(8)::abs_vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)   ! full level
    real(8)::abs_vy(ADM_gall,ADM_kall,ADM_lall)            ! full level
    real(8)::abs_vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)   ! full level
    real(8)::abs_vz(ADM_gall,ADM_kall,ADM_lall)            ! full level
    real(8)::abs_vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)   ! full level
   
    ! tensor
    real(8)::uijh(ADM_gall,ADM_kall,ADM_lall,3,3)
    real(8)::uijh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,3,3)  ! partial ui / partial xj
    real(8)::uij(ADM_gall,ADM_kall,ADM_lall,3,3)
    real(8)::uij_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,3,3)  ! partial ui / partial xj
   

    real(8) :: stratos(ADM_GALL,ADM_kall,ADM_LALL)
    real(8) :: stratos_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: stratosh(ADM_GALL,ADM_kall,ADM_LALL)
    real(8) :: stratosh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    !
    real(8) :: fq(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: fq_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: wrk(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: wrk_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,TRC_VMAX)
    real(8) :: wrkh(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: wrkh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,TRC_VMAX)
    !
    real(8) :: rhoh(ADM_gall,ADM_kall,ADM_lall)  ! rho at half level
    real(8) :: rhoh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! rho at half level
    real(8) :: pih(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: pih_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    real(8) :: rho_K_Cp_PI(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rho_K_Cp_PIh(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rho_K_Cp_PI_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rho_K_Cp_PIh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    real(8) :: dummy(ADM_GALL,ADM_kall,ADM_LALL)
    real(8) :: dummy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
!   logical,save:: dbgfirst=.true.
    logical,save:: dbgfirst=.false.

    real(8) :: wrkwrk(ADM_GALL,ADM_kall,ADM_LALL)
    real(8) :: wrkwrk2(ADM_GALL,ADM_kall,ADM_LALL)

    integer :: ierr

    !================================= (1) calculate velocity in cartesian coordinate ==========================

!w=0
!w_pl=0

    !
    ! calculate velocity in absolute (cartesian) coordinate
     !
    abs_vx=0
    abs_vy=0
    abs_vz=0
    abs_vxh=0
    abs_vyh=0
    abs_vzh=0
    abs_vx_pl=0
    abs_vy_pl=0
    abs_vz_pl=0
    abs_vxh_pl=0
    abs_vyh_pl=0
    abs_vzh_pl=0

    ! half level
    !   k=ADM_kmin is ground, and ground velocity is assumed to be zero
    abs_vxh(:,ADM_kmin,:) = 0!vx(:,ADM_kmin,:)/2 + w(:,ADM_kmin,:) * GRD_x(:,ADM_knone,:,GRD_XDIR)/GRD_rscale
    abs_vyh(:,ADM_kmin,:) = 0!vy(:,ADM_kmin,:)/2 + w(:,ADM_kmin,:) * GRD_x(:,ADM_knone,:,GRD_YDIR)/GRD_rscale
    abs_vzh(:,ADM_kmin,:) = 0!vz(:,ADM_kmin,:)/2 + w(:,ADM_kmin,:) * GRD_x(:,ADM_knone,:,GRD_ZDIR)/GRD_rscale
    do k=ADM_kmin+1,ADM_kmax
       abs_vxh(:,k,:) = (vx(:,k,:)+vx(:,k-1,:))/2 + w(:,k,:) * GRD_x(:,ADM_knone,:,GRD_XDIR)/GRD_rscale
       abs_vyh(:,k,:) = (vy(:,k,:)+vy(:,k-1,:))/2 + w(:,k,:) * GRD_x(:,ADM_knone,:,GRD_YDIR)/GRD_rscale
       abs_vzh(:,k,:) = (vz(:,k,:)+vz(:,k-1,:))/2 + w(:,k,:) * GRD_x(:,ADM_knone,:,GRD_ZDIR)/GRD_rscale
    enddo

    if(ADM_prc_me==ADM_prc_pl) then
      abs_vxh_pl(:,ADM_kmin,:) = 0!vx_pl(:,ADM_kmin,:)/2 + w_pl(:,ADM_kmin,:) * GRD_x_pl(:,ADM_knone,:,GRD_XDIR)/GRD_rscale
      abs_vyh_pl(:,ADM_kmin,:) = 0!vy_pl(:,ADM_kmin,:)/2 + w_pl(:,ADM_kmin,:) * GRD_x_pl(:,ADM_knone,:,GRD_YDIR)/GRD_rscale
      abs_vzh_pl(:,ADM_kmin,:) = 0!vz_pl(:,ADM_kmin,:)/2 + w_pl(:,ADM_kmin,:) * GRD_x_pl(:,ADM_knone,:,GRD_ZDIR)/GRD_rscale
      do k=ADM_kmin+1,ADM_kmax
       abs_vxh_pl(:,k,:) = (vx_pl(:,k,:)+vx_pl(:,k-1,:))/2 + w_pl(:,k,:) * GRD_x_pl(:,ADM_knone,:,GRD_XDIR)/GRD_rscale
       abs_vyh_pl(:,k,:) = (vy_pl(:,k,:)+vy_pl(:,k-1,:))/2 + w_pl(:,k,:) * GRD_x_pl(:,ADM_knone,:,GRD_YDIR)/GRD_rscale
       abs_vzh_pl(:,k,:) = (vz_pl(:,k,:)+vz_pl(:,k-1,:))/2 + w_pl(:,k,:) * GRD_x_pl(:,ADM_knone,:,GRD_ZDIR)/GRD_rscale
      enddo
    endif
    !
    ! full level
    !   w(k=ADM_kmax+1) is top  and I don't know about its treatment
    do k=ADM_kmin,ADM_kmax-1
       abs_vx(:,k,:) = vx(:,k,:) + (w(:,k,:)+w(:,k+1,:))/2 * GRD_x(:,ADM_knone,:,GRD_XDIR)/GRD_rscale
       abs_vy(:,k,:) = vy(:,k,:) + (w(:,k,:)+w(:,k+1,:))/2 * GRD_x(:,ADM_knone,:,GRD_YDIR)/GRD_rscale
       abs_vz(:,k,:) = vz(:,k,:) + (w(:,k,:)+w(:,k+1,:))/2 * GRD_x(:,ADM_knone,:,GRD_ZDIR)/GRD_rscale
    enddo
    k=ADM_kmax
    abs_vx(:,k,:) = vx(:,k,:) + (w(:,k,:)+w(:,k,:))/2 * GRD_x(:,ADM_knone,:,GRD_XDIR)/GRD_rscale
    abs_vy(:,k,:) = vy(:,k,:) + (w(:,k,:)+w(:,k,:))/2 * GRD_x(:,ADM_knone,:,GRD_YDIR)/GRD_rscale
    abs_vz(:,k,:) = vz(:,k,:) + (w(:,k,:)+w(:,k,:))/2 * GRD_x(:,ADM_knone,:,GRD_ZDIR)/GRD_rscale

    if(ADM_prc_me==ADM_prc_pl) then 
      do k=ADM_kmin,ADM_kmax-1
       abs_vx_pl(:,k,:) = vx_pl(:,k,:) + (w_pl(:,k,:)+w_pl(:,k+1,:))/2 * GRD_x_pl(:,ADM_knone,:,GRD_XDIR)/GRD_rscale
       abs_vy_pl(:,k,:) = vy_pl(:,k,:) + (w_pl(:,k,:)+w_pl(:,k+1,:))/2 * GRD_x_pl(:,ADM_knone,:,GRD_YDIR)/GRD_rscale
       abs_vz_pl(:,k,:) = vz_pl(:,k,:) + (w_pl(:,k,:)+w_pl(:,k+1,:))/2 * GRD_x_pl(:,ADM_knone,:,GRD_ZDIR)/GRD_rscale
      enddo
      k=ADM_kmax
      abs_vx_pl(:,k,:) = vx_pl(:,k,:) + (w_pl(:,k,:)+w_pl(:,k,:))/2 * GRD_x_pl(:,ADM_knone,:,GRD_XDIR)/GRD_rscale
      abs_vy_pl(:,k,:) = vy_pl(:,k,:) + (w_pl(:,k,:)+w_pl(:,k,:))/2 * GRD_x_pl(:,ADM_knone,:,GRD_YDIR)/GRD_rscale
      abs_vz_pl(:,k,:) = vz_pl(:,k,:) + (w_pl(:,k,:)+w_pl(:,k,:))/2 * GRD_x_pl(:,ADM_knone,:,GRD_ZDIR)/GRD_rscale
    endif

!    call dbgmx('abs_vx',abs_vx)
!    call dbgmx('abs_vxh',abs_vxh)

    do l=1,ADM_lall
       if (dbgfirst)then

       call history_in('vx',vx(:,:,l))
       call history_in('vy',vy(:,:,l))
       call history_in('vz',vz(:,:,l))
       call history_in('w',w(:,:,l))
       call history_in('tem',tem(:,:,l))
       call history_in('rho',rho(:,:,l))
       call history_in('abs_vx',abs_vx(:,:,l))
       call history_in('abs_vxh',abs_vxh(:,:,l))
       call history_in('abs_vy',abs_vy(:,:,l))
       call history_in('abs_vyh',abs_vyh(:,:,l))
       call history_in('abs_vz',abs_vz(:,:,l))
       call history_in('abs_vzh',abs_vzh(:,:,l))

       call history_in('smg_oprt_cx',smg_oprt_cx(:,:,l))
       call history_in('smg_oprt_cy',smg_oprt_cy(:,:,l))
       call history_in('smg_oprt_cz',smg_oprt_cz(:,:,l))
       call history_in('smg_oprt_cxh',smg_oprt_cxh(:,:,l))
       call history_in('smg_oprt_cyh',smg_oprt_cyh(:,:,l))
       call history_in('smg_oprt_czh',smg_oprt_czh(:,:,l))

       call history_in('smg_oprt_gsqrt',smg_oprt_gsqrt(:,:,l))!1
       call history_in('smg_oprt_gsqrth',smg_oprt_gsqrth(:,:,l))!1
       call history_in('smg_oprt_gam',smg_oprt_gam(:,:,l))!1
       call history_in('smg_oprt_gamh',smg_oprt_gamh(:,:,l))!1
       call history_in('smg_oprt_rgam',smg_oprt_rgam(:,:,l))!1
       call history_in('smg_oprt_rgamh',smg_oprt_rgamh(:,:,l))!1

       endif
    enddo

    !=============================== (2) calculate Sij and K_coefh at both full and half level ====================
!    call dbgmx('aa oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('aa oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))

    sijh=0.d0
    sijh_pl=0.d0
    uijh=0.d0
    uijh_pl=0.d0

    sij=0.d0
    sij_pl=0.d0
    uij=0.d0
    uij_pl=0.d0

    ! x
!    write(*,*) 'x'
    call Gradient3dfh(        &
         uij(:,:,:,1,1),uij_pl(:,:,:,1,1),&
         uij(:,:,:,1,2),uij_pl(:,:,:,1,2),&
         uij(:,:,:,1,3),uij_pl(:,:,:,1,3),&
         uijh(:,:,:,1,1),uijh_pl(:,:,:,1,1),&
         uijh(:,:,:,1,2),uijh_pl(:,:,:,1,2),&
         uijh(:,:,:,1,3),uijh_pl(:,:,:,1,3),&
         abs_vx(:,:,:), abs_vx_pl(:,:,:),          &
         abs_vxh(:,:,:), abs_vxh_pl(:,:,:),          &
         .true.&
          )

    ! y
!    write(*,*) 'y'
    call Gradient3dfh(        &
         uij(:,:,:,2,1),uij_pl(:,:,:,2,1),&
         uij(:,:,:,2,2),uij_pl(:,:,:,2,2),&
         uij(:,:,:,2,3),uij_pl(:,:,:,2,3),&
         uijh(:,:,:,2,1),uijh_pl(:,:,:,2,1),&
         uijh(:,:,:,2,2),uijh_pl(:,:,:,2,2),&
         uijh(:,:,:,2,3),uijh_pl(:,:,:,2,3),&
         abs_vy(:,:,:), abs_vy_pl(:,:,:),          &
         abs_vyh(:,:,:), abs_vyh_pl(:,:,:),          &
         .true.&
          )

    ! z
!    write(*,*) 'z'
    call Gradient3dfh(        &
         uij(:,:,:,3,1),uij_pl(:,:,:,3,1),&
         uij(:,:,:,3,2),uij_pl(:,:,:,3,2),&
         uij(:,:,:,3,3),uij_pl(:,:,:,3,3),&
         uijh(:,:,:,3,1),uijh_pl(:,:,:,3,1),&
         uijh(:,:,:,3,2),uijh_pl(:,:,:,3,2),&
         uijh(:,:,:,3,3),uijh_pl(:,:,:,3,3),&
         abs_vz(:,:,:), abs_vz_pl(:,:,:),          &
         abs_vzh(:,:,:), abs_vzh_pl(:,:,:),          &
         .true.&
          )

!    call dbgmx('bb oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('bb oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))
    do l=1,ADM_lall
       if (dbgfirst)then

          call history_in('uij11',uij(:,:,l,1,1))
          call history_in('uijh11',uijh(:,:,l,1,1))
          call history_in('uij21',uij(:,:,l,2,1))
          call history_in('uijh21',uijh(:,:,l,2,1))
          call history_in('uij31',uij(:,:,l,3,1))
          call history_in('uijh31',uijh(:,:,l,3,1))
          call history_in('uij13',uij(:,:,l,1,3))
          call history_in('uijh13',uijh(:,:,l,1,3))

       endif
    enddo

    do i=1,3
       do j=1,3
          sij(:,ADM_kmin:ADM_kmax,:,i,j)=(uij(:,ADM_kmin:ADM_kmax,:,i,j)+uij(:,ADM_kmin:ADM_kmax,:,j,i))/2.d0
          sij_pl(:,ADM_kmin:ADM_kmax,:,i,j)=(uij_pl(:,ADM_kmin:ADM_kmax,:,i,j)+uij_pl(:,ADM_kmin:ADM_kmax,:,j,i))/2.d0
          sijh(:,ADM_kmin:ADM_kmax,:,i,j)=(uijh(:,ADM_kmin:ADM_kmax,:,i,j)+uijh(:,ADM_kmin:ADM_kmax,:,j,i))/2.d0
          sijh_pl(:,ADM_kmin:ADM_kmax,:,i,j)=(uijh_pl(:,ADM_kmin:ADM_kmax,:,i,j)+uijh_pl(:,ADM_kmin:ADM_kmax,:,j,i))/2.d0
       enddo
    enddo
    sijsij(:,:,:) = 0.d0
    sijsij_pl(:,:,:) = 0.d0
    sijsijh(:,:,:) = 0.d0
    sijsijh_pl(:,:,:) = 0.d0
    do i=1,3
       do j=1,3
          sijsij(:,ADM_kmin:ADM_kmax,:)=sijsij(:,ADM_kmin:ADM_kmax,:)+&
               sij(:,ADM_kmin:ADM_kmax,:,i,j)*sij(:,ADM_kmin:ADM_kmax,:,i,j)
          sijsij_pl(:,ADM_kmin:ADM_kmax,:)=sijsij_pl(:,ADM_kmin:ADM_kmax,:)+&
               sij_pl(:,ADM_kmin:ADM_kmax,:,i,j)*sij_pl(:,ADM_kmin:ADM_kmax,:,i,j)
          sijsijh(:,ADM_kmin:ADM_kmax,:)=sijsijh(:,ADM_kmin:ADM_kmax,:)+&
               sijh(:,ADM_kmin:ADM_kmax,:,i,j)*sijh(:,ADM_kmin:ADM_kmax,:,i,j)
          sijsijh_pl(:,ADM_kmin:ADM_kmax,:)=sijsijh_pl(:,ADM_kmin:ADM_kmax,:)+&
               sijh_pl(:,ADM_kmin:ADM_kmax,:,i,j)*sijh_pl(:,ADM_kmin:ADM_kmax,:,i,j)
       enddo
    enddo

    SMG_CS2 = SMG_CS**2
    LENGTH_maxlim2=LENGTH_maxlim**2

!!$    ! stratos effect
!!$    if (stratos_effect) then
!!$       !
!!$       ql=0
!!$       if (RAIN_TYPE.eq.'WARM')  ql(:,:,:)=q(:,:,:,IQC)+q(:,:,:,IQR)
!!$       if (MP_TYPE.eq.'NSW5')  ql(:,:,:)=q(:,:,:,IQC)+q(:,:,:,IQR)+q(:,:,:,IQI)+q(:,:,:,IQS) ! ( ice effect is omitted)
!!$       if (MP_TYPE.eq.'NSW6')  ql(:,:,:)=q(:,:,:,IQC)+q(:,:,:,IQR)+q(:,:,:,IQI)+q(:,:,:,IQS)+q(:,:,:,IQG)! ( ice effect is omitted)
!!$       qtot(:,:,:)=ql(:,:,:)+qv(:,:,:)
!!$       !
!!$       do k = ADM_kmin+1, ADM_kmax
!!$          potemh(:,k,:)=(potem(:,k,:)+potem(:,k-1,:))/2
!!$       enddo
!!$       k=ADM_kmin
!!$       potemh(:,k,:)=potem(:,k,:)
!!$       !
!!$       beta_theta = 
!!$
!!$
!!$
!!$       !half
!!$       do l = 1, ADM_lall_pl
!!$          do k = ADM_kmin+1, ADM_kmax
!!$             do n = 1, ADM_gall_pl
!!$                n2= CNST_EGRAV/(potemh(n,k,l) * (potem(n,k,l)-potem(n,k-1,l))/ ( GRD_vz(n,k,l,GRD_Z)-GRD_vz(n,k-1,l,GRD_Z) ))
!!$                ri= n2/max(  (vx(n,k,l)-vx(n,k-1,l))**2+(vy(n,k,l)-vy(n,k-1,l))**2+(vz(n,k,l)-vz(n,k-1,l))**2, small) *&
!!$                     ( GRD_vz(n,k,l,GRD_Z)-GRD_vz(n,k-1,l,GRD_Z) )**2
!!$                qtot = 
!!$
!!$                stratos_pl(n,k,l)=
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$    else
    stratos=0
    stratosh=0
    stratos_pl=0
    stratosh_pl=0
!!$   endif

!    call dbgmx('cc oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('cc oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))


!    call dbgmx('sijsij',sijsij)

    wrkwrk=0.0
    wrkwrk2=0.0

    ! half --------------
    K_coefh= 0.0d0
    K_coefh_pl= 0.0d0
    K_coefh(:,ADM_kmin-1,:) = 0.0d0 ! im not sure
    K_coefh(:,ADM_kmax+1,:) = 0.0d0
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
          do n = 1, ADM_gall
             del_xyz2=GMTR_area(n,l)+( GRD_vz(n,k,l,GRD_Z)-GRD_vz(n,k-1,l,GRD_Z) )**2 ! (length scale)**2 @ half level, (actual height is used)
             del_xyz2=min(del_xyz2,LENGTH_maxlim2)
!            write(*,*) l,k,n,del_xyz2,GMTR_area(n,l),( GRD_vz(n,k,l,GRD_Z)-GRD_vz(n,k-1,l,GRD_Z) )**2
             K_coefh(n,k,l)   = min(max(SMG_CS2*del_xyz2 * sqrt(max(2.0d0*sijsijh(n,k,l)+stratosh(n,k,l),SMALL)),&
                  K_coef_minlim),K_coef_maxlim)
             wrkwrk(n,k,l)=del_xyz2
             wrkwrk2(n,l,l)=n
          enddo
       enddo
    enddo

    !hogehoge
    write(adm_log_fid,*)'chksmg1',nl,grd_rscale

    if(nl==1)then

       call dbgmx('K_coefh',K_coefh(:,:,:))

       do l=1,ADM_lall
          call history_in('K_coefh',K_coefh(:,:,l)) ! sonouchi kesu
          call history_in('K_coefh_smg',K_coefh(:,:,l))
          call history_in('length',wrkwrk(:,:,l))
          call history_in('vx',vx(:,:,l))
          call history_in('vy',vy(:,:,l))
          call history_in('vz',vz(:,:,l))
          call history_in('abs_vx',abs_vx(:,:,l))
          call history_in('uij11',uij(:,:,l,1,1))
          call history_in('region',wrkwrk2(:,:,l))
       enddo

    endif

    call dbgmx('length',wrkwrk(:,:,:))
    call dbgmx('vx',vx(:,:,:))
    call dbgmx('w',w(:,:,:))
    call dbgmx('grdx',GRD_x(:,adm_knone,:,grd_xdir:grd_xdir))
    call dbgmx('absvx',abs_vx(:,:,:))
    call dbgmx('uij11',uij(:,:,:,1,1))
    call dbgmx('uij12',uij(:,:,:,1,2))
    call dbgmx('uij13',uij(:,:,:,1,3))
    call dbgmx('uij21',uij(:,:,:,2,1))
    call dbgmx('uij22',uij(:,:,:,2,2))
    call dbgmx('uij23',uij(:,:,:,2,3))
    call dbgmx('uij31',uij(:,:,:,3,1))
    call dbgmx('uij32',uij(:,:,:,3,2))
    call dbgmx('uij33',uij(:,:,:,3,3))
    call dbgmx('sijsijh',sqrt(max(2.0d0*sijsijh(:,:,:)+stratosh(:,:,:),SMALL)))
    call dbgmx('sijsijha',sqrt(max(2.0d0*sijsijh(:,:,:),SMALL)))
    call dbgmx('sijsijhb',sqrt(max(stratosh(:,:,:),SMALL)))

!write(adm_log_fid,*)'chksmg1', &
!   maxval(sqrt(max(2.0d0*sijsijh(:,:,:)+stratosh(:,:,:),SMALL))), &
!   minval(sqrt(max(2.0d0*sijsijh(:,:,:)+stratosh(:,:,:),SMALL)))
!write(adm_log_fid,*)'chksmg2', maxval(wrkwrk(:,:,:)),minval(wrkwrk(:,:,:))
!write(adm_log_fid,*)'chksmg3', maxval(K_coefh(:,:,:)),minval(K_coefh(:,:,:))

    if(ADM_prc_me==ADM_prc_pl) then
       K_coefh_pl(:,ADM_kmin-1,:) = 0.0d0 ! im not sure
       K_coefh_pl(:,ADM_kmax+1,:) = 0.0d0
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             do n = 1, ADM_gall_pl
                del_xyz2=GMTR_area_pl(n,l)+( GRD_vz_pl(n,k,l,GRD_Z)-GRD_vz_pl(n,k-1,l,GRD_Z) )**2 ! (length scale)**2 @ half level, (actual height is used)
                del_xyz2=min(del_xyz2,LENGTH_maxlim2)
                K_coefh_pl(n,k,l)   = min(max(SMG_CS2*del_xyz2 * sqrt(max(2.0d0*sijsijh_pl(n,k,l)+stratosh_pl(n,k,l),SMALL)),&
                     K_coef_minlim),K_coef_maxlim)
             enddo
          enddo
       enddo
    endif


    ! full --------------
    K_coef= 0.0d0
    K_coef_pl= 0.0d0
    K_coef(:,ADM_kmin-1,:) = 0.0d0
    K_coef(:,ADM_kmax+1,:) = 0.0d0
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
          do n = 1, ADM_gall
             del_xyz2=GMTR_area(n,l)+( GRD_vz(n,k+1,l,GRD_ZH)-GRD_vz(n,k,l,GRD_ZH) )**2 ! (length scale)**2 @ full level, (actual height is used)
             del_xyz2=min(del_xyz2,LENGTH_maxlim2)
             K_coef(n,k,l)   = min(max(SMG_CS2*del_xyz2 * sqrt(max(2.0d0*sijsij(n,k,l)+stratos(n,k,l),SMALL)),&
                  K_coef_minlim),K_coef_maxlim)
          enddo
       enddo
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
       K_coef_pl(:,ADM_kmin-1,:) = 0.0d0
       K_coef_pl(:,ADM_kmax+1,:) = 0.0d0
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             do n = 1, ADM_gall_pl
                del_xyz2=GMTR_area_pl(n,l)+( GRD_vz_pl(n,k+1,l,GRD_ZH)-GRD_vz_pl(n,k,l,GRD_ZH) )**2 ! (length scale)**2 @ full level, (actual height is used)
                del_xyz2=min(del_xyz2,LENGTH_maxlim2)
                K_coef_pl(n,k,l)  = min(max(SMG_CS2*del_xyz2 * sqrt(max(2.0d0*sijsij_pl(n,k,l)+stratos_pl(n,k,l),SMALL))&
                  ,K_coef_minlim),K_coef_maxlim)
             enddo
          enddo
       enddo
    endif

!    call dbgmx('dd oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('dd oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))

    do l=1,ADM_lall
       if (dbgfirst)then

          call history_in('K_coef',K_coef(:,:,l))
          call history_in('K_coefh',K_coefh(:,:,l))
          call history_in('sijsij',sijsij(:,:,l))
          call history_in('sijsijh',sijsijh(:,:,l))

       endif
    enddo

    !============================= (3) multiply K_coef(h) by rho ============================
    ! hereafter, K_coef(h) means rho*K
    !
!    call dbgmx('before_K_coef',K_coef(:,:,:))
!    call dbgmx('before_K_coefh',K_coefh(:,:,:))

    K_coef(:,ADM_kmin:ADM_kmax,:)=K_coef(:,ADM_kmin:ADM_kmax,:)*rho(:,ADM_kmin:ADM_kmax,:)
    do k=ADM_kmin+1,ADM_kmax
       K_coefh(:,k,:)    = K_coefh(:,k,:)*(rho(:,k,:)+rho(:,k-1,:))*0.5
    enddo
    k=ADM_kmin ! ground.  i'm not sure it is needed
    K_coefh(:,k,:)    = K_coefh(:,k,:)*rho(:,k,:)
    K_coefh(:,ADM_kmin-1,:) = 0.0d0
    K_coefh(:,ADM_kmax+1,:) = 0.0d0

    if(ADM_prc_me==ADM_prc_pl) then
       K_coef_pl(:,ADM_kmin:ADM_kmax,:)=K_coef_pl(:,ADM_kmin:ADM_kmax,:)*rho_pl(:,ADM_kmin:ADM_kmax,:)
       do k=ADM_kmin+1,ADM_kmax
          K_coefh_pl(:,k,:) = K_coefh_pl(:,k,:)*(rho_pl(:,k,:)+rho_pl(:,k-1,:))*0.5
       enddo
       k=ADM_kmin ! ground.  i'm not sure it is needed
       K_coefh_pl(:,k,:) = K_coefh_pl(:,k,:)*rho_pl(:,k,:)
       K_coefh_pl(:,ADM_kmin-1,:) = 0.0d0
       K_coefh_pl(:,ADM_kmax+1,:) = 0.0d0
    endif

!    call dbgmx('K_coef',K_coef(:,:,:))
!    call dbgmx('K_coefh',K_coefh(:,:,:))

!    call dbgmx('ee oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('ee oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))

    !============================= (4) calculate gradient3D for q and potem ============================
    ! here, i assume d rhog = 0 (i.e. nothing to do with frhog)
    !
    ! q
    do nq=1, TRC_VMAX

!       write(*,*) 'nq=',nq, 'nq+IMAX=',nq+IMAX
!       write(*,*) q(:,:,:,nq)
       if (adm_prc_me.eq.1) then
          do l=1,ADM_lall_pl
             do k=adm_kmin,ADM_kmax
!                write(*,*) ADM_prc_me, l,k,q_pl(:,k,l,nq)
             enddo
          enddo
       endif
       do l=1,ADM_lall
          do k=adm_kmin,ADM_kmax
!                write(*,*) ADM_prc_me, l,k,q(:,k,l,nq)
          enddo
       enddo

       call Gradient3dfh(        &
            var(:,:,:,GRD_XDIR,nq+IMAX),var_pl(:,:,:,GRD_XDIR,nq+IMAX),&  ! output
            var(:,:,:,GRD_YDIR,nq+IMAX),var_pl(:,:,:,GRD_YDIR,nq+IMAX),&  ! output
            var(:,:,:,GRD_ZDIR,nq+IMAX),var_pl(:,:,:,GRD_ZDIR,nq+IMAX),&  ! output
            varh(:,:,:,GRD_XDIR,nq+IMAX),varh_pl(:,:,:,GRD_XDIR,nq+IMAX),&  ! output
            varh(:,:,:,GRD_YDIR,nq+IMAX),varh_pl(:,:,:,GRD_YDIR,nq+IMAX),&  ! output
            varh(:,:,:,GRD_ZDIR,nq+IMAX),varh_pl(:,:,:,GRD_ZDIR,nq+IMAX),&  ! output
            q(:,:,:,nq), q_pl(:,:,:,nq),          & ! input
            q(:,:,:,nq), q_pl(:,:,:,nq),          & ! dummy(not used)
            input_sclh=.false.                &
            )
!       stop

    enddo

!    call dbgmx('ff oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('ff oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))

    ! theta
!       write(*,*) 'pot'
    call Gradient3dfh(        &
         var(:,:,:,GRD_XDIR,ipotem),var_pl(:,:,:,GRD_XDIR,ipotem),&  ! output
         var(:,:,:,GRD_YDIR,ipotem),var_pl(:,:,:,GRD_YDIR,ipotem),&  ! output
         var(:,:,:,GRD_ZDIR,ipotem),var_pl(:,:,:,GRD_ZDIR,ipotem),&  ! output
         varh(:,:,:,GRD_XDIR,ipotem),varh_pl(:,:,:,GRD_XDIR,ipotem),&  ! output
         varh(:,:,:,GRD_YDIR,ipotem),varh_pl(:,:,:,GRD_YDIR,ipotem),&  ! output
         varh(:,:,:,GRD_ZDIR,ipotem),varh_pl(:,:,:,GRD_ZDIR,ipotem),&  ! output
         potem(:,:,:), potem_pl(:,:,:),          & ! input
         potem(:,:,:), potem_pl(:,:,:),          & ! dummy(not used)
         input_sclh=.false.               &
          )
!    call dbgmx('gg oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('gg oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))

    !============================= (5) calculate inside the divergence and communicate ============================
    !
!    pi(:,ADM_kmin:ADM_kmax,:)=temd(:,ADM_kmin:ADM_kmax,:)/potem(:,ADM_kmin:ADM_kmax,:)
    pi(:,ADM_kmin:ADM_kmax,:)=tem(:,ADM_kmin:ADM_kmax,:)/potem(:,ADM_kmin:ADM_kmax,:)
    do k=ADM_kmin+1,ADM_kmax
       pih(:,k,:)    = (pi(:,k,:)+pi(:,k-1,:))*0.5
    enddo
    pih(:,ADM_kmin,:)=pi(:,ADM_kmin,:) ! im not sure
    !
    if(ADM_prc_me==ADM_prc_pl) then
!      pi_pl(:,ADM_kmin:ADM_kmax,:)=temd_pl(:,ADM_kmin:ADM_kmax,:)/potem_pl(:,ADM_kmin:ADM_kmax,:)
       pi_pl(:,ADM_kmin:ADM_kmax,:)=tem_pl(:,ADM_kmin:ADM_kmax,:)/potem_pl(:,ADM_kmin:ADM_kmax,:)
       do k=ADM_kmin+1,ADM_kmax
          pih_pl(:,k,:)    = (pi_pl(:,k,:)+pi_pl(:,k-1,:))*0.5
       enddo
       pih_pl(:,ADM_kmin,:)=pi_pl(:,ADM_kmin,:) ! im not sure
    endif
    !
    !
!    call dbgmx('sijvx',sij(:,:,:,grd_xdir,ivx))
!    call dbgmx('sijvy',sij(:,:,:,grd_ydir,ivx))
!    call dbgmx('sijvz',sij(:,:,:,grd_zdir,ivx))
    do idir=GRD_XDIR,GRD_ZDIR
       do ivar=IVX,IVZ
          var(:,ADM_kmin:ADM_kmax,:,idir,ivar)=Sij(:,ADM_kmin:ADM_kmax,:,idir,ivar) * 2* K_coef(:,ADM_kmin:ADM_kmax,:) * GAMMA
          varh(:,ADM_kmin:ADM_kmax,:,idir,ivar)=Sijh(:,ADM_kmin:ADM_kmax,:,idir,ivar) * 2*K_coefH(:,ADM_kmin:ADM_kmax,:) * GAMMA
       enddo
       var(:,ADM_kmin:ADM_kmax,:,idir,IPOTEM)=var(:,ADM_kmin:ADM_kmax,:,idir,IPOTEM)*K_coef(:,ADM_kmin:ADM_kmax,:)&
            / Pr * CNST_CP * pi(:,ADM_kmin:ADM_kmax,:) ! energy
       varh(:,ADM_kmin:ADM_kmax,:,idir,IPOTEM)=varh(:,ADM_kmin:ADM_kmax,:,idir,IPOTEM)*K_coefH(:,ADM_kmin:ADM_kmax,:)&
            / Pr * CNST_CP * pih(:,ADM_kmin:ADM_kmax,:) ! energy
       do ivar=IMAX+1,IMAX+TRC_VMAX
          var(:,ADM_kmin:ADM_kmax,:,idir,ivar)=var(:,ADM_kmin:ADM_kmax,:,idir,ivar)*K_coef(:,ADM_kmin:ADM_kmax,:) / Pr
          varh(:,ADM_kmin:ADM_kmax,:,idir,ivar)=varh(:,ADM_kmin:ADM_kmax,:,idir,ivar)*K_coefh(:,ADM_kmin:ADM_kmax,:) / Pr
       enddo
       if(ADM_prc_me==ADM_prc_pl) then
          do ivar=IVX,IVZ
             var_pl(:,ADM_kmin:ADM_kmax,:,idir,ivar)=Sij_pl(:,ADM_kmin:ADM_kmax,:,idir,ivar)&
                  * 2*K_coef_pl(:,ADM_kmin:ADM_kmax,:) * GAMMA
             varh_pl(:,ADM_kmin:ADM_kmax,:,idir,ivar)=Sijh_pl(:,ADM_kmin:ADM_kmax,:,idir,ivar)&
                  * 2*K_coefH_pl(:,ADM_kmin:ADM_kmax,:) * GAMMA
          enddo
          var_pl(:,ADM_kmin:ADM_kmax,:,idir,IPOTEM)=var_pl(:,ADM_kmin:ADM_kmax,:,idir,IPOTEM)*&
               K_coef_pl(:,ADM_kmin:ADM_kmax,:) / Pr * CNST_CP * pi_pl(:,ADM_kmin:ADM_kmax,:) ! energy
          varh_pl(:,ADM_kmin:ADM_kmax,:,idir,IPOTEM)=varh_pl(:,ADM_kmin:ADM_kmax,:,idir,IPOTEM)&
               *K_coefH_pl(:,ADM_kmin:ADM_kmax,:) / Pr * CNST_CP * pih_pl(:,ADM_kmin:ADM_kmax,:) ! energy
          do ivar=IMAX+1,IMAX+TRC_VMAX
             var_pl(:,ADM_kmin:ADM_kmax,:,idir,ivar)=var_pl(:,ADM_kmin:ADM_kmax,:,idir,ivar)&
                  *K_coef_pl(:,ADM_kmin:ADM_kmax,:) / Pr
             varh_pl(:,ADM_kmin:ADM_kmax,:,idir,ivar)=varh_pl(:,ADM_kmin:ADM_kmax,:,idir,ivar)&
                  *K_coefh_pl(:,ADM_kmin:ADM_kmax,:) / Pr
          enddo
       endif
    enddo

!    write(*,*) reshape(var(:,:,:,:,:), (/ADM_gall,ADM_kall,ADM_lall,3*(IMAX+TRC_VMAX)/)) , &
!         reshape(var_pl(:,:,:,:,:),(/ADM_gall_pl,ADM_kall,ADM_lall_pl,3*(IMAX+TRC_VMAX)/)) 

!    call comm_data_transfer(var(:,:,:,:,1),var_pl(:,:,:,:,1))
!    call comm_data_transfer(reshape(var(:,:,:,:,:), (/ADM_gall,ADM_kall,ADM_lall,3*(IMAX+TRC_VMAX)/)) , &
!         reshape(var_pl(:,:,:,:,:),(/ADM_gall_pl,ADM_kall,ADM_lall_pl,3*(IMAX+TRC_VMAX)/)) )
!    call comm_data_transfer(reshape(var, (/size(var,1),size(var,2),size(var,3),3*(IMAX+TRC_VMAX) /)), &
!         reshape(var_pl, (/size(var_pl,1),size(var_pl,2),size(var_pl,3),3*(IMAX+TRC_VMAX) /) ))
!    call comm_data_transfer(reshape(varh, (/size(varh,1),size(varh,2),size(varh,3),3*(IMAX+TRC_VMAX) /)), &
!         reshape(varh_pl, (/size(varh_pl,1),size(varh_pl,2),size(varh_pl,3),3*(IMAX+TRC_VMAX) /) ))

!    write(*,*) shape(var)
    do ivar=IVX,IMAX+TRC_VMAX
!       call dbgmx('vxbefore'//char(ivar+48),var(:,:,:,grd_xdir,ivar))
!       call dbgmx('vxbefore_pl'//char(ivar+48),var_pl(:,:,:,grd_xdir,ivar))

       call comm_data_transfer(var(:,:,:,:,ivar),var_pl(:,:,:,:,ivar))
       call comm_data_transfer(varh(:,:,:,:,ivar),varh_pl(:,:,:,:,ivar))

!       call dbgmx('vxafter'//char(ivar+48),var(:,:,:,grd_xdir,ivar))
!       call dbgmx('vxafter_pl'//char(ivar+48),var_pl(:,:,:,grd_xdir,ivar))
    enddo

    !============================= (6) calculate divergence  ============================
    grhogvx=0
    grhogvy=0
    grhogvz=0
    grhogvx_pl=0
    grhogvy_pl=0
    grhogvz_pl=0
    grhogq=0
    grhoge=0

    ! velocity
!write(*,*)'a',ADM_prc_me
    call gsqrt_Div3dfh(  &
         var(:,:,:,GRD_XDIR,IVX), var_pl(:,:,:,GRD_XDIR,IVX),&
         var(:,:,:,GRD_YDIR,IVX), var_pl(:,:,:,GRD_YDIR,IVX),&
         var(:,:,:,GRD_ZDIR,IVX), var_pl(:,:,:,GRD_ZDIR,IVX),&
         varh(:,:,:,GRD_XDIR,IVX), varh_pl(:,:,:,GRD_XDIR,IVX),&
         varh(:,:,:,GRD_YDIR,IVX), varh_pl(:,:,:,GRD_YDIR,IVX),&
         varh(:,:,:,GRD_ZDIR,IVX), varh_pl(:,:,:,GRD_ZDIR,IVX),&
         grhogvx, grhogvx_pl,                  & !out (full level velocity)
         grhogvxh, grhogvxh_pl,                  & !out (half level velocity)
         output_sclh=.true.                & !in optional ( input half level scl, or not )
       )

!write(*,*)'b',ADM_prc_me
    call gsqrt_Div3dfh(  &
         var(:,:,:,GRD_XDIR,IVY), var_pl(:,:,:,GRD_XDIR,IVY),&
         var(:,:,:,GRD_YDIR,IVY), var_pl(:,:,:,GRD_YDIR,IVY),&
         var(:,:,:,GRD_ZDIR,IVY), var_pl(:,:,:,GRD_ZDIR,IVY),&
         varh(:,:,:,GRD_XDIR,IVY), varh_pl(:,:,:,GRD_XDIR,IVY),&
         varh(:,:,:,GRD_YDIR,IVY), varh_pl(:,:,:,GRD_YDIR,IVY),&
         varh(:,:,:,GRD_ZDIR,IVY), varh_pl(:,:,:,GRD_ZDIR,IVY),&
         grhogvy, grhogvy_pl,                  & !out (full level velocity)
         grhogvyh, grhogvyh_pl,                  & !out (half level velocity)
         output_sclh=.true.                & !in optional ( input half level scl, or not )
       )

!write(*,*)'c',ADM_prc_me
    call gsqrt_Div3dfh(  &
         var(:,:,:,GRD_XDIR,IVZ), var_pl(:,:,:,GRD_XDIR,IVZ),&
         var(:,:,:,GRD_YDIR,IVZ), var_pl(:,:,:,GRD_YDIR,IVZ),&
         var(:,:,:,GRD_ZDIR,IVZ), var_pl(:,:,:,GRD_ZDIR,IVZ),&
         varh(:,:,:,GRD_XDIR,IVZ), varh_pl(:,:,:,GRD_XDIR,IVZ),&
         varh(:,:,:,GRD_YDIR,IVZ), varh_pl(:,:,:,GRD_YDIR,IVZ),&
         varh(:,:,:,GRD_ZDIR,IVZ), varh_pl(:,:,:,GRD_ZDIR,IVZ),&
         grhogvz, grhogvz_pl,                  & !out (full level velocity)
         grhogvzh, grhogvzh_pl,                  & !out (half level velocity)
         output_sclh=.true.                & !in optional ( input half level scl, or not )
       )

    ! vx,vy,vz
    call OPRT_horizontalize_vec(&
         grhogvx, grhogvx_pl,   & !--- inout
         grhogvy, grhogvy_pl,   & !--- inout
         grhogvz, grhogvz_pl)     !--- inout

    frhogvx(:,ADM_kmin:ADM_kmax,:)=frhogvx(:,ADM_kmin:ADM_kmax,:)+grhogvx(:,ADM_kmin:ADM_kmax,:)
    frhogvy(:,ADM_kmin:ADM_kmax,:)=frhogvy(:,ADM_kmin:ADM_kmax,:)+grhogvy(:,ADM_kmin:ADM_kmax,:)
    frhogvz(:,ADM_kmin:ADM_kmax,:)=frhogvz(:,ADM_kmin:ADM_kmax,:)+grhogvz(:,ADM_kmin:ADM_kmax,:)
    frhogvx_pl(:,ADM_kmin:ADM_kmax,:)=frhogvx_pl(:,ADM_kmin:ADM_kmax,:)+grhogvx_pl(:,ADM_kmin:ADM_kmax,:)
    frhogvy_pl(:,ADM_kmin:ADM_kmax,:)=frhogvy_pl(:,ADM_kmin:ADM_kmax,:)+grhogvy_pl(:,ADM_kmin:ADM_kmax,:)
    frhogvz_pl(:,ADM_kmin:ADM_kmax,:)=frhogvz_pl(:,ADM_kmin:ADM_kmax,:)+grhogvz_pl(:,ADM_kmin:ADM_kmax,:)

    !w
    do k=ADM_kmin+1,ADM_kmax
       frhogw(:,k,:) = frhogw(:,k,:) + ( &
            grhogvxh(:,k,:) * GRD_x(:,ADM_knone,:,GRD_XDIR)+&
            grhogvyh(:,k,:) * GRD_x(:,ADM_knone,:,GRD_YDIR)+&
            grhogvzh(:,k,:) * GRD_x(:,ADM_knone,:,GRD_ZDIR)  ) / GRD_rscale
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin+1,ADM_kmax
          frhogw_pl(:,k,:) = frhogw_pl(:,k,:) + ( &
               grhogvxh_pl(:,k,:) * GRD_x_pl(:,ADM_knone,:,GRD_XDIR)+&
               grhogvyh_pl(:,k,:) * GRD_x_pl(:,ADM_knone,:,GRD_YDIR)+&
               grhogvzh_pl(:,k,:) * GRD_x_pl(:,ADM_knone,:,GRD_ZDIR)  ) / GRD_rscale
       enddo
    endif

    ! q
    do nq=1, TRC_VMAX

       call gsqrt_Div3dfh(  &
         var(:,:,:,GRD_XDIR,IMAX+nq), var_pl(:,:,:,GRD_XDIR,IMAX+nq),&
         var(:,:,:,GRD_YDIR,IMAX+nq), var_pl(:,:,:,GRD_YDIR,IMAX+nq),&
         var(:,:,:,GRD_ZDIR,IMAX+nq), var_pl(:,:,:,GRD_ZDIR,IMAX+nq),&
         varh(:,:,:,GRD_XDIR,IMAX+nq), varh_pl(:,:,:,GRD_XDIR,IMAX+nq),&
         varh(:,:,:,GRD_YDIR,IMAX+nq), varh_pl(:,:,:,GRD_YDIR,IMAX+nq),&
         varh(:,:,:,GRD_ZDIR,IMAX+nq), varh_pl(:,:,:,GRD_ZDIR,IMAX+nq),&
         grhogq(:,:,:,nq), grhogq_pl(:,:,:,nq),                  & !out (full level velocity)
         dummy,dummy_pl,                      & !dummy
         output_sclh=.false.                & !in optional ( input half level scl, or not )
       )

    enddo
    frhogq(:,ADM_kmin:ADM_kmax,:,:)=frhogq(:,ADM_kmin:ADM_kmax,:,:)+grhogq(:,ADM_kmin:ADM_kmax,:,:)

    ! energy
    call gsqrt_Div3dfh(  &
         var(:,:,:,GRD_XDIR,ipotem), var_pl(:,:,:,GRD_XDIR,ipotem),&
         var(:,:,:,GRD_YDIR,ipotem), var_pl(:,:,:,GRD_YDIR,ipotem),&
         var(:,:,:,GRD_ZDIR,ipotem), var_pl(:,:,:,GRD_ZDIR,ipotem),&
         varh(:,:,:,GRD_XDIR,ipotem), varh_pl(:,:,:,GRD_XDIR,ipotem),&
         varh(:,:,:,GRD_YDIR,ipotem), varh_pl(:,:,:,GRD_YDIR,ipotem),&
         varh(:,:,:,GRD_ZDIR,ipotem), varh_pl(:,:,:,GRD_ZDIR,ipotem),&
         grhoge(:,:,:), grhoge_pl(:,:,:),                  & !out (full level)
         dummy,dummy_pl,                      & !dummy
         output_sclh=.false.                & !in optional ( input half level scl, or not )
       )

    frhoge(:,ADM_kmin:ADM_kmax,:)=frhoge(:,ADM_kmin:ADM_kmax,:)+grhoge(:,ADM_kmin:ADM_kmax,:)

    ! I don't know about grhogetot. more precisely, change in kinetic energy sould be concidered.
    frhogetot(:,ADM_kmin:ADM_kmax,:)=frhogetot(:,ADM_kmin:ADM_kmax,:)+grhoge(:,ADM_kmin:ADM_kmax,:)

    do l=1,ADM_lall
       if (dbgfirst)then

          call history_in('gvx',grhogvx(:,:,l))
          call history_in('gvy',grhogvy(:,:,l))
          call history_in('gvz',grhogvz(:,:,l))
          call history_in('gvxh',grhogvxh(:,:,l))
          call history_in('gvyh',grhogvyh(:,:,l))
          call history_in('gvzh',grhogvzh(:,:,l))
          call history_in('gtem',grhoge(:,:,l))
          call history_in('gq1',grhogq(:,:,l,1))
          call history_in('var11',var(:,:,l,1,1))
          call history_in('var21',var(:,:,l,2,1))
          call history_in('var31',var(:,:,l,3,1))
          call history_in('var11h',varh(:,:,l,1,1))
          call history_in('var21h',varh(:,:,l,2,1))
          call history_in('var31h',varh(:,:,l,3,1))

       endif
    enddo
!    dbgfirst=.false.

    return
  end subroutine tb_smg_driver
  !-------------------------------------------------------------------------------
  subroutine tb_smg_oprt_init(  &
       )
  
    use mod_grd, only: &
         GRD_xdir,     &
         GRD_ydir,     &
         GRD_zdir,     &
         GRD_rdgz,     &!101201
         GRD_rdgzh,    &!
         GRD_vz,       &!101201
         GRD_Z,        &!101201
         GRD_x,        &
         GRD_x_pl,     &
         GRD_rscale!=cnst_eradisu
  
    use mod_runconf, only :         &
         TRC_VMAX
    implicit none
    integer::k
                                                                                       

!    write(*,*) GRD_x/GRD_rscale


    if(first)then  ! --> why is it needed?
       allocate(smg_oprt_cxh(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_cxh_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_cyh(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_cyh_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_czh(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_czh_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_cx(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_cx_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_cy(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_cy_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_cz(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_cz_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_rgamH(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_rgamH_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_rgam(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_rgam_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
!       allocate(smg_oprt_GzGz(ADM_gall,ADM_kall,ADM_lall))
!       allocate(smg_oprt_GzGz_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
!       allocate(smg_oprt_GzGzh(ADM_gall,ADM_kall,ADM_lall))
!       allocate(smg_oprt_GzGzh_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_gsqrt(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_gsqrt_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_gsqrtH(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_gsqrtH_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_GAM(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_GAM_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))
       allocate(smg_oprt_GAMH(ADM_gall,ADM_kall,ADM_lall))
       allocate(smg_oprt_GAMH_pl(ADM_gall_PL,ADM_kall,ADM_lall_PL))

       allocate(var(ADM_gall,ADM_kall,ADM_lall,GRD_XDIR:GRD_ZDIR,IMAX+TRC_VMAX))
       allocate(varh(ADM_gall,ADM_kall,ADM_lall,GRD_XDIR:GRD_ZDIR,IMAX+TRC_VMAX))
       allocate(var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_XDIR:GRD_ZDIR,IMAX+TRC_VMAX))
       allocate(varh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_XDIR:GRD_ZDIR,IMAX+TRC_VMAX))
       first=.false.
    endif



    ! ad hoc
    smg_oprt_gsqrt(:,ADM_kmin:ADM_kmax,:) = VMTR_GSGAM2(:,ADM_kmin:ADM_kmax,:)/VMTR_GAM2(:,ADM_kmin:ADM_kmax,:)
    smg_oprt_gsqrtH(:,ADM_kmin:ADM_kmax,:) = VMTR_GSGAM2H(:,ADM_kmin:ADM_kmax,:)/VMTR_GAM2H(:,ADM_kmin:ADM_kmax,:)
    smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:) = sqrt(VMTR_GAM2(:,ADM_kmin:ADM_kmax,:))
    smg_oprt_GAMH(:,ADM_kmin:ADM_kmax,:) = sqrt(VMTR_GAM2H(:,ADM_kmin:ADM_kmax,:))

    if(ADM_prc_me==ADM_prc_pl) then
       smg_oprt_gsqrt_pl(:,ADM_kmin:ADM_kmax,:) = VMTR_GSGAM2_pl(:,ADM_kmin:ADM_kmax,:)/VMTR_GAM2_pl(:,ADM_kmin:ADM_kmax,:)
       smg_oprt_gsqrtH_pl(:,ADM_kmin:ADM_kmax,:) = VMTR_GSGAM2H_pl(:,ADM_kmin:ADM_kmax,:)/VMTR_GAM2H_pl(:,ADM_kmin:ADM_kmax,:)
       smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:) = sqrt(VMTR_GAM2_pl(:,ADM_kmin:ADM_kmax,:))
       smg_oprt_GAMH_pl(:,ADM_kmin:ADM_kmax,:) = sqrt(VMTR_GAM2H_pl(:,ADM_kmin:ADM_kmax,:))
    endif

    
!    do k=ADM_kmin,ADM_kmax
!       if (ADM_prc_me.eq.1) write(*,*) k,smg_oprt_GAM_pl(:,k,:),VMTR_GAM2_pl(:,k,:)
!    enddo
!    stop

    ! 1/gamma at half level
    smg_oprt_rgamH(:,ADM_kmin:ADM_kmax,:)= 1.0d0/sqrt(VMTR_GAM2H(:,ADM_kmin:ADM_kmax,:))
    smg_oprt_rgam(:,ADM_kmin:ADM_kmax,:)= 1.0d0/sqrt(VMTR_GAM2(:,ADM_kmin:ADM_kmax,:))
    if(ADM_prc_me==ADM_prc_pl) then
      smg_oprt_rgamH_pl(:,ADM_kmin:ADM_kmax,:)= 1.0d0/sqrt(VMTR_GAM2H_pl(:,ADM_kmin:ADM_kmax,:))
      smg_oprt_rgam_pl(:,ADM_kmin:ADM_kmax,:)= 1.0d0/sqrt(VMTR_GAM2_pl(:,ADM_kmin:ADM_kmax,:))
    endif

!    call dbgmx('smg_oprt_rgam',smg_oprt_rgam)
!    call dbgmx('smg_oprt_rgamH',smg_oprt_rgamH)


    ! coef for grad3d
    !half
    do k=ADM_kmin,ADM_kmax
       smg_oprt_cxh(:,k,:) =  GRD_rdgz(k)*   (GRD_x(:,ADM_knone,:,GRD_XDIR) &
         *VMTR_RGSH(:,k,:)/GRD_rscale  + VMTR_GZXH(:,k,:)*smg_oprt_rgamH(:,k,:) )
       smg_oprt_cyh(:,k,:) =  GRD_rdgz(k)*   (GRD_x(:,ADM_knone,:,GRD_YDIR) &
         *VMTR_RGSH(:,k,:)/GRD_rscale + VMTR_GZYH(:,k,:)*smg_oprt_rgamH(:,k,:) )
       smg_oprt_czh(:,k,:) =  GRD_rdgz(k)*   (GRD_x(:,ADM_knone,:,GRD_ZDIR) &
         *VMTR_RGSH(:,k,:)/GRD_rscale + VMTR_GZZH(:,k,:)*smg_oprt_rgamH(:,k,:) )
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
      do k=ADM_kmin,ADM_kmax
        smg_oprt_cxh_pl(:,k,:) =  GRD_rdgz(k)*  (GRD_x_pl(:,ADM_knone,:,GRD_XDIR) &
         *VMTR_RGSH_pl(:,k,:)/GRD_rscale + VMTR_GZXH_pl(:,k,:)*smg_oprt_rgamH_pl(:,k,:)  )
        smg_oprt_cyh_pl(:,k,:) =  GRD_rdgz(k)*  (GRD_x_pl(:,ADM_knone,:,GRD_YDIR) &
         *VMTR_RGSH_pl(:,k,:)/GRD_rscale + VMTR_GZYH_pl(:,k,:)*smg_oprt_rgamH_pl(:,k,:) )
        smg_oprt_czh_pl(:,k,:) =  GRD_rdgz(k)*  (GRD_x_pl(:,ADM_knone,:,GRD_ZDIR)  &
         *VMTR_RGSH_pl(:,k,:)/GRD_rscale + VMTR_GZZH_pl(:,k,:)*smg_oprt_rgamH_pl(:,k,:) )
      enddo
    endif

    !full
    do k=ADM_kmin,ADM_kmax
       smg_oprt_cx(:,k,:) =  GRD_rdgzh(k)*   (GRD_x(:,ADM_knone,:,GRD_XDIR) &
         /smg_oprt_gsqrt(:,k,:)/GRD_rscale  + VMTR_GZX(:,k,:)*smg_oprt_rgam(:,k,:) )
       smg_oprt_cy(:,k,:) =  GRD_rdgzh(k)*   (GRD_x(:,ADM_knone,:,GRD_YDIR) &
         /smg_oprt_gsqrt(:,k,:)/GRD_rscale + VMTR_GZY(:,k,:)*smg_oprt_rgam(:,k,:) )
       smg_oprt_cz(:,k,:) =  GRD_rdgzh(k)*   (GRD_x(:,ADM_knone,:,GRD_ZDIR) &
         /smg_oprt_gsqrt(:,k,:)/GRD_rscale + VMTR_GZZ(:,k,:)*smg_oprt_rgam(:,k,:) )
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
      do k=ADM_kmin,ADM_kmax
        smg_oprt_cx_pl(:,k,:) =  GRD_rdgzh(k)*  (GRD_x_pl(:,ADM_knone,:,GRD_XDIR) &
         /smg_oprt_gsqrt_pl(:,k,:)/GRD_rscale + VMTR_GZX_pl(:,k,:)*smg_oprt_rgam_pl(:,k,:)  )
        smg_oprt_cy_pl(:,k,:) =  GRD_rdgzh(k)*  (GRD_x_pl(:,ADM_knone,:,GRD_YDIR) &
         /smg_oprt_gsqrt_pl(:,k,:)/GRD_rscale + VMTR_GZY_pl(:,k,:)*smg_oprt_rgam_pl(:,k,:) )
        smg_oprt_cz_pl(:,k,:) =  GRD_rdgzh(k)*  (GRD_x_pl(:,ADM_knone,:,GRD_ZDIR)  &
         /smg_oprt_gsqrt_pl(:,k,:)/GRD_rscale + VMTR_GZZ_pl(:,k,:)*smg_oprt_rgam_pl(:,k,:) )
      enddo
    endif




!    call dbgmx('oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))

!!$!    ! 1/gamma at half level
!!$!    smg_oprt_rgamH(:,ADM_kmin:ADM_kmax,:)= 1d0/sqrt(VMTR_GAM2H(:,ADM_kmin:ADM_kmax,:))
!!$!    smg_oprt_rgamH_pl(:,ADM_kmin:ADM_kmax,:)= 1d0/sqrt(VMTR_GAM2H_pl(:,ADM_kmin:ADM_kmax,:))
!!$
!!$
!!$!koko
!!$!    smg_oprt_GzGzh(:,:,:)= VMTR_GXXH(:,:,:)**2+VMTR_GYXH(:,:,:)**2+VMTR_GZXH(:,:,:)**2
!!$!    smg_oprt_GzGzh_pl(:,:,:)= VMTR_GXXH_pl(:,:,:)**2+VMTR_GYXH_pl(:,:,:)**2+VMTR_GZXH_pl(:,:,:)**2
!!$!    smg_oprt_GzGz(:,:,:)= VMTR_GXX(:,:,:)**2+VMTR_GYX(:,:,:)**2+VMTR_GZX(:,:,:)**2
!!$!    smg_oprt_GzGz_pl(:,:,:)= VMTR_GXX_pl(:,:,:)**2+VMTR_GYX_pl(:,:,:)**2+VMTR_GZX_pl(:,:,:)**2
!!$
!!$    smg_oprt_GzGzh(:,:,:)= VMTR_GZXH(:,:,:)**2+VMTR_GZYH(:,:,:)**2+VMTR_GZZH(:,:,:)**2
!!$    smg_oprt_GzGz(:,:,:)= VMTR_GZX(:,:,:)**2+VMTR_GZY(:,:,:)**2+VMTR_GZZ(:,:,:)**2
!!$    if(ADM_prc_me==ADM_prc_pl) then
!!$      smg_oprt_GzGzh_pl(:,:,:)= VMTR_GZXH_pl(:,:,:)**2+VMTR_GZYH_pl(:,:,:)**2+VMTR_GZZH_pl(:,:,:)**2
!!$      smg_oprt_GzGz_pl(:,:,:)= VMTR_GZX_pl(:,:,:)**2+VMTR_GZY_pl(:,:,:)**2+VMTR_GZZ_pl(:,:,:)**2
!!$    endif





    return
  end subroutine tb_smg_oprt_init
  !-----------------------------------------------------------------------------------------
  ! gradient3d hybrid routine (full and half)
  subroutine Gradient3dfh(  &
       vx, vx_pl,                  & !out (full level tensor component)
       vy, vy_pl,                  & !out (full level tensor component)
       vz, vz_pl,                  & !out (full level tensor component)
       vxh, vxh_pl,                  & !out (half level tensor component)
       vyh, vyh_pl,                  & !out (half level tensor component)
       vzh, vzh_pl,                  & !out (half level tensor component)
       scl, scl_pl,                  & !in (full level velocity)
       sclh, sclh_pl,              & !in (half level velocity)
       input_sclh                & !in optional ( input half level scl, or not )
       )
    ! note that vx, vy, vz are not velocity but just vector.
  
    use mod_oprt, only:&
         OPRT_gradient
  
    implicit none
    real(8),intent(in)::sclh(ADM_gall,ADM_kall,ADM_lall)          ! half level
    real(8),intent(in)::sclh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL) ! half level
    real(8),intent(in)::scl(ADM_gall,ADM_kall,ADM_lall)            ! full level
    real(8),intent(in)::scl_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)   ! full level
    !
    ! horizontal gradient
    real(8),intent(out)::vx(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(out)::vy(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(out)::vz(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(out)::vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    real(8),intent(out)::vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    real(8),intent(out)::vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    !
    real(8),intent(out)::vxh(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(out)::vyh(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(out)::vzh(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(out)::vxh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    real(8),intent(out)::vyh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    real(8),intent(out)::vzh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    !
    logical,intent(in),optional::input_sclh
    logical ::input_sclh_in!=.true.

    integer::k,i,l

    if (.not.present(input_sclh)) input_sclh_in=.true.
    if (present(input_sclh)) input_sclh_in=input_sclh
  
    !initial
    vx=0 ! needed because halo region is not written in oprt_gradient
    vy=0 ! needed because halo region is not written in oprt_gradient
    vz=0 ! needed because halo region is not written in oprt_gradient
    vx_pl=0 ! needed because halo region is not written in oprt_gradient
    vy_pl=0 ! needed because halo region is not written in oprt_gradient
    vz_pl=0 ! needed because halo region is not written in oprt_gradient
    vxh=0 ! needed because halo region is not written in oprt_gradient
    vyh=0 ! needed because halo region is not written in oprt_gradient
    vzh=0 ! needed because halo region is not written in oprt_gradient
    vxh_pl=0 ! needed because halo region is not written in oprt_gradient
    vyh_pl=0 ! needed because halo region is not written in oprt_gradient
    vzh_pl=0 ! needed because halo region is not written in oprt_gradient


    ! horizontal gradient for full
    call OPRT_gradient(&
       vx, vx_pl,                  &
       vy, vy_pl,                  &
       vz, vz_pl,                  &
       scl, scl_pl              &
       )
    ! horizontal gradient for half
    if (input_sclh_in) then
!       write(*,*) 'input!'
       call OPRT_gradient(&
            vxh, vxh_pl,                  &
            vyh, vyh_pl,                  &
            vzh, vzh_pl,                  &
            sclh, sclh_pl              &
            )
    else
       do k=ADM_kmin+1,ADM_kmax
          vxh(:,k,:)=(vx(:,k,:)+vx(:,k-1,:))/2
          vyh(:,k,:)=(vy(:,k,:)+vy(:,k-1,:))/2
          vzh(:,k,:)=(vz(:,k,:)+vz(:,k-1,:))/2
       enddo
       k=ADM_kmin
       vxh(:,k,:)=(vx(:,k,:))! not used effectively
       vyh(:,k,:)=(vy(:,k,:))! not used effectively
       vzh(:,k,:)=(vz(:,k,:))! not used effectively
       if(ADM_prc_me==ADM_prc_pl) then
          do k=ADM_kmin+1,ADM_kmax
             vxh_pl(:,k,:)=(vx_pl(:,k,:)+vx_pl(:,k-1,:))/2
             vyh_pl(:,k,:)=(vy_pl(:,k,:)+vy_pl(:,k-1,:))/2
             vzh_pl(:,k,:)=(vz_pl(:,k,:)+vz_pl(:,k-1,:))/2
          enddo
          k=ADM_kmin
          vxh_pl(:,k,:)=(vx_pl(:,k,:))! not used effectively
          vyh_pl(:,k,:)=(vy_pl(:,k,:))! not used effectively
          vzh_pl(:,k,:)=(vz_pl(:,k,:))! not used effectively
       endif
    endif


    !full level ---------------------------------------------------------------------
    vx(:,ADM_kmin:ADM_kmax,:)=vx(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgam(:,ADM_kmin:ADM_kmax,:)
    vy(:,ADM_kmin:ADM_kmax,:)=vy(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgam(:,ADM_kmin:ADM_kmax,:)
    vz(:,ADM_kmin:ADM_kmax,:)=vz(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgam(:,ADM_kmin:ADM_kmax,:)


    do l=1,ADM_lall
           if (input_sclh_in)call history_in('scl',scl(:,:,l))
           if (input_sclh_in)call history_in('hogevx',vx(:,:,l))
    enddo

    do k=ADM_kmin,ADM_kmax-1
       ! add vertical gradient
!       if (ADM_prc_me.eq.1) write(*,*) 'fb',k, vx(200,k,1)
       vx(:,k,:)=vx(:,k,:)+(sclh(:,k+1,:)-sclh(:,k,:)) * smg_oprt_cx(:,k,:)
       vy(:,k,:)=vy(:,k,:)+(sclh(:,k+1,:)-sclh(:,k,:)) * smg_oprt_cy(:,k,:)
       vz(:,k,:)=vz(:,k,:)+(sclh(:,k+1,:)-sclh(:,k,:)) * smg_oprt_cz(:,k,:)
!       if (ADM_prc_me.eq.1) write(*,*) 'f',k, (sclh(200,k+1,1)-sclh(200,k,1)) * smg_oprt_cx(200,k,1),vx(200,k,1)
    enddo
    k=ADM_kmax ! vertical gradient at k=kmax is ignored
    vx(:,k,:)=vx(:,k,:)!+(sclh(:,k,:)-0) * smg_oprt_cx(:,k,:)
    vy(:,k,:)=vy(:,k,:)!+(sclh(:,k,:)-0) * smg_oprt_cy(:,k,:)
    vz(:,k,:)=vz(:,k,:)!+(sclh(:,k,:)-0) * smg_oprt_cz(:,k,:)


!    do l=1,ADM_lall
!           if (input_sclh_in)call history_in('dscl',(cshift(sclh(:,:,l),1,2)-sclh(:,:,l)) * smg_oprt_cx(:,k,:))
!    enddo


    !
    if(ADM_prc_me==ADM_prc_pl) then
       vx_pl(:,ADM_kmin:ADM_kmax,:)=vx_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgam_pl(:,ADM_kmin:ADM_kmax,:)
       vy_pl(:,ADM_kmin:ADM_kmax,:)=vy_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgam_pl(:,ADM_kmin:ADM_kmax,:)
       vz_pl(:,ADM_kmin:ADM_kmax,:)=vz_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgam_pl(:,ADM_kmin:ADM_kmax,:)
       do k=ADM_kmin,ADM_kmax-1
          vx_pl(:,k,:)=vx_pl(:,k,:)+(sclh_pl(:,k+1,:)-sclh_pl(:,k,:)) * smg_oprt_cx_pl(:,k,:)
          vy_pl(:,k,:)=vy_pl(:,k,:)+(sclh_pl(:,k+1,:)-sclh_pl(:,k,:)) * smg_oprt_cy_pl(:,k,:)
          vz_pl(:,k,:)=vz_pl(:,k,:)+(sclh_pl(:,k+1,:)-sclh_pl(:,k,:)) * smg_oprt_cz_pl(:,k,:)
       enddo
       k=ADM_kmax  ! vertical gradient at k=kmax is ignored
       vx_pl(:,k,:)=vx_pl(:,k,:)!+(sclh_pl(:,k,:)-0) * smg_oprt_cx_pl(:,k,:)
       vy_pl(:,k,:)=vy_pl(:,k,:)!+(sclh_pl(:,k,:)-0) * smg_oprt_cy_pl(:,k,:)
       vz_pl(:,k,:)=vz_pl(:,k,:)!+(sclh_pl(:,k,:)-0) * smg_oprt_cz_pl(:,k,:)
    endif


    !half level ------------------------------------------------------------------
    !
    vxh(:,ADM_kmin:ADM_kmax,:)=vxh(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgamH(:,ADM_kmin:ADM_kmax,:)
    vyh(:,ADM_kmin:ADM_kmax,:)=vyh(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgamH(:,ADM_kmin:ADM_kmax,:)
    vzh(:,ADM_kmin:ADM_kmax,:)=vzh(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgamH(:,ADM_kmin:ADM_kmax,:)
    do l=1,ADM_lall
           if (input_sclh_in)call history_in('sclh',sclh(:,:,l))
           if (input_sclh_in)call history_in('hogevxh',vxh(:,:,l))
    enddo

    do k=ADM_kmin+1,ADM_kmax
       ! add vertical gradient
!       if (ADM_prc_me.eq.1) write(*,*) 'hb',k, vxh(200,k,1)
       vxh(:,k,:)=vxh(:,k,:)+(scl(:,k,:)-scl(:,k-1,:)) * smg_oprt_cxh(:,k,:)
       vyh(:,k,:)=vyh(:,k,:)+(scl(:,k,:)-scl(:,k-1,:)) * smg_oprt_cyh(:,k,:)
       vzh(:,k,:)=vzh(:,k,:)+(scl(:,k,:)-scl(:,k-1,:)) * smg_oprt_czh(:,k,:)
!       if (ADM_prc_me.eq.1) write(*,*) 'h',k, (scl(200,k+1,1)-scl(200,k,1)) * smg_oprt_cx(200,k,1),vxh(200,k,1)
    enddo
    k=ADM_kmin ! ground level (not used effectively)
    vxh(:,k,:)=vxh(:,k,:)!+(scl(:,k,:)-0) * smg_oprt_cxh(:,k,:)
    vyh(:,k,:)=vyh(:,k,:)!+(scl(:,k,:)-0) * smg_oprt_cyh(:,k,:)
    vzh(:,k,:)=vzh(:,k,:)!+(scl(:,k,:)-0) * smg_oprt_czh(:,k,:)
    !
    do l=1,ADM_lall
!           if (input_sclh_in)call history_in('dsclh',(cshift(scl(:,:,l),1,2)-scl(:,:,l)) * smg_oprt_cxh(:,k,:))
    enddo

    if(ADM_prc_me==ADM_prc_pl) then
       vxh_pl(:,ADM_kmin:ADM_kmax,:)=vxh_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgamH_pl(:,ADM_kmin:ADM_kmax,:)
       vyh_pl(:,ADM_kmin:ADM_kmax,:)=vyh_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgamH_pl(:,ADM_kmin:ADM_kmax,:)
       vzh_pl(:,ADM_kmin:ADM_kmax,:)=vzh_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_rgamH_pl(:,ADM_kmin:ADM_kmax,:)
       do k=ADM_kmin+1,ADM_kmax
          vxh_pl(:,k,:)=vxh_pl(:,k,:)+(scl_pl(:,k,:)-scl_pl(:,k-1,:)) * smg_oprt_cxh_pl(:,k,:)
          vyh_pl(:,k,:)=vyh_pl(:,k,:)+(scl_pl(:,k,:)-scl_pl(:,k-1,:)) * smg_oprt_cyh_pl(:,k,:)
          vzh_pl(:,k,:)=vzh_pl(:,k,:)+(scl_pl(:,k,:)-scl_pl(:,k-1,:)) * smg_oprt_czh_pl(:,k,:)
       enddo
       k=ADM_kmin ! ground level (not used effectively)
       vxh_pl(:,k,:)=vxh_pl(:,k,:)!+(scl_pl(:,k,:)-0) * smg_oprt_cxh_pl(:,k,:)
       vyh_pl(:,k,:)=vyh_pl(:,k,:)!+(scl_pl(:,k,:)-0) * smg_oprt_cyh_pl(:,k,:)
       vzh_pl(:,k,:)=vzh_pl(:,k,:)!+(scl_pl(:,k,:)-0) * smg_oprt_czh_pl(:,k,:)
    endif

!    call dbgmx('grad3d_vx',vx)
!    call dbgmx('grad3d_vxh',vxh)
!    call dbgmx('grad3d_vy',vy)
!    call dbgmx('grad3d_vyh',vyh)
!    call dbgmx('grad3d_vz',vz)
!    call dbgmx('grad3d_vzh',vzh)
!    call dbgmx('grad3d_vx_pl',vx_pl)

    return
  end subroutine Gradient3dfh
  !-------------------------------------------------------------------------------
  ! gsqrt * gam^2 * div3d  (full and half)
  subroutine gsqrt_Div3dfh(  &
       vx, vx_pl,                  & !in (full level) tensor component)
       vy, vy_pl,                  & !in (full level)
       vz, vz_pl,                  & !in (full level)
       vxh, vxh_pl,                  & !in (half level)
       vyh, vyh_pl,                  & !in (half level)
       vzh, vzh_pl,                  & !in (half level)
       scl, scl_pl,                  & !out (full level velocity)
       sclh, sclh_pl,              & !out (half level velocity)
       output_sclh                & !optional, in ( input half level scl, or not )
       )
    ! note that vx, vy, vz are not velocity but just vector.
  
    use mod_oprt, only:&
         OPRT_divergence
    use mod_grd,only:&
         GRD_RDGZ,&
         GRD_RDGZH,&
         GRD_x,&
         GRD_x_pl,&
         GRD_rscale,&
         GRD_XDIR,&
         GRD_YDIR,&
         GRD_ZDIR
!    use mod_vmtr,only:&
!         vmtr_gam,&
!         vmtr_gam_pl

    implicit none
    real(8),intent(out)::sclh(ADM_gall,ADM_kall,ADM_lall)          ! half level
    real(8),intent(out)::sclh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL) ! half level
    real(8),intent(out)::scl(ADM_gall,ADM_kall,ADM_lall)            ! full level
    real(8),intent(out)::scl_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)   ! full level
    !
    ! horizontal gradient
    real(8),intent(in)::vx(ADM_gall,ADM_kall,ADM_lall)           ! full level
    real(8),intent(in)::vy(ADM_gall,ADM_kall,ADM_lall)           ! full level
    real(8),intent(in)::vz(ADM_gall,ADM_kall,ADM_lall)           ! full level
    real(8),intent(in)::vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! full level
    real(8),intent(in)::vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! full level
    real(8),intent(in)::vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! full level
    real(8),intent(in)::vxh(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(in)::vyh(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(in)::vzh(ADM_gall,ADM_kall,ADM_lall)           ! half level
    real(8),intent(in)::vxh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    real(8),intent(in)::vyh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    real(8),intent(in)::vzh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)  ! half level
    !
    logical,intent(in),optional::output_sclh
    logical ::output_sclh_in=.true.

    integer::k,i,l

    real(8)::tmp(ADM_gall,ADM_kall,ADM_lall)         
    real(8)::tmp_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8)::tmph(ADM_gall,ADM_kall,ADM_lall)         
    real(8)::tmph_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8)::tmp_vx(ADM_gall,ADM_kall,ADM_lall)         
    real(8)::tmp_vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8)::tmp_vy(ADM_gall,ADM_kall,ADM_lall)         
    real(8)::tmp_vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8)::tmp_vz(ADM_gall,ADM_kall,ADM_lall)         
    real(8)::tmp_vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8)::tmp_vxh(ADM_gall,ADM_kall,ADM_lall)         
    real(8)::tmp_vxh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8)::tmp_vyh(ADM_gall,ADM_kall,ADM_lall)         
    real(8)::tmp_vyh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8)::tmp_vzh(ADM_gall,ADM_kall,ADM_lall)         
    real(8)::tmp_vzh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    if (present(output_sclh)) output_sclh_in=output_sclh


    ! initiation because halo is not written
    tmp=0
    tmp_vx=0
    tmp_vy=0
    tmp_vz=0
    tmph=0
    tmp_vxh=0
    tmp_vyh=0
    tmp_vzh=0

    tmp_pl=0
    tmp_vx_pl=0
    tmp_vy_pl=0
    tmp_vz_pl=0
    tmph_pl=0
    tmp_vxh_pl=0
    tmp_vyh_pl=0
    tmp_vzh_pl=0

    scl=0
    sclh=0
    scl_pl=0
    sclh_pl=0

    !-------------------------- for full level -------------------------
    tmp_vx(:,ADM_kmin:ADM_kmax,:)= vx(:,ADM_kmin:ADM_kmax,:)* smg_oprt_gsqrt(:,ADM_kmin:ADM_kmax,:)
    tmp_vy(:,ADM_kmin:ADM_kmax,:)= vy(:,ADM_kmin:ADM_kmax,:)* smg_oprt_gsqrt(:,ADM_kmin:ADM_kmax,:)
    tmp_vz(:,ADM_kmin:ADM_kmax,:)= vz(:,ADM_kmin:ADM_kmax,:)* smg_oprt_gsqrt(:,ADM_kmin:ADM_kmax,:)
    if(ADM_prc_me==ADM_prc_pl) then
       tmp_vx_pl(:,ADM_kmin:ADM_kmax,:)= vx_pl(:,ADM_kmin:ADM_kmax,:)* smg_oprt_gsqrt_pl(:,ADM_kmin:ADM_kmax,:)
       tmp_vy_pl(:,ADM_kmin:ADM_kmax,:)= vy_pl(:,ADM_kmin:ADM_kmax,:)* smg_oprt_gsqrt_pl(:,ADM_kmin:ADM_kmax,:)
       tmp_vz_pl(:,ADM_kmin:ADM_kmax,:)= vz_pl(:,ADM_kmin:ADM_kmax,:)* smg_oprt_gsqrt_pl(:,ADM_kmin:ADM_kmax,:)
    endif
    !
    call OPRT_divergence(&
         tmp,tmp_pl,&
         tmp_vx, tmp_vx_pl,   &
         tmp_vy, tmp_vy_pl,   &
         tmp_vz, tmp_vz_pl    &
         )

!       write(*,*) 'a',smg_oprt_GAM(1144,9,1),tmp(1144,9,1)
!       write(*,*) 'b',tmp_vx(1144,9,1),tmp_vy(1144,9,1),tmp_vz(1144,9,1)
!       write(*,*) 'a',smg_oprt_GAM(1128,9,1),tmp(1128,9,1)
!       write(*,*) 'b',tmp_vx(1128,9,1),tmp_vy(1128,9,1),tmp_vz(1128,9,1)

    do l=1,ADM_lall
    do k=ADM_kmin,ADM_kmax
    do i=1,ADM_gall
!       write(*,*) ADM_prc_me,l,k,i,smg_oprt_GAM(i,k,l),tmp(i,k,l)
    enddo
    enddo
    enddo

    scl(:,ADM_kmin:ADM_kmax,:) = tmp(:,ADM_kmin:ADM_kmax,:) * smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:)
    if(ADM_prc_me==ADM_prc_pl) then
       scl_pl(:,ADM_kmin:ADM_kmax,:) = tmp_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:)
    endif
!    call dbgmx('div',tmp)
!    call dbgmx('div_pl',tmp_pl)
!    call dbgmx('tmpvx',tmp_vx)
!    call dbgmx('tmpvx_pl',tmp_vx_pl)
!    call dbgmx('tmpvy',tmp_vy)
!    call dbgmx('tmpvy_pl',tmp_vy_pl)
!    call dbgmx('tmpvz',tmp_vz)
!    call dbgmx('tmpvz_pl',tmp_vz_pl)
!    call dbgmx('oprt_GAM',smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:))
!    if (ADM_prc_me.eq.1) call dbgmx('oprt_GAM_pl',smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:))
    !-----------
 
    do k=ADM_kmin,ADM_kmax
      tmph(:,k,:) = smg_oprt_gsqrtH(:,k,:) * (&
            VMTR_GZXH(:,k,:)*vxh(:,k,:) &
          + VMTR_GZYH(:,k,:)*vyh(:,k,:) &
          + VMTR_GZZH(:,k,:)*vzh(:,k,:) )
    enddo
    do k=ADM_kmin,ADM_kmax-1
       tmp(:,k,:)=(tmph(:,k+1,:)-tmph(:,k,:))*GRD_rdgz(k)
    enddo
    k=ADM_kmax
    tmp(:,k,:)=(0-tmph(:,k,:))*GRD_rdgz(k) ! flux from top is zero
    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin,ADM_kmax
          tmph_pl(:,k,:) = smg_oprt_gsqrtH_pl(:,k,:) * (&
          + VMTR_GZXH_pl(:,k,:)*vxh_pl(:,k,:) &
          + VMTR_GZYH_pl(:,k,:)*vyh_pl(:,k,:) &
          + VMTR_GZZH_pl(:,k,:)*vzh_pl(:,k,:) )
      enddo
      do k=ADM_kmin,ADM_kmax-1
         tmp_pl(:,k,:)=(tmph_pl(:,k+1,:)-tmph_pl(:,k,:))*GRD_rdgz(k)
      enddo
      k=ADM_kmax
      tmp_pl(:,k,:)=(0-tmph_pl(:,k,:))*GRD_rdgz(k) ! flux from top is zero
    endif
    scl(:,ADM_kmin:ADM_kmax,:)=scl(:,ADM_kmin:ADM_kmax,:) + &
         tmp(:,ADM_kmin:ADM_kmax,:) * smg_oprt_GAM(:,ADM_kmin:ADM_kmax,:)
    if(ADM_prc_me==ADM_prc_pl) then
       scl_pl(:,ADM_kmin:ADM_kmax,:) = scl_pl(:,ADM_kmin:ADM_kmax,:) + &
            tmp_pl(:,ADM_kmin:ADM_kmax,:)*smg_oprt_GAM_pl(:,ADM_kmin:ADM_kmax,:)
    endif

    !-----------

    do k=ADM_kmin,ADM_kmax
      tmph(:,k,:) = VMTR_GAM2H(:,k,:) * (&
            GRD_X(:,ADM_KNONE,:,GRD_XDIR)*vxh(:,k,:) &
          + GRD_X(:,ADM_KNONE,:,GRD_YDIR)*vyh(:,k,:) &
          + GRD_X(:,ADM_KNONE,:,GRD_ZDIR)*vzh(:,k,:) ) / GRD_rscale
    enddo
    do k=ADM_kmin,ADM_kmax-1
       tmp(:,k,:)=(tmph(:,k+1,:)-tmph(:,k,:))*GRD_rdgz(k)
    enddo
    k=ADM_kmax
    tmp(:,k,:)=(0-tmph(:,k,:))*GRD_rdgz(k)  ! flux from top is zero
    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin,ADM_kmax
          tmph_pl(:,k,:) = VMTR_GAM2H_pl(:,k,:) * (&
          + GRD_X_pl(:,ADM_KNONE,:,GRD_XDIR)*vxh_pl(:,k,:) &
          + GRD_X_pl(:,ADM_KNONE,:,GRD_YDIR)*vyh_pl(:,k,:) &
          + GRD_X_pl(:,ADM_KNONE,:,GRD_ZDIR)*vzh_pl(:,k,:) ) / GRD_rscale
      enddo
      do k=ADM_kmin,ADM_kmax-1
         tmp_pl(:,k,:)=(tmph_pl(:,k+1,:)-tmph_pl(:,k,:))*GRD_rdgz(k)
      enddo
      k=ADM_kmax
      tmp_pl(:,k,:)=(0-tmph_pl(:,k,:))*GRD_rdgz(k) ! flux from top is zero
    endif
    scl(:,ADM_kmin:ADM_kmax,:)=scl(:,ADM_kmin:ADM_kmax,:) + tmp(:,ADM_kmin:ADM_kmax,:)
    if(ADM_prc_me==ADM_prc_pl) then
       scl_pl(:,ADM_kmin:ADM_kmax,:) = scl_pl(:,ADM_kmin:ADM_kmax,:) + tmp_pl(:,ADM_kmin:ADM_kmax,:)
    endif

    
    !-------------------------- for half level -------------------------
    if (.not.output_sclh_in) then
       return
    endif
    tmp_vxh(:,ADM_kmin:ADM_kmax,:)= vxh(:,ADM_kmin:ADM_kmax,:)*smg_oprt_gsqrtH(:,ADM_kmin:ADM_kmax,:) ! half level !!
    tmp_vyh(:,ADM_kmin:ADM_kmax,:)= vyh(:,ADM_kmin:ADM_kmax,:)*smg_oprt_gsqrtH(:,ADM_kmin:ADM_kmax,:)
    tmp_vzh(:,ADM_kmin:ADM_kmax,:)= vzh(:,ADM_kmin:ADM_kmax,:)*smg_oprt_gsqrtH(:,ADM_kmin:ADM_kmax,:)
    if(ADM_prc_me==ADM_prc_pl) then
       tmp_vxh_pl(:,ADM_kmin:ADM_kmax,:)= vxh_pl(:,ADM_kmin:ADM_kmax,:)*smg_oprt_gsqrtH_pl(:,ADM_kmin:ADM_kmax,:)
       tmp_vyh_pl(:,ADM_kmin:ADM_kmax,:)= vyh_pl(:,ADM_kmin:ADM_kmax,:)*smg_oprt_gsqrtH_pl(:,ADM_kmin:ADM_kmax,:)
       tmp_vzh_pl(:,ADM_kmin:ADM_kmax,:)= vzh_pl(:,ADM_kmin:ADM_kmax,:)*smg_oprt_gsqrtH_pl(:,ADM_kmin:ADM_kmax,:)
    endif
    !
    call OPRT_divergence(& ! half level !!
         tmph,tmph_pl,&
         tmp_vxh, tmp_vxh_pl,   &
         tmp_vyh, tmp_vyh_pl,   &
         tmp_vzh, tmp_vzh_pl    &
         )
    sclh(:,ADM_kmin:ADM_kmax,:) = tmph(:,ADM_kmin:ADM_kmax,:) * smg_oprt_GAMH(:,ADM_kmin:ADM_kmax,:)
    if(ADM_prc_me==ADM_prc_pl) then
       sclh_pl(:,ADM_kmin:ADM_kmax,:) = tmph_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_GAMH_pl(:,ADM_kmin:ADM_kmax,:)
    endif

    !-----------
 
    do k=ADM_kmin,ADM_kmax
      tmp(:,k,:) = smg_oprt_gsqrt(:,k,:) * (&
            VMTR_GZX(:,k,:)*vx(:,k,:) &
          + VMTR_GZY(:,k,:)*vy(:,k,:) &
          + VMTR_GZZ(:,k,:)*vz(:,k,:) )
    enddo
    do k=ADM_kmin+1,ADM_kmax
       tmph(:,k,:)=(tmp(:,k,:)-tmp(:,k-1,:))*GRD_rdgzH(k)
    enddo
    k=ADM_kmin
    tmph(:,k,:)=(tmp(:,k,:)-0)*GRD_rdgzH(k) ! flux from ground is zero 
    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin,ADM_kmax
          tmp_pl(:,k,:) = smg_oprt_gsqrt_pl(:,k,:) * (&
          + VMTR_GZX_pl(:,k,:)*vx_pl(:,k,:) &
          + VMTR_GZY_pl(:,k,:)*vy_pl(:,k,:) &
          + VMTR_GZZ_pl(:,k,:)*vz_pl(:,k,:) )
      enddo
      do k=ADM_kmin+1,ADM_kmax
         tmph_pl(:,k,:)=(tmp_pl(:,k,:)-tmp_pl(:,k-1,:))*GRD_rdgzh(k)
      enddo
      k=ADM_kmin
      tmph_pl(:,k,:)=(tmp_pl(:,k,:)-0)*GRD_rdgzh(k) ! flux from ground is zero
    endif
    sclh(:,ADM_kmin:ADM_kmax,:)=sclh(:,ADM_kmin:ADM_kmax,:) + tmph(:,ADM_kmin:ADM_kmax,:) * smg_oprt_GAMH(:,ADM_kmin:ADM_kmax,:)
    if(ADM_prc_me==ADM_prc_pl) then
       sclh_pl(:,ADM_kmin:ADM_kmax,:) = sclh_pl(:,ADM_kmin:ADM_kmax,:) + &
            tmph_pl(:,ADM_kmin:ADM_kmax,:) * smg_oprt_GAMH_pl(:,ADM_kmin:ADM_kmax,:)
    endif

    !-----------

    do k=ADM_kmin,ADM_kmax
      tmp(:,k,:) = VMTR_GAM2(:,k,:) * (&
            GRD_X(:,ADM_KNONE,:,GRD_XDIR)*vx(:,k,:) &
          + GRD_X(:,ADM_KNONE,:,GRD_YDIR)*vy(:,k,:) &
          + GRD_X(:,ADM_KNONE,:,GRD_ZDIR)*vz(:,k,:) ) / GRD_rscale
    enddo
    do k=ADM_kmin+1,ADM_kmax
       tmph(:,k,:)=(tmp(:,k,:)-tmp(:,k-1,:))*GRD_rdgzH(k)
    enddo
    k=ADM_kmin
    tmph(:,k,:)=(tmp(:,k,:)-0)*GRD_rdgzH(k) ! flux from ground is zero
    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin,ADM_kmax
          tmp_pl(:,k,:) = VMTR_GAM2_pl(:,k,:) * ( &
          + GRD_X_pl(:,ADM_KNONE,:,GRD_XDIR)*vx_pl(:,k,:) &
          + GRD_X_pl(:,ADM_KNONE,:,GRD_YDIR)*vy_pl(:,k,:) &
          + GRD_X_pl(:,ADM_KNONE,:,GRD_ZDIR)*vz_pl(:,k,:) ) / GRD_rscale
      enddo
      do k=ADM_kmin+1,ADM_kmax
         tmph_pl(:,k,:)=(tmp_pl(:,k,:)-tmp_pl(:,k-1,:))*GRD_rdgzh(k)
      enddo
      k=ADM_kmin
      tmph_pl(:,k,:)=(tmp_pl(:,k,:)-0)*GRD_rdgzh(k) ! flux from ground is zero
    endif
    sclh(:,ADM_kmin:ADM_kmax,:)=sclh(:,ADM_kmin:ADM_kmax,:) + tmph(:,ADM_kmin:ADM_kmax,:)
    if(ADM_prc_me==ADM_prc_pl) then
       sclh_pl(:,ADM_kmin:ADM_kmax,:) = sclh_pl(:,ADM_kmin:ADM_kmax,:) + tmph_pl(:,ADM_kmin:ADM_kmax,:)
    endif

    return
  end subroutine gsqrt_Div3dfh
  !-------------------------------------------------------------------------------
  subroutine dbgmx(cha,var)
    character(*):: cha
    real(8):: var(:,:,:)
!    if (ADM_prc_me.eq.1) then
!      write(*,*) trim(cha),maxval(var),minval(var), maxloc(var),minloc(var)
       write(adm_log_fid,*) trim(cha),maxval(var),minval(var), maxloc(var),minloc(var)
!    endif
  end subroutine dbgmx



end module mod_tb_smg
!-------------------------------------------------------------------------------
