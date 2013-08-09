!-------------------------------------------------------------------------------
!
!+  Vertical Implicit module
!
!-------------------------------------------------------------------------------
module mod_vi
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module is for the vertical implicit scheme of non-hydorostatic
  !       model.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                06-08-11   Add averaged rhog for tracer advection.
  !                11-05-07   Y.Yamada: Implementation of ES tuning cord by NEC.
  !                             Modified line: (20110405 NEC)
  !                                                or (!ftr< vi_small_step.r???)
  !                11-11-28   Y.Yamada: Merge Terai-san timer code
  !                                                    into the original code.
  !                11-12-29   Y.Yamada: Delete ES tuning and merge fjtimer
  !                12-3-9    S.Iga: tuned (phase4-1)
  !
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
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
  public :: vi_small_step

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: vi_path2

  private :: vi_rhow_update_matrix
  private :: vi_rhow

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: vmax_split2 = 6
  integer, private, parameter :: vmax_mean_c = 5

  integer, private, parameter :: I_rhogvx = 1
  integer, private, parameter :: I_rhogvy = 2
  integer, private, parameter :: I_rhogvz = 3
  integer, private, parameter :: I_rhog   = 4
  integer, private, parameter :: I_rhogw  = 5
  integer, private, parameter :: I_rhoge  = 6

  real(8), private, allocatable, save :: Mc   (:,:,:)
  real(8), private, allocatable, save :: Mc_pl(:,:,:)
  real(8), private, allocatable, save :: Ml   (:,:,:)
  real(8), private, allocatable, save :: Ml_pl(:,:,:)
  real(8), private, allocatable, save :: Mu   (:,:,:)
  real(8), private, allocatable, save :: Mu_pl(:,:,:)

  real(8), private, allocatable, save :: A2_o     (:,:,:)
  real(8), private, allocatable, save :: A2_o_pl  (:,:,:)
  real(8), private, allocatable, save :: CooCip   (:,:,:)
  real(8), private, allocatable, save :: CooCip_pl(:,:,:)
  real(8), private, allocatable, save :: CooCim   (:,:,:)
  real(8), private, allocatable, save :: CooCim_pl(:,:,:)
  real(8), private, allocatable, save :: D2       (:,:,:)
  real(8), private, allocatable, save :: D2_pl    (:,:,:)

  real(8), private, save :: alfa

  logical, private, save :: iflag = .true.

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine vi_small_step( &
       rhog,         rhog_pl,         & !--- [INOUT]
       rhogvx,       rhogvx_pl,       & !--- [INOUT]
       rhogvy,       rhogvy_pl,       & !--- [INOUT]
       rhogvz,       rhogvz_pl,       & !--- [INOUT]
       rhogw,        rhogw_pl,        & !--- [INOUT]
       rhoge,        rhoge_pl,        & !--- [INOUT]
       vx,           vx_pl,           & !--- [IN]    : Vx
       vy,           vy_pl,           & !--- [IN]    : Vy
       vz,           vz_pl,           & !--- [IN]    : Vz
       eth,          eth_pl,          & !--- [IN]    : enthalpy
       rhogd,        rhogd_pl,        & !--- [IN]    : perturb dens.( gam2 X G^{1/2} )
       pregd,        pregd_pl,        & !--- [IN]    : perturb pres.( gam2 X G^{1/2} )
       grhog0,       grhog0_pl,       & !--- [IN]
       grhogvx0,     grhogvx0_pl,     & !--- [IN]
       grhogvy0,     grhogvy0_pl,     & !--- [IN]
       grhogvz0,     grhogvz0_pl,     & !--- [IN]
       grhogw0,      grhogw0_pl,      & !--- [IN]
       grhoge0,      grhoge0_pl,      & !--- [IN]
       grhogetot0,   grhogetot0_pl,   & !--- [IN]
       rhog_split,   rhog_split_pl,   & !--- [INOUT] : rhog_split
       rhogvx_split, rhogvx_split_pl, & !--- [INOUT] : rhogvx_split
       rhogvy_split, rhogvy_split_pl, & !--- [INOUT] : rhogvy_split
       rhogvz_split, rhogvz_split_pl, & !--- [INOUT] : rhogvz_split
       rhogw_split,  rhogw_split_pl,  & !--- [INOUT] : rhogw_split
       rhoge_split,  rhoge_split_pl,  & !--- [INOUT] : rhoge_split
       v_mean_c,     v_mean_c_pl,     & !--- [OUT]
       num_of_itr,                    & !--- [IN]
       dt                             ) !--- [IN]
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gmax,    &
       ADM_gmin,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_comm, only: &
       COMM_data_transfer
    use mod_time, only: &
       TIME_SPLIT
    use mod_cnst, only: &
       CNST_EGRAV, &
       CNST_RAIR,  &
       CNST_CV
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac, &
       GRD_cfac, &
       GRD_dfac
    use mod_oprt, only: &
       OPRT_horizontalize_vec
    use mod_vmtr, only: &
       VMTR_GSGAM2H,     &
       VMTR_GSGAM2H_pl,  &
       VMTR_GSGAM2,      &
       VMTR_GSGAM2_pl,   &
       VMTR_RGSGAM2H,    &
       VMTR_RGSGAM2H_pl, &
       VMTR_RGSGAM2,     &
       VMTR_RGSGAM2_pl
    use mod_runconf, only: &
       NON_HYDRO_ALPHA 
    use mod_bndcnd, only: &
       BNDCND_rhovxvyvz
    use mod_numfilter, only: &
       numfilter_divdamp,   &
       numfilter_divdamp_2d
    use mod_src, only: &
       src_advection_convergence, &
       src_flux_convergence,      &
       src_gradient,              &
       src_buoyancy,              &
       I_SRC_default,             &
       I_SRC_horizontal
    implicit none

    real(8), intent(inout) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   ) ! density ( gam2 X G^{1/2} )
    real(8), intent(inout) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx  ( gam2 X G^{1/2} )
    real(8), intent(inout) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy  ( gam2 X G^{1/2} )
    real(8), intent(inout) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vz  ( gam2 X G^{1/2} )
    real(8), intent(inout) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w   ( gam2 X G^{1/2} )
    real(8), intent(inout) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*ein ( gam2 X G^{1/2} )
    real(8), intent(inout) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in)    :: vx      (ADM_gall,   ADM_kall,ADM_lall   ) ! diagnostic variables
    real(8), intent(in)    :: vx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: eth     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: eth_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: pregd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: pregd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)    :: grhog0       (ADM_gall,   ADM_kall,ADM_lall   ) ! source term (large step)
    real(8), intent(in)    :: grhog0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogvx0     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogvx0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogvy0     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogvy0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogvz0     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogvz0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogw0      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogw0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhoge0      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhoge0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogetot0   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogetot0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(inout) :: rhog_split     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhog_split_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvx_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhogvx_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvy_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhogvy_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvz_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhogvz_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogw_split    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhogw_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhoge_split    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhoge_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(out)   :: v_mean_c   (ADM_gall,   ADM_kall,ADM_lall   ,vmax_mean_c) ! mean_flux for tracer advection
    real(8), intent(out)   :: v_mean_c_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax_mean_c)

    integer, intent(in)    :: num_of_itr
    real(8), intent(in)    :: dt

    !--- tendency term (large step + small step)
    real(8) :: grhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: grhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendency term 2
    real(8) :: drhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- divergence damping
    real(8) :: gdx       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gdx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gdy       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gdy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gdz       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gdz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gdvz      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gdvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: gdx_2d    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gdx_2d_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gdy_2d    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gdy_2d_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gdz_2d    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gdz_2d_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- pressure gradient force
    real(8) :: gpx       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gpx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gpy       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gpy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gpz       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gpz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gpvz      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gpvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- buoyancy force
    real(8) :: gbvz      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gbvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- pressure work
    real(8) :: drhoge_pw   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhoge_pw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: drhog_h     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhog_h_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gz_tilde    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gz_tilde_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gpzw        (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: gpzw_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogez     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhogez_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: eth_h       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: eth_h_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: pregd_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: pregd_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- for communication
    real(8) :: v_split2   (ADM_gall,   ADM_kall,ADM_lall   ,vmax_split2)
    real(8) :: v_split2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax_split2)
    logical :: comm_flag1(vmax_split2)
    logical :: comm_flag2(vmax_split2)

    real(8) :: rweight_itr

    integer :: ij, k, l, ns

    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    !--- comm_flag1 ( rhogvx_split2, rhogvy_split2, rhogvz_split2 )
    comm_flag1(:) = .false.
    comm_flag1(I_rhogvx) = .true.
    comm_flag1(I_rhogvy) = .true.
    comm_flag1(I_rhogvz) = .true.

    !--- comm_flag2 ( rhog_split2, rhogw_split2, rhoge_split2 )
    comm_flag2(:) = .false.
    comm_flag2(I_rhog)  = .true.
    comm_flag2(I_rhogw) = .true.
    comm_flag2(I_rhoge) = .true.

    !---------------------------------------------------------------------------
    ! vi_path0
    !---------------------------------------------------------------------------

    !---< Calculation of drhog
    call src_flux_convergence( rhogvx, rhogvx_pl, & !--- [IN]
                               rhogvy, rhogvy_pl, & !--- [IN]
                               rhogvz, rhogvz_pl, & !--- [IN]
                               rhogw,  rhogw_pl,  & !--- [IN]
                               drhog,  drhog_pl,  & !--- [OUT]
                               I_SRC_default      ) !--- [IN]

    !---< Calculation of drhogvx, drhogvy, drhogvz, drhogw

    !--- divergence damping
    call numfilter_divdamp( rhogvx, rhogvx_pl, & !--- [IN]
                            rhogvy, rhogvy_pl, & !--- [IN]
                            rhogvz, rhogvz_pl, & !--- [IN]
                            rhogw,  rhogw_pl,  & !--- [IN]
                            gdx,    gdx_pl,    & !--- [OUT]
                            gdy,    gdy_pl,    & !--- [OUT]
                            gdz,    gdz_pl,    & !--- [OUT]
                            gdvz,   gdvz_pl    ) !--- [OUT]

    call numfilter_divdamp_2d( rhogvx, rhogvx_pl, & !--- [IN]
                               rhogvy, rhogvy_pl, & !--- [IN]
                               rhogvz, rhogvz_pl, & !--- [IN]
                               gdx_2d, gdx_2d_pl, & !--- [OUT]
                               gdy_2d, gdy_2d_pl, & !--- [OUT]
                               gdz_2d, gdz_2d_pl  ) !--- [OUT]

    !--- pressure force
    call src_gradient( pregd, pregd_pl, & !--- [IN]
                       gpx,   gpx_pl,   & !--- [OUT]
                       gpy,   gpy_pl,   & !--- [OUT]
                       gpz,   gpz_pl,   & !--- [OUT]
                       gpvz,  gpvz_pl,  & !--- [OUT]
                       I_SRC_default    ) !--- [IN]

    !--- buoyancy force
    call src_buoyancy( rhogd, rhogd_pl, & !--- [IN]
                       gbvz,  gbvz_pl   ) !--- [OUT]

    !---< Calculation of drhoge

    !--- advection convergence for eth
    call src_advection_convergence( rhogvx, rhogvx_pl, & !--- [IN]
                                    rhogvy, rhogvy_pl, & !--- [IN]
                                    rhogvz, rhogvz_pl, & !--- [IN]
                                    rhogw,  rhogw_pl,  & !--- [IN]
                                    eth,    eth_pl,    & !--- [IN]
                                    drhoge, drhoge_pl, & !--- [OUT]
                                    I_SRC_default      ) !--- [IN]

    !--- pressure work (horizontal)
    drhoge_pw(:,:,:) = vx(:,:,:) * gpx(:,:,:) &
                     + vy(:,:,:) * gpy(:,:,:) &
                     + vz(:,:,:) * gpz(:,:,:)

    if ( ADM_prc_me == ADM_prc_pl ) then
       drhoge_pw_pl(:,:,:) = vx_pl(:,:,:) * gpx_pl(:,:,:) &
                           + vy_pl(:,:,:) * gpy_pl(:,:,:) &
                           + vz_pl(:,:,:) * gpz_pl(:,:,:)
    endif


    !--- Caluculation of (-rhogw * g_tilde) and interpolation to half level
    do l = 1, ADM_lall
       !--- Calculation of gz_tilde
       do k  = ADM_kmin, ADM_kmax+1
       do ij = 1, ADM_gall
          drhog_h(ij,k,l) = 0.5D0 * ( GRD_afac(k) * VMTR_RGSGAM2(ij,k  ,l) * rhog(ij,k,  l) &
                                    + GRD_bfac(k) * VMTR_RGSGAM2(ij,k-1,l) * rhog(ij,k-1,l) &
                                    ) * VMTR_GSGAM2H(ij,k,l)
       enddo
       enddo
       do ij = 1, ADM_gall
          drhog_h(ij,ADM_kmin-1,l) = drhog_h(ij,ADM_kmin,l)
       enddo

       do k = 1, ADM_kall
          do ij = 1, ADM_gall
             gz_tilde(ij,k,l) = CNST_EGRAV - ( gpvz(ij,k,l)-gbvz(ij,k,l) ) / drhog_h(ij,k,l)
             gpzw    (ij,k,l) = -gz_tilde(ij,k,l) * rhogw(ij,k,l)
          enddo
       enddo

       do k  = ADM_kmin, ADM_kmax
       do ij = 1, ADM_gall
          drhogez(ij,k,l) = 0.5D0 * ( GRD_cfac(k) * VMTR_RGSGAM2H(ij,k+1,l) * gpzw(ij,k+1,l) &
                                    + GRD_dfac(k) * VMTR_RGSGAM2H(ij,k,  l) * gpzw(ij,k,  l) &
                                    ) * VMTR_GSGAM2(ij,k,l)
       enddo
       enddo
       do ij = 1, ADM_gall
          drhogez(ij,ADM_kmin-1,l) = 0.D0
          drhogez(ij,ADM_kmax+1,l) = 0.D0
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k  = ADM_kmin, ADM_kmax+1
          do ij = 1, ADM_gall_pl
             drhog_h_pl(ij,k,l) = 0.5D0 * ( GRD_afac(k) * VMTR_RGSGAM2_pl(ij,k  ,l)*rhog_pl(ij,k,  l) &
                                          + GRD_bfac(k) * VMTR_RGSGAM2_pl(ij,k-1,l)*rhog_pl(ij,k-1,l) &
                                          ) * VMTR_GSGAM2H_pl(ij,k,l)
          enddo
          enddo
          do ij = 1, ADM_gall_pl
             drhog_h_pl(ij,ADM_kmin-1,l) = drhog_h_pl(ij,ADM_kmin,l)
          enddo

          do k = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                gz_tilde_pl(ij,k,l) = CNST_EGRAV - ( gpvz_pl(ij,k,l) - gbvz_pl(ij,k,l) ) / drhog_h_pl(ij,k,l)
                gpzw_pl    (ij,k,l) = -gz_tilde_pl(ij,k,l) * rhogw_pl(ij,k,l)
             enddo
          enddo

          do k  = ADM_kmin, ADM_kmax
          do ij = 1, ADM_gall_pl
             drhogez_pl(ij,k,l) = 0.5D0 * ( GRD_cfac(k) * gpzw_pl(ij,k+1,l) * VMTR_RGSGAM2H_pl(ij,k+1,l) &
                                          + GRD_dfac(k) * gpzw_pl(ij,k,  l) * VMTR_RGSGAM2H_pl(ij,k,  l) &
                                          ) * VMTR_GSGAM2_pl(ij,k,l)
          enddo
          enddo
          do ij = 1, ADM_gall_pl
             drhogez_pl(ij,ADM_kmin-1,l) = 0.D0
             drhogez_pl(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo
    endif

    !--- sum of tendencies ( large step + pres-grad + div-damp + div-damp_2d + buoyancy )
    do l  = 1, ADM_lall
       do k  = 1, ADM_kall
          do ij = 1, ADM_gall
             grhog  (ij,k,l) = grhog0  (ij,k,l) + drhog(ij,k,l)

             grhogvx(ij,k,l) = grhogvx0(ij,k,l) - gpx (ij,k,l) + gdx(ij,k,l) + gdx_2d(ij,k,l)
             grhogvy(ij,k,l) = grhogvy0(ij,k,l) - gpy (ij,k,l) + gdy(ij,k,l) + gdy_2d(ij,k,l)
             grhogvz(ij,k,l) = grhogvz0(ij,k,l) - gpz (ij,k,l) + gdz(ij,k,l) + gdz_2d(ij,k,l)

             grhogw (ij,k,l) = grhogw0 (ij,k,l) - gpvz(ij,k,l) + gdvz(ij,k,l) * NON_HYDRO_ALPHA &
                                                + gbvz(ij,k,l)

             grhoge (ij,k,l) = grhoge0 (ij,k,l) + drhoge   (ij,k,l) &
                                                + drhoge_pw(ij,k,l) &
                                                + drhogez  (ij,k,l)
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l  = 1, ADM_lall_pl
          do k  = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                grhog_pl  (ij,k,l) = grhog0_pl  (ij,k,l) + drhog_pl(ij,k,l)

                grhogvx_pl(ij,k,l) = grhogvx0_pl(ij,k,l) - gpx_pl (ij,k,l) + gdx_pl (ij,k,l)
                grhogvy_pl(ij,k,l) = grhogvy0_pl(ij,k,l) - gpy_pl (ij,k,l) + gdy_pl (ij,k,l)
                grhogvz_pl(ij,k,l) = grhogvz0_pl(ij,k,l) - gpz_pl (ij,k,l) + gdz_pl (ij,k,l)

                grhogw_pl (ij,k,l) = grhogw0_pl (ij,k,l) - gpvz_pl(ij,k,l) + gdvz_pl(ij,k,l) * NON_HYDRO_ALPHA &
                                                         + gbvz_pl(ij,k,l)

                grhoge_pl (ij,k,l) = grhoge0_pl (ij,k,l) + drhoge_pl   (ij,k,l) &
                                                         + drhoge_pw_pl(ij,k,l) &
                                                         + drhogez_pl  (ij,k,l)
             enddo
          enddo
       enddo
    endif

    !---------------------------------------------------------------------------
    ! END vi_path0
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
          do ij = 1, ADM_gall
             eth_h(ij,k,l) = 0.5D0 * ( GRD_afac(k) * eth(ij,k,  l) &
                                     + GRD_bfac(k) * eth(ij,k-1,l) )
          enddo
       enddo
       do ij = 1, ADM_gall
          eth_h(ij,ADM_kmin-1,l) = eth_h(ij,ADM_kmin,l)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
             do ij = 1, ADM_gall_pl
                eth_h_pl(ij,k,l) = 0.5D0 * ( GRD_afac(k) * eth_pl(ij,k,  l) &
                                           + GRD_bfac(k) * eth_pl(ij,k-1,l) )
             enddo
          enddo
          do ij = 1, ADM_gall_pl
             eth_h_pl(ij,ADM_kmin-1,l) = eth_h_pl(ij,ADM_kmin,l)
          enddo
       enddo
    endif

    ! update working matrix in mod_rhow
    call vi_rhow_update_matrix( eth_h,    eth_h_pl,    & !--- [IN] : enthalpy at the h-lev
                                gz_tilde, gz_tilde_pl, & !--- [IN] : effective gravitation at the h-lev
                                dt                     ) !--- [IN] : delta t

    !--- initialization of mean mass flux
    rweight_itr = 1.D0 / real(num_of_itr,kind=8)

    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do ij = 1, ADM_gall
             v_mean_c(ij,k,l,I_rhogvx) = rhogvx(ij,k,l)
             v_mean_c(ij,k,l,I_rhogvy) = rhogvy(ij,k,l)
             v_mean_c(ij,k,l,I_rhogvz) = rhogvz(ij,k,l)
             v_mean_c(ij,k,l,I_rhogw)  = rhogw (ij,k,l)
             v_mean_c(ij,k,l,I_rhog)   = rhog  (ij,k,l)
          enddo
       enddo
    enddo

    do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          do ij = 1, ADM_gall_pl
             v_mean_c_pl(ij,k,l,I_rhogvx) = rhogvx_pl(ij,k,l)
             v_mean_c_pl(ij,k,l,I_rhogvy) = rhogvy_pl(ij,k,l)
             v_mean_c_pl(ij,k,l,I_rhogvz) = rhogvz_pl(ij,k,l)
             v_mean_c_pl(ij,k,l,I_rhogw)  = rhogw_pl (ij,k,l)
             v_mean_c_pl(ij,k,l,I_rhog)   = rhog_pl  (ij,k,l)
          enddo
       enddo
    enddo

    !---------------------------------------------------------------------------
    !
    !> Start small step iteration
    !
    !---------------------------------------------------------------------------
    do ns = 1, num_of_itr

       !---< calculation of pregd(*) from rhog(*) & rhoge(*)
       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do ij = 1, ADM_gall
                pregd_split(ij,k,l) = rhoge_split(ij,k,l) * ( CNST_RAIR / CNST_CV )
             enddo
          enddo
          do ij = 1, ADM_gall
             pregd_split(ij,ADM_kmin-1,l) = pregd_split(ij,ADM_kmin,l)
             pregd_split(ij,ADM_kmax+1,l) = pregd_split(ij,ADM_kmax,l)
             rhoge_split(ij,ADM_kmin-1,l) = rhoge_split(ij,ADM_kmin,l)
             rhoge_split(ij,ADM_kmax+1,l) = rhoge_split(ij,ADM_kmax,l)
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do ij = 1, ADM_gall_pl
                   pregd_split_pl(ij,k,l) = rhoge_split_pl(ij,k,l) * ( CNST_RAIR / CNST_CV )
                enddo
             enddo
             do ij = 1, ADM_gall_pl
                pregd_split_pl(ij,ADM_kmin-1,l) = pregd_split_pl(ij,ADM_kmin,l)
                pregd_split_pl(ij,ADM_kmax+1,l) = pregd_split_pl(ij,ADM_kmax,l)
                rhoge_split_pl(ij,ADM_kmin-1,l) = rhoge_split_pl(ij,ADM_kmin,l)
                rhoge_split_pl(ij,ADM_kmax+1,l) = rhoge_split_pl(ij,ADM_kmax,l)
             enddo
          enddo
       endif

       if ( TIME_SPLIT ) then
          !---------------------------------------------------------------------
          ! vi_path1
          !---------------------------------------------------------------------

          !---< Calculation of grhogvx, grhogvy, grhogvz, grhogw

          !--- divergence damping
          call numfilter_divdamp( rhogvx_split, rhogvx_split_pl, & !--- [IN]
                                  rhogvy_split, rhogvy_split_pl, & !--- [IN]
                                  rhogvz_split, rhogvz_split_pl, & !--- [IN]
                                  rhogw_split,  rhogw_split_pl,  & !--- [IN]
                                  gdx,          gdx_pl,          & !--- [OUT]
                                  gdy,          gdy_pl,          & !--- [OUT]
                                  gdz,          gdz_pl,          & !--- [OUT]
                                  gdvz,         gdvz_pl          ) !--- [OUT]

          !--- 2d divergence damping
          call numfilter_divdamp_2d( rhogvx_split, rhogvx_split_pl, & !--- [IN]
                                     rhogvy_split, rhogvy_split_pl, & !--- [IN]
                                     rhogvz_split, rhogvz_split_pl, & !--- [IN]
                                     gdx_2d,       gdx_2d_pl,       & !--- [OUT]
                                     gdy_2d,       gdy_2d_pl,       & !--- [OUT]
                                     gdz_2d,       gdz_2d_pl        ) !--- [OUT]

          !--- pressure force
          !--- gpvz=0.D0 becaude of f_type='HORIZONTAL'.
          call src_gradient( pregd_split, pregd_split_pl, & !--- [IN]
                             gpx,         gpx_pl,         & !--- [OUT]
                             gpy,         gpy_pl,         & !--- [OUT]
                             gpz,         gpz_pl,         & !--- [OUT]
                             gpvz,        gpvz_pl,        & !--- [OUT]
                             I_SRC_horizontal             ) !--- [IN]

          !--- buoyancy force
          !--- not calculated, because this term is implicit.

          !---------------------------------------------------------------------
          ! END vi_path1
          !---------------------------------------------------------------------

          !--- sum of tendency ( large step + pres-grad + div-damp + div-damp_2d )
          do l = 1, ADM_lall
             do k = 1, ADM_kall
                do ij = 1, ADM_gall
                   drhogvx(ij,k,l) = grhogvx(ij,k,l) - gpx(ij,k,l) + gdx(ij,k,l) + gdx_2d(ij,k,l)
                   drhogvy(ij,k,l) = grhogvy(ij,k,l) - gpy(ij,k,l) + gdy(ij,k,l) + gdy_2d(ij,k,l)
                   drhogvz(ij,k,l) = grhogvz(ij,k,l) - gpz(ij,k,l) + gdz(ij,k,l) + gdz_2d(ij,k,l)
                   drhogw(ij,k,l)  = grhogw(ij,k,l)                + gdvz(ij,k,l) * NON_HYDRO_ALPHA
                enddo
             enddo
          enddo

          if ( ADM_prc_me == ADM_prc_pl ) then       
             do l = 1, ADM_lall_pl
                do k = 1, ADM_kall
                   do ij = 1, ADM_gall_pl
                      drhogvx_pl(ij,k,l) = grhogvx_pl(ij,k,l) - gpx_pl(ij,k,l) + gdx_pl(ij,k,l) + gdx_2d_pl(ij,k,l)
                      drhogvy_pl(ij,k,l) = grhogvy_pl(ij,k,l) - gpy_pl(ij,k,l) + gdy_pl(ij,k,l) + gdy_2d_pl(ij,k,l)
                      drhogvz_pl(ij,k,l) = grhogvz_pl(ij,k,l) - gpz_pl(ij,k,l) + gdz_pl(ij,k,l) + gdz_2d_pl(ij,k,l)
                      drhogw_pl(ij,k,l)  = grhogw_pl(ij,k,l)                   + gdvz_pl(ij,k,l) * NON_HYDRO_ALPHA
                   enddo
                enddo
             enddo
          endif

       else !--- NO-SPLITING
          !---------------------------------------------------------------------
          ! vi_path1 is skipped
          !---------------------------------------------------------------------

          !------ sum of tendency ( large step )
          do l = 1, ADM_lall
             do k = 1, ADM_kall
                do ij = 1, ADM_gall
                   drhogvx(ij,k,l) = grhogvx(ij,k,l)
                   drhogvy(ij,k,l) = grhogvy(ij,k,l)
                   drhogvz(ij,k,l) = grhogvz(ij,k,l)
                   drhogw(ij,k,l)  = grhogw(ij,k,l)
                enddo
             enddo
          enddo

          if ( ADM_prc_me == ADM_prc_pl ) then
             do l = 1, ADM_lall_pl
                do k = 1, ADM_kall
                   do ij = 1, ADM_gall_pl
                      drhogvx_pl(ij,k,l) = grhogvx_pl(ij,k,l)
                      drhogvy_pl(ij,k,l) = grhogvy_pl(ij,k,l)
                      drhogvz_pl(ij,k,l) = grhogvz_pl(ij,k,l)
                      drhogw_pl(ij,k,l)  = grhogw_pl(ij,k,l)
                   enddo
                enddo
             enddo
          endif

       endif ! Split/Non-split

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do ij = 1, ADM_gall
                v_split2(ij,k,l,I_rhogvx) = rhogvx_split(ij,k,l) + dt * drhogvx(ij,k,l)
                v_split2(ij,k,l,I_rhogvy) = rhogvy_split(ij,k,l) + dt * drhogvy(ij,k,l)
                v_split2(ij,k,l,I_rhogvz) = rhogvz_split(ij,k,l) + dt * drhogvz(ij,k,l)
             enddo
          enddo

          call BNDCND_rhovxvyvz( ADM_gall,                 & !--- [IN]
                                 rhog(:,:,l),              & !--- [IN]
                                 v_split2(:,:,l,I_rhogvx), & !--- [INOUT]
                                 v_split2(:,:,l,I_rhogvy), & !--- [INOUT]
                                 v_split2(:,:,l,I_rhogvz)  ) !--- [INOUT]
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do ij = 1, ADM_gall_pl
                   v_split2_pl(ij,k,l,I_rhogvx) = rhogvx_split_pl(ij,k,l) + dt * drhogvx_pl(ij,k,l)
                   v_split2_pl(ij,k,l,I_rhogvy) = rhogvy_split_pl(ij,k,l) + dt * drhogvy_pl(ij,k,l)
                   v_split2_pl(ij,k,l,I_rhogvz) = rhogvz_split_pl(ij,k,l) + dt * drhogvz_pl(ij,k,l)
                enddo
             enddo

             call BNDCND_rhovxvyvz( ADM_gall_pl,                 & !--- [IN]
                                    rhog_pl(:,:,l),              & !--- [IN]
                                    v_split2_pl(:,:,l,I_rhogvx), & !--- [INOUT]
                                    v_split2_pl(:,:,l,I_rhogvy), & !--- [INOUT]
                                    v_split2_pl(:,:,l,I_rhogvz)  ) !--- [INOUT]
          enddo
       endif

       !--- communication of horizontal momentum
       call COMM_data_transfer(v_split2,v_split2_pl,trn=comm_flag1)

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             v_split2(suf(ADM_gall_1d,1),k,l,I_rhogvx) = v_split2(suf(ADM_gmax+1,ADM_gmin),k,l,I_rhogvx)
             v_split2(suf(1,ADM_gall_1d),k,l,I_rhogvx) = v_split2(suf(ADM_gmin,ADM_gmax+1),k,l,I_rhogvx)
             v_split2(suf(ADM_gall_1d,1),k,l,I_rhogvy) = v_split2(suf(ADM_gmax+1,ADM_gmin),k,l,I_rhogvy)
             v_split2(suf(1,ADM_gall_1d),k,l,I_rhogvy) = v_split2(suf(ADM_gmin,ADM_gmax+1),k,l,I_rhogvy)
             v_split2(suf(ADM_gall_1d,1),k,l,I_rhogvz) = v_split2(suf(ADM_gmax+1,ADM_gmin),k,l,I_rhogvz)
             v_split2(suf(1,ADM_gall_1d),k,l,I_rhogvz) = v_split2(suf(ADM_gmin,ADM_gmax+1),k,l,I_rhogvz)
          enddo
       enddo

       !--- 2nd small step : vertical implicit : rhog, rhogw, pregd: next step
       call vi_path2( v_split2(:,:,:,I_rhog),   v_split2_pl(:,:,:,I_rhog),   & !--- [OUT]
                      v_split2(:,:,:,I_rhogvx), v_split2_pl(:,:,:,I_rhogvx), & !--- [IN]
                      v_split2(:,:,:,I_rhogvy), v_split2_pl(:,:,:,I_rhogvy), & !--- [IN]
                      v_split2(:,:,:,I_rhogvz), v_split2_pl(:,:,:,I_rhogvz), & !--- [IN]
                      v_split2(:,:,:,I_rhogw),  v_split2_pl(:,:,:,I_rhogw),  & !--- [OUT]
                      v_split2(:,:,:,I_rhoge),  v_split2_pl(:,:,:,I_rhoge),  & !--- [OUT]
                      rhog_split,   rhog_split_pl,                           & !--- [IN]
                      rhogvx_split, rhogvx_split_pl,                         & !--- [IN]
                      rhogvy_split, rhogvy_split_pl,                         & !--- [IN]
                      rhogvz_split, rhogvz_split_pl,                         & !--- [IN]
                      rhogw_split,  rhogw_split_pl,                          & !--- [IN]
                      rhoge_split,  rhoge_split_pl,                          & !--- [IN]
                      pregd_split,  pregd_split_pl,                          & !--- [IN]
                      rhog,         rhog_pl,                                 & !--- [IN]
                      rhogvx,       rhogvx_pl,                               & !--- [IN]
                      rhogvy,       rhogvy_pl,                               & !--- [IN]
                      rhogvz,       rhogvz_pl,                               & !--- [IN]
                      rhogw,        rhogw_pl,                                & !--- [IN]
                      eth,          eth_pl,                                  & !--- [IN]
                      grhog,        grhog_pl,                                & !--- [IN]
                      drhogw,       drhogw_pl,                               & !--- [IN]
                      grhoge,       grhoge_pl,                               & !--- [IN]
                      grhogetot0,   grhogetot0_pl,                           & !--- [IN]
                      dt                                                     ) !--- [IN]

       !--- communication of rhog, rhogw, and rhoge
       call COMM_data_transfer(v_split2,v_split2_pl,trn=comm_flag2)

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             v_split2(suf(ADM_gall_1d,1),k,l,I_rhog ) = v_split2(suf(ADM_gmax+1,ADM_gmin),k,l,I_rhog )
             v_split2(suf(1,ADM_gall_1d),k,l,I_rhog ) = v_split2(suf(ADM_gmin,ADM_gmax+1),k,l,I_rhog )
             v_split2(suf(ADM_gall_1d,1),k,l,I_rhogw) = v_split2(suf(ADM_gmax+1,ADM_gmin),k,l,I_rhogw)
             v_split2(suf(1,ADM_gall_1d),k,l,I_rhogw) = v_split2(suf(ADM_gmin,ADM_gmax+1),k,l,I_rhogw)
             v_split2(suf(ADM_gall_1d,1),k,l,I_rhoge) = v_split2(suf(ADM_gmax+1,ADM_gmin),k,l,I_rhoge)
             v_split2(suf(1,ADM_gall_1d),k,l,I_rhoge) = v_split2(suf(ADM_gmin,ADM_gmax+1),k,l,I_rhoge)
          enddo
       enddo

       !--- update for next step
       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do ij = 1, ADM_gall
                rhog_split(ij,k,l)   = v_split2(ij,k,l,I_rhog)
                rhogvx_split(ij,k,l) = v_split2(ij,k,l,I_rhogvx)
                rhogvy_split(ij,k,l) = v_split2(ij,k,l,I_rhogvy)
                rhogvz_split(ij,k,l) = v_split2(ij,k,l,I_rhogvz)
                rhogw_split(ij,k,l)  = v_split2(ij,k,l,I_rhogw)
                rhoge_split(ij,k,l)  = v_split2(ij,k,l,I_rhoge)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do ij = 1, ADM_gall_pl
                   rhog_split_pl(ij,k,l)   = v_split2_pl(ij,k,l,I_rhog)
                   rhogvx_split_pl(ij,k,l) = v_split2_pl(ij,k,l,I_rhogvx)
                   rhogvy_split_pl(ij,k,l) = v_split2_pl(ij,k,l,I_rhogvy)
                   rhogvz_split_pl(ij,k,l) = v_split2_pl(ij,k,l,I_rhogvz)
                   rhogw_split_pl(ij,k,l)  = v_split2_pl(ij,k,l,I_rhogw)
                   rhoge_split_pl(ij,k,l)  = v_split2_pl(ij,k,l,I_rhoge)
                enddo
             enddo
          enddo
       endif

       !--- calculation of mean mass flux ( for tracers )
       do l = 1, ADM_lall
          do k  = 1, ADM_kall
             do ij = 1, ADM_gall
                v_mean_c(ij,k,l,I_rhog)   = v_mean_c(ij,k,l,I_rhog)   + rhog_split  (ij,k,l) * rweight_itr
                v_mean_c(ij,k,l,I_rhogvx) = v_mean_c(ij,k,l,I_rhogvx) + rhogvx_split(ij,k,l) * rweight_itr
                v_mean_c(ij,k,l,I_rhogvy) = v_mean_c(ij,k,l,I_rhogvy) + rhogvy_split(ij,k,l) * rweight_itr
                v_mean_c(ij,k,l,I_rhogvz) = v_mean_c(ij,k,l,I_rhogvz) + rhogvz_split(ij,k,l) * rweight_itr
                v_mean_c(ij,k,l,I_rhogw)  = v_mean_c(ij,k,l,I_rhogw)  + rhogw_split (ij,k,l) * rweight_itr
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k  = 1, ADM_kall
                do ij = 1, ADM_gall_pl
                   v_mean_c_pl(ij,k,l,I_rhog)   = v_mean_c_pl(ij,k,l,I_rhog)   + rhog_split_pl  (ij,k,l) * rweight_itr
                   v_mean_c_pl(ij,k,l,I_rhogvx) = v_mean_c_pl(ij,k,l,I_rhogvx) + rhogvx_split_pl(ij,k,l) * rweight_itr
                   v_mean_c_pl(ij,k,l,I_rhogvy) = v_mean_c_pl(ij,k,l,I_rhogvy) + rhogvy_split_pl(ij,k,l) * rweight_itr
                   v_mean_c_pl(ij,k,l,I_rhogvz) = v_mean_c_pl(ij,k,l,I_rhogvz) + rhogvz_split_pl(ij,k,l) * rweight_itr
                   v_mean_c_pl(ij,k,l,I_rhogw)  = v_mean_c_pl(ij,k,l,I_rhogw)  + rhogw_split_pl (ij,k,l) * rweight_itr
                enddo
             enddo
          enddo
       endif

    enddo  ! small step end

    !--- update prognostic variables
    do l = 1, ADM_lall
       do k  = 1, ADM_kall
       do ij = 1, ADM_gall
             rhog  (ij,k,l) = rhog  (ij,k,l) + rhog_split  (ij,k,l)
             rhogvx(ij,k,l) = rhogvx(ij,k,l) + rhogvx_split(ij,k,l)
             rhogvy(ij,k,l) = rhogvy(ij,k,l) + rhogvy_split(ij,k,l)
             rhogvz(ij,k,l) = rhogvz(ij,k,l) + rhogvz_split(ij,k,l)
             rhogw (ij,k,l) = rhogw (ij,k,l) + rhogw_split (ij,k,l)
             rhoge (ij,k,l) = rhoge (ij,k,l) + rhoge_split (ij,k,l)
       enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k  = 1, ADM_kall
          do ij = 1, ADM_gall_pl
             rhog_pl  (ij,k,l) = rhog_pl  (ij,k,l) + rhog_split_pl  (ij,k,l)
             rhogvx_pl(ij,k,l) = rhogvx_pl(ij,k,l) + rhogvx_split_pl(ij,k,l)
             rhogvy_pl(ij,k,l) = rhogvy_pl(ij,k,l) + rhogvy_split_pl(ij,k,l)
             rhogvz_pl(ij,k,l) = rhogvz_pl(ij,k,l) + rhogvz_split_pl(ij,k,l)
             rhogw_pl (ij,k,l) = rhogw_pl (ij,k,l) + rhogw_split_pl (ij,k,l)
             rhoge_pl (ij,k,l) = rhoge_pl (ij,k,l) + rhoge_split_pl (ij,k,l)
          enddo
          enddo
       enddo
    endif

    call OPRT_horizontalize_vec( rhogvx, rhogvx_pl, & !--- [INOUT]
                                 rhogvy, rhogvy_pl, & !--- [INOUT]
                                 rhogvz, rhogvz_pl  ) !--- [INOUT]

    !--- communication of mean velocity
    call COMM_data_transfer( v_mean_c, v_mean_c_pl )

    return
  end subroutine vi_small_step

  !-----------------------------------------------------------------------------
  subroutine vi_path2( &
       rhog_split2,   rhog_split2_pl,   & !--- [OUT]
       rhogvx_split2, rhogvx_split2_pl, & !--- [IN]
       rhogvy_split2, rhogvy_split2_pl, & !--- [IN]
       rhogvz_split2, rhogvz_split2_pl, & !--- [IN]
       rhogw_split2,  rhogw_split2_pl,  & !--- [OUT]
       rhoge_split2,  rhoge_split2_pl,  & !--- [OUT]
       rhog_split,    rhog_split_pl,    & !--- [IN]
       rhogvx_split,  rhogvx_split_pl,  & !--- [IN]
       rhogvy_split,  rhogvy_split_pl,  & !--- [IN]
       rhogvz_split,  rhogvz_split_pl,  & !--- [IN]
       rhogw_split,   rhogw_split_pl,   & !--- [IN]
       rhoge_split,   rhoge_split_pl,   & !--- [IN]
       pregd_split,   pregd_split_pl,   & !--- [IN]
       rhog,          rhog_pl,          & !--- [IN]
       rhogvx,        rhogvx_pl,        & !--- [IN]
       rhogvy,        rhogvy_pl,        & !--- [IN]
       rhogvz,        rhogvz_pl,        & !--- [IN]
       rhogw,         rhogw_pl,         & !--- [IN]
       eth,           eth_pl,           & !--- [IN]
       grhog,         grhog_pl,         & !--- [IN]
       grhogw,        grhogw_pl,        & !--- [IN]
       grhoge,        grhoge_pl,        & !--- [IN]
       grhogetot,     grhogetot_pl,     & !--- [IN]
       dt                               ) !--- [IN]
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_cnst, only: &
       CNST_RAIR, &
       CNST_CV
    use mod_time, only: &
       TIME_SPLIT
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    use mod_src, only: &
       src_flux_convergence,      &
       src_advection_convergence, &
       I_SRC_horizontal,          &
       I_SRC_default
    use mod_cnvvar, only: &
       cnvvar_rhokin, &
       cnvvar_rhokin_ijkl
    use mod_bsstate, only: &
       phi, phi_pl
    use mod_bndcnd, only: &
       BNDCND_rhow
    implicit none

    real(8), intent(out) :: rhog_split2     (ADM_gall,   ADM_kall,ADM_lall   ) ! prognostic vars (split, at n step)
    real(8), intent(out) :: rhog_split2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx_split2   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_split2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy_split2   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_split2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz_split2   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_split2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogw_split2    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogw_split2_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhoge_split2    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhoge_split2_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: rhog_split      (ADM_gall,   ADM_kall,ADM_lall   ) ! prognostic vars (split, at n step)
    real(8), intent(in)  :: rhog_split_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx_split    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy_split    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz_split    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw_split     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogw_split_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhoge_split     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhoge_split_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: pregd_split     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: pregd_split_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: rhog     (ADM_gall,   ADM_kall,ADM_lall   ) ! prognostic vars
    real(8), intent(in)  :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: eth      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: eth_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: grhog       (ADM_gall,   ADM_kall,ADM_lall   ) ! large step tendency
    real(8), intent(in)  :: grhog_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: grhogw      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: grhogw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: grhoge      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: grhoge_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: grhogetot   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: grhogetot_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: dt

    real(8) :: drhog_split       (ADM_gall,   ADM_kall,ADM_lall   ) ! source term at t=n+1
    real(8) :: drhog_split_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhoge_split      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhoge_split_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogetot_split   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhogetot_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: drhog1    (ADM_gall,   ADM_kall,ADM_lall   ) ! source term ( large step + t=n+1 )
    real(8) :: drhog1_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhoge1   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: drhoge1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dpre1     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: dpre1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: wk_rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: wk_rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: wk_rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: wk_rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: wk_rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: wk_rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: wk_rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: wk_rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: wk_rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: wk_rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: rhogkin0   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogkin0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogkin1   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogkin1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogkin2   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogkin2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ethtot     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: ethtot_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++vi_path2')

    ! calc rhogkin ( previous )
    call cnvvar_rhokin_ijkl( rhog,     rhog_pl,    & !--- [IN]
                             rhogvx,   rhogvx_pl,  & !--- [IN]
                             rhogvy,   rhogvy_pl,  & !--- [IN]
                             rhogvz,   rhogvz_pl,  & !--- [IN]
                             rhogw,    rhogw_pl,   & !--- [IN]
                             rhogkin0, rhogkin0_pl ) !--- [OUT]

    ! prognostic variables ( large step + split (t=n) )
    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do ij = 1, ADM_gall
             wk_rhog  (ij,k,l) = rhog  (ij,k,l) + rhog_split  (ij,k,l)
             wk_rhogvx(ij,k,l) = rhogvx(ij,k,l) + rhogvx_split(ij,k,l)
             wk_rhogvy(ij,k,l) = rhogvy(ij,k,l) + rhogvy_split(ij,k,l)
             wk_rhogvz(ij,k,l) = rhogvz(ij,k,l) + rhogvz_split(ij,k,l)
             wk_rhogw (ij,k,l) = rhogw (ij,k,l) + rhogw_split (ij,k,l)
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                wk_rhog_pl  (ij,k,l) = rhog_pl  (ij,k,l) + rhog_split_pl  (ij,k,l)
                wk_rhogvx_pl(ij,k,l) = rhogvx_pl(ij,k,l) + rhogvx_split_pl(ij,k,l)
                wk_rhogvy_pl(ij,k,l) = rhogvy_pl(ij,k,l) + rhogvy_split_pl(ij,k,l)
                wk_rhogvz_pl(ij,k,l) = rhogvz_pl(ij,k,l) + rhogvz_split_pl(ij,k,l)
                wk_rhogw_pl (ij,k,l) = rhogw_pl (ij,k,l) + rhogw_split_pl (ij,k,l)
             enddo
          enddo
       enddo
    endif

    ! calc rhogkin ( previous + split(t=n) )
    call cnvvar_rhokin_ijkl( wk_rhog,   wk_rhog_pl,   & !--- [IN]
                             wk_rhogvx, wk_rhogvx_pl, & !--- [IN]
                             wk_rhogvy, wk_rhogvy_pl, & !--- [IN]
                             wk_rhogvz, wk_rhogvz_pl, & !--- [IN]
                             wk_rhogw,  wk_rhogw_pl,  & !--- [IN]
                             rhogkin1,  rhogkin1_pl   ) !--- [OUT]

    !---------------------------------------------------------------------------
    ! update drhog & drhoge
    !---------------------------------------------------------------------------

    if ( TIME_SPLIT ) then
       !--- horizontal flux convergence
       call src_flux_convergence( rhogvx_split2, rhogvx_split2_pl, & !--- [IN]
                                  rhogvy_split2, rhogvy_split2_pl, & !--- [IN]
                                  rhogvz_split2, rhogvz_split2_pl, & !--- [IN]
                                  rhogw_split,   rhogw_split_pl,   & !--- [IN]
                                  drhog_split,   drhog_split_pl,   & !--- [OUT]
                                  I_SRC_horizontal                 ) !--- [IN]
    else
       drhog_split   (:,:,:) = 0.D0
       drhog_split_pl(:,:,:) = 0.D0
    endif

    if ( TIME_SPLIT ) then
       !--- horizontal advection convergence
       call src_advection_convergence( rhogvx_split2, rhogvx_split2_pl, & !--- [IN]
                                       rhogvy_split2, rhogvy_split2_pl, & !--- [IN]
                                       rhogvz_split2, rhogvz_split2_pl, & !--- [IN]
                                       rhogw_split,   rhogw_split_pl,   & !--- [IN]
                                       eth,           eth_pl,           & !--- [IN]
                                       drhoge_split,  drhoge_split_pl,  & !--- [OUT]
                                       I_SRC_horizontal                 ) !--- [IN]
    else
       drhoge_split   (:,:,:) = 0.D0
       drhoge_split_pl(:,:,:) = 0.D0
    endif

    !--- update drhog, drhoge, and calc source term of pressure
    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do ij = 1, ADM_gall
             drhog1(ij,k,l)  = grhog(ij,k,l) + drhog_split(ij,k,l)
             drhoge1(ij,k,l) = grhoge(ij,k,l) + drhoge_split(ij,k,l)
             dpre1(ij,k,l)   = drhoge1(ij,k,l) * CNST_RAIR / CNST_CV
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                drhog1_pl(ij,k,l)  = grhog_pl(ij,k,l) + drhog_split_pl(ij,k,l)
                drhoge1_pl(ij,k,l) = grhoge_pl(ij,k,l) + drhoge_split_pl(ij,k,l)
                dpre1_pl(ij,k,l)   = drhoge1_pl(ij,k,l) * CNST_RAIR / CNST_CV
             enddo
          enddo
       enddo
    endif

    !---------------------------------------------------------------------------
    ! verical implict calculation
    !---------------------------------------------------------------------------

    !------ boundary condition for rhogw_split2
    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do ij = 1, ADM_gall
             rhogw_split2(ij,k,l) = 0.D0 ! [add] S.Iga 090827
          enddo
       enddo

       call BNDCND_rhow( ADM_gall,               & !--- [IN]
                         rhogvx_split2(:,:,l),   & !--- [IN]
                         rhogvy_split2(:,:,l),   & !--- [IN]
                         rhogvz_split2(:,:,l),   & !--- [IN]
                         rhogw_split2 (:,:,l),   & !--- [INOUT]
                         VMTR_C2Wfact (:,:,:,l)  ) !--- [IN]
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                rhogw_split2_pl(ij,k,l) = 0.D0 ! [add] S.Iga 090827
             enddo
          enddo

          call BNDCND_rhow( ADM_gall_pl,               & !--- [IN]
                            rhogvx_split2_pl(:,:,l),   & !--- [IN]
                            rhogvy_split2_pl(:,:,l),   & !--- [IN]
                            rhogvz_split2_pl(:,:,l),   & !--- [IN]
                            rhogw_split2_pl (:,:,l),   & !--- [INOUT]
                            VMTR_C2Wfact_pl (:,:,:,l)  ) !--- [IN]
       enddo
    endif

    !---< vertical implicit : solved by tridiagonal matrix
    call vi_rhow( rhogw_split2, rhogw_split2_pl, & !--- [INOUT]
                  rhogw_split,  rhogw_split_pl,  & !--- [IN]
                  pregd_split,  pregd_split_pl,  & !--- [IN]
                  rhog_split,   rhog_split_pl,   & !--- [IN]
                  drhog1,       drhog1_pl,       & !--- [IN]
                  grhogw,       grhogw_pl,       & !--- [IN]
                  dpre1,        dpre1_pl,        & !--- [IN]
                  dt                             ) !--- [IN]

    !--- < rhog integration > ---
    call src_flux_convergence( rhogvx_split2, rhogvx_split2_pl, & !--- [IN]
                               rhogvy_split2, rhogvy_split2_pl, & !--- [IN]
                               rhogvz_split2, rhogvz_split2_pl, & !--- [IN]
                               rhogw_split2,  rhogw_split2_pl,  & !--- [IN]
                               drhog_split,   drhog_split_pl,   & !--- [OUT]
                               I_SRC_default                    ) !--- [IN] [mod] H.Yashiro 20120530

    !------ update drhog & rhog_split2
    rhog_split2(:,:,:) = rhog_split(:,:,:) + ( grhog(:,:,:) + drhog_split(:,:,:) ) * dt

    if ( ADM_prc_me == ADM_prc_pl ) then
       rhog_split2_pl(:,:,:) = rhog_split_pl(:,:,:) + ( grhog_pl(:,:,:) + drhog_split_pl(:,:,:) ) * dt
    endif

    !---------------------------------------------------------------------------
    ! energy correction by ETOT scheme (Satoh,2002)
    !---------------------------------------------------------------------------

    ! prognostic variables ( large step + split (t=n+1) )
    wk_rhog  (:,:,:) = rhog  (:,:,:) + rhog_split2  (:,:,:)
    wk_rhogvx(:,:,:) = rhogvx(:,:,:) + rhogvx_split2(:,:,:)
    wk_rhogvy(:,:,:) = rhogvy(:,:,:) + rhogvy_split2(:,:,:)
    wk_rhogvz(:,:,:) = rhogvz(:,:,:) + rhogvz_split2(:,:,:)
    wk_rhogw (:,:,:) = rhogw (:,:,:) + rhogw_split2 (:,:,:)

    if ( ADM_prc_me == ADM_prc_pl ) then
       wk_rhog_pl  (:,:,:) = rhog_pl  (:,:,:) + rhog_split2_pl  (:,:,:)
       wk_rhogvx_pl(:,:,:) = rhogvx_pl(:,:,:) + rhogvx_split2_pl(:,:,:)
       wk_rhogvy_pl(:,:,:) = rhogvy_pl(:,:,:) + rhogvy_split2_pl(:,:,:)
       wk_rhogvz_pl(:,:,:) = rhogvz_pl(:,:,:) + rhogvz_split2_pl(:,:,:)
       wk_rhogw_pl (:,:,:) = rhogw_pl (:,:,:) + rhogw_split2_pl (:,:,:)
    endif

    ! calc rhogkin ( previous + split(t=n+1) )
    call cnvvar_rhokin_ijkl( wk_rhog,   wk_rhog_pl,   & !--- [IN]
                             wk_rhogvx, wk_rhogvx_pl, & !--- [IN]
                             wk_rhogvy, wk_rhogvy_pl, & !--- [IN]
                             wk_rhogvz, wk_rhogvz_pl, & !--- [IN]
                             wk_rhogw,  wk_rhogw_pl,  & !--- [IN]
                             rhogkin2,  rhogkin2_pl   ) !--- [OUT]

    !--- calculate ( h + v^{2}/2 + phi )
    ethtot(:,:,:) = eth(:,:,:)                    &
                  + rhogkin0(:,:,:) / rhog(:,:,:) &
                  + phi(:,:,:)

    if ( ADM_prc_me == ADM_prc_pl ) then    
       ethtot_pl(:,:,:) = eth_pl(:,:,:)                       &
                        + rhogkin0_pl(:,:,:) / rhog_pl(:,:,:) &
                        + phi_pl(:,:,:)
    endif

    !--- advection convergence for eth + kin + phi
    call src_advection_convergence( wk_rhogvx,       wk_rhogvx_pl,       & !--- [IN]
                                    wk_rhogvy,       wk_rhogvy_pl,       & !--- [IN]
                                    wk_rhogvz,       wk_rhogvz_pl,       & !--- [IN]
                                    wk_rhogw,        wk_rhogw_pl,        & !--- [IN]
                                    ethtot,          ethtot_pl,          & !--- [IN]
                                    drhogetot_split, drhogetot_split_pl, & !--- [OUT]
                                    I_SRC_default                        ) !--- [IN]

    rhoge_split2(:,:,:) = rhoge_split      (:,:,:)                 & ! t=n
                        + ( grhogetot      (:,:,:)                 & ! tendency of total energy (num.diff+smg+nudge)
                          + drhogetot_split(:,:,:) ) * dt          & ! tendency of total energy (adv.conv.)
                        + ( rhogkin1       (:,:,:)                 & ! kinetic   energy (t=n)
                          - rhogkin2       (:,:,:) )               & ! kinetic   energy (t=n+1)
                        + ( rhog_split     (:,:,:)                 & ! potential energy (diff,t=n)
                          - rhog_split2    (:,:,:) ) * phi(:,:,:)    ! potential energy (diff,t=n+1)

    if ( ADM_prc_me == ADM_prc_pl ) then       
       rhoge_split2_pl(:,:,:) = rhoge_split_pl      (:,:,:)                    &
                              + ( grhogetot_pl      (:,:,:)                    &
                                + drhogetot_split_pl(:,:,:) ) * dt             &
                              + ( rhogkin1_pl       (:,:,:)                    &
                                - rhogkin2_pl       (:,:,:) )                  &
                              + ( rhog_split_pl     (:,:,:)                    &
                                - rhog_split2_pl    (:,:,:) ) * phi_pl(:,:,:)
    endif

    call DEBUG_rapend('++++vi_path2')

    return
  end subroutine vi_path2

  !-----------------------------------------------------------------------------
  subroutine vi_rhow_update_matrix( &
       eth,     eth_pl,     & !--- [IN]
       g_tilde, g_tilde_pl, & !--- [IN]
       dt                   ) !--- [IN]
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_EGRAV, &
       CNST_CV,    &
       CNST_RAIR
    use mod_grd, only: &
       GRD_rdgzh, &
       GRD_rdgz,  &
       GRD_cfac,  &
       GRD_dfac
    use mod_vmtr, only: &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_GAM2H,      &
       VMTR_GAM2H_pl,   &
       VMTR_RGAM2H,     &
       VMTR_RGAM2H_pl,  &
       VMTR_RGSH,       &
       VMTR_RGSH_pl
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    implicit none

    real(8), intent(in) :: eth       (ADM_gall,   ADM_kall,ADM_lall   ) ! enthalpy at the h-lev
    real(8), intent(in) :: eth_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: g_tilde   (ADM_gall,   ADM_kall,ADM_lall   ) ! effective gravitation at the h-lev
    real(8), intent(in) :: g_tilde_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: dt

    real(8), save :: GCVovR ! g * Cv / R
    real(8), save :: ACVovR ! alfa * Cv / R

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++vi_rhow_update_matrix')

! Original concept
!
!    A_o(:,:,:) = VMTR_RGSGAM2(:,:,:)
!    A_i(:,:,:) = VMTR_GAM2H(:,:,:) * eth(:,:,:) ! [debug] 20120727 H.Yashiro
!    B  (:,:,:) = g_tilde(:,:,:)
!    C_o(:,:,:) = VMTR_RGAM2H (:,:,:) * ( CNST_CV / CNST_RAIR * CNST_EGRAV )
!    C_i(:,:,:) = 1.D0 / VMTR_RGAM2H(:,:,:)
!    D  (:,:,:) = CNST_CV / CNST_RAIR / ( dt*dt ) / VMTR_RGSH(:,:,:)
!
!    do k = ADM_kmin+1, ADM_kmax
!       Mc(:,k,:) = dble(NON_HYDRO_ALPHA)*D(:,k,:)&
!            + GRD_rdgzh(k)*(                     &
!            + GRD_rdgz(k)*A_o(:,k,:)*A_i(:,k,:)      &
!            + GRD_rdgz(k-1)*A_o(:,k-1,:)*A_i(:,k,:)  &
!            -0.5D0*(GRD_dfac(k)-GRD_cfac(k-1))*(B(:,k,:)+C_o(:,k,:)*C_i(:,k,:))&
!            )
!       Mu(:,k,:) = -GRD_rdgzh(k)*GRD_rdgz(k)*A_o(:,k,:)*A_i(:,k+1,:)&
!            -0.5D0*GRD_rdgzh(k)*GRD_cfac(k)*(B(:,k+1,:)+C_o(:,k,:)*C_i(:,k+1,:))
!       Ml(:,k,:) = -GRD_rdgzh(k)*GRD_rdgz(k)*A_o(:,k,:)*A_i(:,k-1,:)&
!            +0.5D0*GRD_rdgzh(k)*GRD_dfac(k-1)*(B(:,k-1,:)+C_o(:,k,:)*C_i(:,k-1,:))
!    end do

    if ( iflag ) then
       iflag = .false.

       allocate( Mc   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
       allocate( Mu   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
       allocate( Ml   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
       Mc   (:,:,:) = 0.D0
       Mc_pl(:,:,:) = 0.D0
       Mu   (:,:,:) = 0.D0
       Mu_pl(:,:,:) = 0.D0
       Ml   (:,:,:) = 0.D0
       Ml_pl(:,:,:) = 0.D0

       allocate( A2_o     (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( A2_o_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
       allocate( CooCip   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( CooCip_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
       allocate( CooCim   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( CooCim_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
       allocate( D2       (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( D2_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

       GCVovR = CNST_EGRAV * CNST_CV / CNST_RAIR
       alfa   = real(NON_HYDRO_ALPHA,kind=8)
       ACVovR = alfa * CNST_CV / CNST_RAIR

       do l = 1, ADM_lall
          do k = ADM_kmin, ADM_kmax
             do ij = 1, ADM_gall
                A2_o  (ij,k,l) = VMTR_RGSGAM2(ij,k,l) * GRD_rdgz(k)
                CooCip(ij,k,l) = GCVovR * VMTR_RGAM2H(ij,k,l) * VMTR_GAM2H(ij,k+1,l)
                CooCim(ij,k,l) = GCVovR * VMTR_RGAM2H(ij,k,l) * VMTR_GAM2H(ij,k-1,l)
                D2    (ij,k,l) = ACVovR / ( dt*dt ) / VMTR_RGSH(ij,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin, ADM_kmax
                do ij = 1, ADM_gall_pl
                   A2_o_pl  (ij,k,l) = VMTR_RGSGAM2_pl(ij,k,l) * GRD_rdgz(k)
                   CooCip_pl(ij,k,l) = GCVovR * VMTR_RGAM2H_pl(ij,k,l) * VMTR_GAM2H_pl(ij,k+1,l)
                   CooCim_pl(ij,k,l) = GCVovR * VMTR_RGAM2H_pl(ij,k,l) * VMTR_GAM2H_pl(ij,k-1,l)
                   D2_pl    (ij,k,l) = ACVovR / ( dt*dt ) / VMTR_RGSH_pl(ij,k,l)
                enddo
             enddo
          enddo
       endif

    endif ! only once

    do l = 1, ADM_lall
       do k  = ADM_kmin+1, ADM_kmax
       do ij = 1, ADM_gall
          Mc(ij,k,l) = D2(ij,k,l) &
                     + GRD_rdgzh(k) * ( ( A2_o(ij,k,l)+A2_o(ij,k-1,l) )       &
                                      * VMTR_GAM2H(ij,k,l) * eth(ij,k,l)      &
                                      - 0.5D0 * ( GRD_dfac(k)-GRD_cfac(k-1) ) &
                                      * ( g_tilde(ij,k,l)+GCVovR )            )

          Mu(ij,k,l) = -GRD_rdgzh(k) * ( A2_o(ij,k,l)                         &
                                       * VMTR_GAM2H(ij,k+1,l) * eth(ij,k+1,l) &
                                       + 0.5D0 * GRD_cfac(k)                  &
                                       * ( g_tilde(ij,k+1,l)+CooCip(ij,k,l) ) )

          Ml(ij,k,l) = -GRD_rdgzh(k) * ( A2_o(ij,k,l)                         &
                                       * VMTR_GAM2H(ij,k-1,l) * eth(ij,k-1,l) &
                                       - 0.5D0 * GRD_dfac(k-1)                &
                                       * ( g_tilde(ij,k-1,l)+CooCim(ij,k,l) ) )
       enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k  = ADM_kmin+1, ADM_kmax
          do ij = 1, ADM_gall_pl
             Mc_pl(ij,k,l) = D2_pl(ij,k,l) &
                           + GRD_rdgzh(k) * ( ( A2_o_pl(ij,k,l)+A2_o_pl(ij,k-1,l) )  &
                                            * VMTR_GAM2H_pl(ij,k,l) * eth_pl(ij,k,l) &
                                            - 0.5D0 * ( GRD_dfac(k)-GRD_cfac(k-1) )  &
                                            * ( g_tilde_pl(ij,k,l)+GCVovR )          )

             Mu_pl(ij,k,l) = -GRD_rdgzh(k) * ( A2_o_pl(ij,k,l)                            &
                                             * VMTR_GAM2H_pl(ij,k+1,l) * eth_pl(ij,k+1,l) &
                                             + 0.5D0 * GRD_cfac(k)                        &
                                             * ( g_tilde_pl(ij,k+1,l)+CooCip_pl(ij,k,l) ) )

             Ml_pl(ij,k,l) = -GRD_rdgzh(k) * ( A2_o_pl(ij,k,l)                            &
                                             * VMTR_GAM2H_pl(ij,k-1,l) * eth_pl(ij,k-1,l) &
                                             - 0.5D0 * GRD_dfac(k-1)                      &
                                             * ( g_tilde_pl(ij,k-1,l)+CooCim_pl(ij,k,l) ) )
          enddo
          enddo
       enddo
    endif

    call DEBUG_rapend('++++vi_rhow_update_matrix')

    return
  end subroutine vi_rhow_update_matrix

  !-----------------------------------------------------------------------------
  subroutine vi_rhow( & 
       rhogw_new, rhogw_new_pl, & !--- [INOUT]
       rhogw,     rhogw_pl,     & !--- [IN]
       preg,      preg_pl,      & !--- [IN]
       rhog,      rhog_pl,      & !--- [IN]
       Sr,        Sr_pl,        & !--- [IN]
       Sw,        Sw_pl,        & !--- [IN]
       Sp,        Sp_pl,        & !--- [IN]
       dt                       )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgzh, &
       GRD_afac,  &
       GRD_bfac
    use mod_cnst, only: &
       CNST_EGRAV, &
       CNST_CV,    &
       CNST_RAIR
    use mod_vmtr, only: &
       VMTR_RGSGAM2,     &
       VMTR_RGSGAM2_pl,  &
       VMTR_RGSGAM2H,    &
       VMTR_RGSGAM2H_pl, &
       VMTR_RGAM2H,      &
       VMTR_RGAM2H_pl,   &
       VMTR_RGAM2,       &
       VMTR_RGAM2_pl,    &
       VMTR_GSGAM2H,     &
       VMTR_GSGAM2H_pl
    implicit none

    real(8), intent(inout) :: rhogw_new   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w (split) ( gam2 X G^{1/2} )
    real(8), intent(inout) :: rhogw_new_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)    :: rhogw   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w    (split) ( gam2 X G^{1/2} )
    real(8), intent(in)    :: rhogw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: preg    (ADM_gall,   ADM_kall,ADM_lall   ) ! pertub p (split) ( gam2 X G^{1/2} )
    real(8), intent(in)    :: preg_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhog    (ADM_gall,   ADM_kall,ADM_lall   ) ! rhod     (split) ( gam2 X G^{1/2} )
    real(8), intent(in)    :: rhog_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: Sr      (ADM_gall,   ADM_kall,ADM_lall   ) ! source term for rho at the int-lev
    real(8), intent(in)    :: Sr_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: Sw      (ADM_gall,   ADM_kall,ADM_lall   ) ! source term for rhow at the h-lev
    real(8), intent(in)    :: Sw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: Sp      (ADM_gall,   ADM_kall,ADM_lall   ) ! source term for p at the int-lev
    real(8), intent(in)    :: Sp_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: dt

    real(8) :: Sall   (ADM_gall,   ADM_kall)
    real(8) :: Sall_pl(ADM_gall_pl,ADM_kall)

    real(8) :: beta   (ADM_gall   )
    real(8) :: beta_pl(ADM_gall_pl)

    real(8) :: gamma   (ADM_gall,   ADM_kall)
    real(8) :: gamma_pl(ADM_gall_pl,ADM_kall)

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++vi_rhow')

    do l = 1, ADM_lall

       !--- boundary conditions
       do ij = 1, ADM_gall
          rhogw_new(ij,ADM_kmin,  l) = rhogw_new(ij,ADM_kmin,  l) * VMTR_RGSGAM2H(ij,ADM_kmin,  l)
          rhogw_new(ij,ADM_kmax+1,l) = rhogw_new(ij,ADM_kmax+1,l) * VMTR_RGSGAM2H(ij,ADM_kmax+1,l)
       enddo

       do k  = ADM_kmin+1, ADM_kmax
          do ij = 1, ADM_gall
             Sall(ij,k) = (   VMTR_RGAM2H (ij,k,  l)               * ( rhogw(ij,k,  l)*alfa + dt*Sw(ij,k,  l) ) &
                          - ( VMTR_RGSGAM2(ij,k,  l)               * ( preg (ij,k,  l)      + dt*Sp(ij,k,  l) ) &
                            - VMTR_RGSGAM2(ij,k-1,l)               * ( preg (ij,k-1,l)      + dt*Sp(ij,k-1,l) ) &
                            ) * dt * GRD_rdgzh(k)                                                               &
                          - ( VMTR_RGAM2  (ij,k,  l) * GRD_afac(k) * ( rhog (ij,k,  l)      + dt*Sr(ij,k,  l) ) &
                            + VMTR_RGAM2  (ij,k-1,l) * GRD_bfac(k) * ( rhog (ij,k-1,l)      + dt*Sr(ij,k-1,l) ) &
                            ) * dt * 0.5D0 * CNST_EGRAV                                                         &
                          ) * ( CNST_CV / CNST_RAIR / (dt*dt) )
          enddo
       enddo

       do ij = 1, ADM_gall
          Sall(ij,ADM_kmin+1) = Sall(ij,ADM_kmin+1) - Ml(ij,ADM_kmin+1,l) * rhogw_new(ij,ADM_kmin,  l)
          Sall(ij,ADM_kmax  ) = Sall(ij,ADM_kmax  ) - Mu(ij,ADM_kmax,  l) * rhogw_new(ij,ADM_kmax+1,l)
       enddo

       !
       !--- < solve tri-daigonal matrix > ---
       !

       ! condition at ADM_kmin+1
       k = ADM_kmin+1
       do ij = 1, ADM_gall
          beta(ij) = Mc(ij,k,l)
          rhogw_new(ij,k,l) = Sall(ij,k) / beta(ij)
       enddo

       !--- forward
       do k  = ADM_kmin+2, ADM_kmax
          do ij = 1, ADM_gall
             gamma(ij,k) = Mu(ij,k-1,l) / beta(ij)

             beta(ij)    = Mc(ij,k,l) - Ml(ij,k,l) * gamma(ij,k)

             rhogw_new(ij,k,l) = ( Sall(ij,k) - Ml(ij,k,l) * rhogw_new(ij,k-1,l) ) / beta(ij)
          enddo
       enddo

       !--- backward
       do k  = ADM_kmax-1, ADM_kmin+1, -1
          do ij = 1, ADM_gall
             rhogw_new(ij,k,l) = rhogw_new(ij,k,l) - gamma(ij,k+1) * rhogw_new(ij,k+1,l)
          enddo
       enddo

       !--- return value ( gam2 X G^{1/2} )
       do k  = ADM_kmin, ADM_kmax+1
          do ij = 1, ADM_gall
             rhogw_new(ij,k,l) = rhogw_new(ij,k,l) * VMTR_GSGAM2H(ij,k,l)
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          !--- boundary conditions
          do ij = 1, ADM_gall_pl
             rhogw_new_pl(ij,ADM_kmin,  l) = rhogw_new_pl(ij,ADM_kmin,  l) * VMTR_RGSGAM2H_pl(ij,ADM_kmin,  l)
             rhogw_new_pl(ij,ADM_kmax+1,l) = rhogw_new_pl(ij,ADM_kmax+1,l) * VMTR_RGSGAM2H_pl(ij,ADM_kmax+1,l)
          enddo

          do k  = ADM_kmin+1, ADM_kmax
             do ij = 1, ADM_gall_pl
                Sall_pl(ij,k) = (   VMTR_RGAM2H_pl (ij,k,  l)               * ( rhogw_pl(ij,k,  l)*alfa + dt*Sw_pl(ij,k,  l) ) &
                                - ( VMTR_RGSGAM2_pl(ij,k,  l)               * ( preg_pl (ij,k,  l)      + dt*Sp_pl(ij,k,  l) ) &
                                  - VMTR_RGSGAM2_pl(ij,k-1,l)               * ( preg_pl (ij,k-1,l)      + dt*Sp_pl(ij,k-1,l) ) &
                                  ) * dt * GRD_rdgzh(k)                                                               &
                                - ( VMTR_RGAM2_pl  (ij,k,  l) * GRD_afac(k) * ( rhog_pl (ij,k,  l)      + dt*Sr_pl(ij,k,  l) ) &
                                  + VMTR_RGAM2_pl  (ij,k-1,l) * GRD_bfac(k) * ( rhog_pl (ij,k-1,l)      + dt*Sr_pl(ij,k-1,l) ) &
                                  ) * dt * 0.5D0 * CNST_EGRAV                                                                  &
                                ) * ( CNST_CV / CNST_RAIR / (dt*dt) )
             enddo
          enddo

          do ij = 1, ADM_gall_pl
             Sall_pl(ij,ADM_kmin+1) = Sall_pl(ij,ADM_kmin+1) - Ml_pl(ij,ADM_kmin+1,l) * rhogw_new_pl(ij,ADM_kmin,  l)
             Sall_pl(ij,ADM_kmax  ) = Sall_pl(ij,ADM_kmax  ) - Mu_pl(ij,ADM_kmax,  l) * rhogw_new_pl(ij,ADM_kmax+1,l)
          enddo

          !
          !--- < solve tri-daigonal matrix > ---
          !

          ! condition at ADM_kmin+1
          k = ADM_kmin+1
          do ij = 1, ADM_gall_pl
             beta_pl(ij) = Mc_pl(ij,k,l)
             rhogw_new_pl(ij,k,l) = Sall_pl(ij,k) / beta_pl(ij)
          enddo

          !--- forward
          do k  = ADM_kmin+2, ADM_kmax
             do ij = 1, ADM_gall_pl
                gamma_pl(ij,k) = Mu_pl(ij,k-1,l) / beta_pl(ij)

                beta_pl(ij)    = Mc_pl(ij,k,l) - Ml_pl(ij,k,l) * gamma_pl(ij,k)

                rhogw_new_pl(ij,k,l) = ( Sall_pl(ij,k) - Ml_pl(ij,k,l) * rhogw_new_pl(ij,k-1,l) ) / beta_pl(ij)
             enddo
          enddo

          !--- backward
          do k  = ADM_kmax-1, ADM_kmin+1, -1
             do ij = 1, ADM_gall_pl
                rhogw_new_pl(ij,k,l) = rhogw_new_pl(ij,k,l) - gamma_pl(ij,k+1) * rhogw_new_pl(ij,k+1,l)
             enddo
          enddo

          !--- return value ( gam2 X G^{1/2} )
          do k  = ADM_kmin, ADM_kmax+1
             do ij = 1, ADM_gall_pl
                rhogw_new_pl(ij,k,l) = rhogw_new_pl(ij,k,l) * VMTR_GSGAM2H_pl(ij,k,l)
             enddo
          enddo
       enddo

    endif

    call DEBUG_rapend('++++vi_rhow')

    return
  end subroutine vi_rhow

end module mod_vi
!-------------------------------------------------------------------------------

