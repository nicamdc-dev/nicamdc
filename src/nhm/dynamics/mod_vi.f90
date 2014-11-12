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
  public :: vi_setup
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
  real(8), private, allocatable :: Mc   (:,:,:)
  real(8), private, allocatable :: Mc_pl(:,:,:)
  real(8), private, allocatable :: Ml   (:,:,:)
  real(8), private, allocatable :: Ml_pl(:,:,:)
  real(8), private, allocatable :: Mu   (:,:,:)
  real(8), private, allocatable :: Mu_pl(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine vi_setup
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    implicit none
    !---------------------------------------------------------------------------

    allocate( Mc   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Mu   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Ml   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    return
  end subroutine vi_setup

  !-----------------------------------------------------------------------------
  subroutine vi_small_step( &
       rhog,         rhog_pl,         &
       rhogvx,       rhogvx_pl,       &
       rhogvy,       rhogvy_pl,       &
       rhogvz,       rhogvz_pl,       &
       rhogw,        rhogw_pl,        &
       rhoge,        rhoge_pl,        &
       vx,           vx_pl,           &
       vy,           vy_pl,           &
       vz,           vz_pl,           &
       eth,          eth_pl,          &
       rhog_prim,    rhog_prim_pl,    &
       preg_prim,    preg_prim_pl,    &
       grhog0,       grhog0_pl,       &
       grhogvx0,     grhogvx0_pl,     &
       grhogvy0,     grhogvy0_pl,     &
       grhogvz0,     grhogvz0_pl,     &
       grhogw0,      grhogw0_pl,      &
       grhoge0,      grhoge0_pl,      &
       grhogetot0,   grhogetot0_pl,   &
       rhog_split,   rhog_split_pl,   &
       rhogvx_split, rhogvx_split_pl, &
       rhogvy_split, rhogvy_split_pl, &
       rhogvz_split, rhogvz_split_pl, &
       rhogw_split,  rhogw_split_pl,  &
       rhoge_split,  rhoge_split_pl,  &
       v_mean_c,     v_mean_c_pl,     &
       num_of_itr,                    &
       dt                             )
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
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_afac, &
       GRD_bfac
    use mod_oprt, only: &
       OPRT_horizontalize_vec
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl, &
       VMTR_W2Cfact,    &
       VMTR_W2Cfact_pl
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
       src_pres_gradient,         &
       src_buoyancy,              &
       I_SRC_default,             &
       I_SRC_horizontal
    implicit none

    real(8), intent(inout) :: rhog           (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho     ( G^1/2 x gam2 )
    real(8), intent(inout) :: rhog_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvx         (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vx  ( G^1/2 x gam2 )
    real(8), intent(inout) :: rhogvx_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvy         (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vy  ( G^1/2 x gam2 )
    real(8), intent(inout) :: rhogvy_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvz         (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vz  ( G^1/2 x gam2 )
    real(8), intent(inout) :: rhogvz_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogw          (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w   ( G^1/2 x gam2 )
    real(8), intent(inout) :: rhogw_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhoge          (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Ein ( G^1/2 x gam2 )
    real(8), intent(inout) :: rhoge_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in)    :: vx             (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_x
    real(8), intent(in)    :: vx_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy             (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_y
    real(8), intent(in)    :: vy_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz             (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_z
    real(8), intent(in)    :: vz_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: eth            (ADM_gall   ,ADM_kall,ADM_lall   ) ! enthalpy
    real(8), intent(in)    :: eth_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhog_prim      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho prime ( G^1/2 x gam2 )
    real(8), intent(in)    :: rhog_prim_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: preg_prim      (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(8), intent(in)    :: preg_prim_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)    :: grhog0         (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term (large step)
    real(8), intent(in)    :: grhog0_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogvx0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogvx0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogvy0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogvy0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogvz0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogvz0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogw0        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogw0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhoge0        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhoge0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: grhogetot0     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: grhogetot0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(inout) :: rhog_split     (ADM_gall   ,ADM_kall,ADM_lall   ) ! splited value
    real(8), intent(inout) :: rhog_split_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvx_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhogvx_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvy_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhogvy_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogvz_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhogvz_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhogw_split    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhogw_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: rhoge_split    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: rhoge_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(out)   :: v_mean_c       (ADM_gall   ,ADM_kall,ADM_lall   ,5) ! mean_flux for tracer advection
    real(8), intent(out)   :: v_mean_c_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,5)

    integer, intent(in)    :: num_of_itr
    real(8), intent(in)    :: dt

    !--- tendency term (large step + small step)
    real(8) :: grhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: grhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: grhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: grhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: grhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: grhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: grhoge    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: grhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- tendency term 2
    real(8) :: drhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhoge    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- divergence damping
    real(8) :: ddivdvx      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: ddivdvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ddivdvy      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: ddivdvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ddivdvz      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: ddivdvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ddivdw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: ddivdw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: ddivdvx_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: ddivdvx_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ddivdvy_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: ddivdvy_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ddivdvz_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: ddivdvz_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- pressure gradient force
    real(8) :: dpgrad       (ADM_gall   ,ADM_kall,ADM_lall   ,GRD_XDIR:GRD_ZDIR)
    real(8) :: dpgrad_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_XDIR:GRD_ZDIR)
    real(8) :: dpgradw      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dpgradw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- buoyancy force
    real(8) :: dbuoiw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dbuoiw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- pressure work
    real(8) :: drhoge_pw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhoge_pw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhoge_pwh   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhoge_pwh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: gz_tilde     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: gz_tilde_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: rhog_h       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: rhog_h_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: eth_h        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: eth_h_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: preg_prim_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: preg_prim_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- for communication
    real(8) :: v_split2_vh   (ADM_gall   ,ADM_kall,ADM_lall   ,3)
    real(8) :: v_split2_vh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,3)
    real(8) :: v_split2_we   (ADM_gall   ,ADM_kall,ADM_lall   ,3)
    real(8) :: v_split2_we_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,3)

    real(8) :: rweight_itr

    integer :: g, k, l, ns

    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    !--- full level -> half level
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          rhog_h(g,k,l) = ( VMTR_C2Wfact(1,g,k,l) * rhog(g,k,  l) &
                          + VMTR_C2Wfact(2,g,k,l) * rhog(g,k-1,l) )
          eth_h (g,k,l) = 0.5D0 * ( GRD_afac(k) * eth(g,k,  l) &
                                  + GRD_bfac(k) * eth(g,k-1,l) )
       enddo
       enddo
       do g = 1, ADM_gall
          rhog_h(g,ADM_kmin-1,l) = rhog_h(g,ADM_kmin,l)
          eth_h (g,ADM_kmin-1,l) = eth_h (g,ADM_kmin,l)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,k,l) = ( VMTR_C2Wfact_pl(1,g,k,l) * rhog_pl(g,k,  l) &
                                + VMTR_C2Wfact_pl(2,g,k,l) * rhog_pl(g,k-1,l) )
             eth_h_pl (g,k,l) = 0.5D0 * ( GRD_afac(k) * eth_pl(g,k,  l) &
                                        + GRD_bfac(k) * eth_pl(g,k-1,l) )
          enddo
          enddo
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,ADM_kmin-1,l) = rhog_h_pl(g,ADM_kmin,l)
             eth_h_pl (g,ADM_kmin-1,l) = eth_h_pl (g,ADM_kmin,l)
          enddo
       enddo
    endif

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
    call numfilter_divdamp( rhogvx,  rhogvx_pl,  & !--- [IN]
                            rhogvy,  rhogvy_pl,  & !--- [IN]
                            rhogvz,  rhogvz_pl,  & !--- [IN]
                            rhogw,   rhogw_pl,   & !--- [IN]
                            ddivdvx, ddivdvx_pl, & !--- [OUT]
                            ddivdvy, ddivdvy_pl, & !--- [OUT]
                            ddivdvz, ddivdvz_pl, & !--- [OUT]
                            ddivdw,  ddivdw_pl   ) !--- [OUT]

    call numfilter_divdamp_2d( rhogvx,     rhogvx_pl,     & !--- [IN]
                               rhogvy,     rhogvy_pl,     & !--- [IN]
                               rhogvz,     rhogvz_pl,     & !--- [IN]
                               ddivdvx_2d, ddivdvx_2d_pl, & !--- [OUT]
                               ddivdvy_2d, ddivdvy_2d_pl, & !--- [OUT]
                               ddivdvz_2d, ddivdvz_2d_pl  ) !--- [OUT]

    !--- pressure force
    call src_pres_gradient( preg_prim, preg_prim_pl, & !--- [IN]
                            dpgrad,    dpgrad_pl,    & !--- [OUT]
                            dpgradw,   dpgradw_pl,   & !--- [OUT]
                            I_SRC_default            ) !--- [IN]

    !--- buoyancy force
    call src_buoyancy( rhog_prim, rhog_prim_pl, & !--- [IN]
                       dbuoiw,    dbuoiw_pl     ) !--- [OUT]

    !---< Calculation of drhoge

    !--- advection convergence for eth
    call src_advection_convergence( rhogvx, rhogvx_pl, & !--- [IN]
                                    rhogvy, rhogvy_pl, & !--- [IN]
                                    rhogvz, rhogvz_pl, & !--- [IN]
                                    rhogw,  rhogw_pl,  & !--- [IN]
                                    eth,    eth_pl,    & !--- [IN]
                                    drhoge, drhoge_pl, & !--- [OUT]
                                    I_SRC_default      ) !--- [IN]

    !--- pressure work
    do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          gz_tilde  (g,k,l) = CNST_EGRAV - ( dpgradw(g,k,l)-dbuoiw(g,k,l) ) / rhog_h(g,k,l)
          drhoge_pwh(g,k,l) = -gz_tilde(g,k,l) * rhogw(g,k,l)
       enddo
       enddo

       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          drhoge_pw(g,k,l) = ( vx(g,k,l) * dpgrad(g,k,l,GRD_XDIR)          &
                             + vy(g,k,l) * dpgrad(g,k,l,GRD_YDIR)          &
                             + vz(g,k,l) * dpgrad(g,k,l,GRD_ZDIR)          ) &
                           + ( VMTR_W2Cfact(1,g,k,l) * drhoge_pwh(g,k+1,l) &
                             + VMTR_W2Cfact(2,g,k,l) * drhoge_pwh(g,k,  l) )
       enddo
       enddo
       do g = 1, ADM_gall
          drhoge_pw(g,ADM_kmin-1,l) = 0.D0
          drhoge_pw(g,ADM_kmax+1,l) = 0.D0
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             gz_tilde_pl  (g,k,l) = CNST_EGRAV - ( dpgradw_pl(g,k,l)-dbuoiw_pl(g,k,l) ) / rhog_h_pl(g,k,l)
             drhoge_pwh_pl(g,k,l) = -gz_tilde_pl(g,k,l) * rhogw_pl(g,k,l)
          enddo
          enddo

          do k  = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             drhoge_pw_pl(g,k,l) = ( vx_pl(g,k,l) * dpgrad_pl(g,k,l,GRD_XDIR)          &
                                   + vy_pl(g,k,l) * dpgrad_pl(g,k,l,GRD_YDIR)          &
                                   + vz_pl(g,k,l) * dpgrad_pl(g,k,l,GRD_ZDIR)          ) &
                                 + ( VMTR_W2Cfact_pl(1,g,k,l) * drhoge_pwh_pl(g,k+1,l) &
                                   + VMTR_W2Cfact_pl(2,g,k,l) * drhoge_pwh_pl(g,k,  l) )
          enddo
          enddo
          do g = 1, ADM_gall_pl
             drhoge_pw_pl(g,ADM_kmin-1,l) = 0.D0
             drhoge_pw_pl(g,ADM_kmax+1,l) = 0.D0
          enddo
       enddo
    endif

    !--- sum of tendencies ( large step + pres-grad + div-damp + div-damp_2d + buoyancy )
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       grhog  (g,k,l) = grhog0  (g,k,l) + drhog(g,k,l)
       grhogvx(g,k,l) = grhogvx0(g,k,l) - dpgrad(g,k,l,GRD_XDIR) + ddivdvx(g,k,l) + ddivdvx_2d(g,k,l)
       grhogvy(g,k,l) = grhogvy0(g,k,l) - dpgrad(g,k,l,GRD_YDIR) + ddivdvy(g,k,l) + ddivdvy_2d(g,k,l)
       grhogvz(g,k,l) = grhogvz0(g,k,l) - dpgrad(g,k,l,GRD_ZDIR) + ddivdvz(g,k,l) + ddivdvz_2d(g,k,l)
       grhogw (g,k,l) = grhogw0 (g,k,l) - dpgradw(g,k,l) + dbuoiw(g,k,l) + ddivdw(g,k,l) * NON_HYDRO_ALPHA
       grhoge (g,k,l) = grhoge0 (g,k,l) + drhoge(g,k,l) + drhoge_pw(g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do g = 1, ADM_gall_pl
          grhog_pl  (g,k,l) = grhog0_pl  (g,k,l) + drhog_pl(g,k,l)
          grhogvx_pl(g,k,l) = grhogvx0_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_XDIR) + ddivdvx_pl(g,k,l)
          grhogvy_pl(g,k,l) = grhogvy0_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_YDIR) + ddivdvy_pl(g,k,l)
          grhogvz_pl(g,k,l) = grhogvz0_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_ZDIR) + ddivdvz_pl(g,k,l)
          grhogw_pl (g,k,l) = grhogw0_pl (g,k,l) - dpgradw_pl(g,k,l) + dbuoiw_pl(g,k,l) + ddivdw_pl(g,k,l) * NON_HYDRO_ALPHA
          grhoge_pl (g,k,l) = grhoge0_pl (g,k,l) + drhoge_pl(g,k,l) + drhoge_pw_pl(g,k,l)
       enddo
       enddo
       enddo
    endif

    !---------------------------------------------------------------------------
    ! END vi_path0
    !---------------------------------------------------------------------------

    ! update working matrix in mod_rhow
    call vi_rhow_update_matrix( eth_h,    eth_h_pl,    & !--- [IN] : enthalpy at the h-lev
                                gz_tilde, gz_tilde_pl, & !--- [IN] : effective gravitation at the h-lev
                                dt                     ) !--- [IN] : delta t

    !--- initialization of mean mass flux
    rweight_itr = 1.D0 / real(num_of_itr,kind=8)

    v_mean_c(:,:,:,1) = rhogvx(:,:,:)
    v_mean_c(:,:,:,2) = rhogvy(:,:,:)
    v_mean_c(:,:,:,3) = rhogvz(:,:,:)
    v_mean_c(:,:,:,4) = rhog  (:,:,:)
    v_mean_c(:,:,:,5) = rhogw (:,:,:)

    v_mean_c_pl(:,:,:,1) = rhogvx_pl(:,:,:)
    v_mean_c_pl(:,:,:,2) = rhogvy_pl(:,:,:)
    v_mean_c_pl(:,:,:,3) = rhogvz_pl(:,:,:)
    v_mean_c_pl(:,:,:,4) = rhog_pl  (:,:,:)
    v_mean_c_pl(:,:,:,5) = rhogw_pl (:,:,:)

    !---------------------------------------------------------------------------
    !
    !> Start small step iteration
    !
    !---------------------------------------------------------------------------
    do ns = 1, num_of_itr

       !---< calculation of preg_prim(*) from rhog(*) & rhoge(*)
       do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             preg_prim_split(g,k,l) = rhoge_split(g,k,l) * ( CNST_RAIR / CNST_CV )
          enddo
          enddo
          do g = 1, ADM_gall
             preg_prim_split(g,ADM_kmin-1,l) = preg_prim_split(g,ADM_kmin,l)
             preg_prim_split(g,ADM_kmax+1,l) = preg_prim_split(g,ADM_kmax,l)
          enddo

          do g = 1, ADM_gall
             rhoge_split(g,ADM_kmin-1,l) = rhoge_split(g,ADM_kmin,l)
             rhoge_split(g,ADM_kmax+1,l) = rhoge_split(g,ADM_kmax,l)
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                preg_prim_split_pl(g,k,l) = rhoge_split_pl(g,k,l) * ( CNST_RAIR / CNST_CV )
             enddo
             enddo
             do g = 1, ADM_gall_pl
                preg_prim_split_pl(g,ADM_kmin-1,l) = preg_prim_split_pl(g,ADM_kmin,l)
                preg_prim_split_pl(g,ADM_kmax+1,l) = preg_prim_split_pl(g,ADM_kmax,l)
             enddo

             do g = 1, ADM_gall_pl
                rhoge_split_pl(g,ADM_kmin-1,l) = rhoge_split_pl(g,ADM_kmin,l)
                rhoge_split_pl(g,ADM_kmax+1,l) = rhoge_split_pl(g,ADM_kmax,l)
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
                                  ddivdvx,      ddivdvx_pl,      & !--- [OUT]
                                  ddivdvy,      ddivdvy_pl,      & !--- [OUT]
                                  ddivdvz,      ddivdvz_pl,      & !--- [OUT]
                                  ddivdw,       ddivdw_pl        ) !--- [OUT]

          !--- 2d divergence damping
          call numfilter_divdamp_2d( rhogvx_split, rhogvx_split_pl, & !--- [IN]
                                     rhogvy_split, rhogvy_split_pl, & !--- [IN]
                                     rhogvz_split, rhogvz_split_pl, & !--- [IN]
                                     ddivdvx_2d,   ddivdvx_2d_pl,   & !--- [OUT]
                                     ddivdvy_2d,   ddivdvy_2d_pl,   & !--- [OUT]
                                     ddivdvz_2d,   ddivdvz_2d_pl    ) !--- [OUT]

          !--- pressure force
          !--- dpgradw=0.D0 becaude of f_type='HORIZONTAL'.
          call src_pres_gradient( preg_prim_split, preg_prim_split_pl, & !--- [IN]
                                  dpgrad,          dpgrad_pl,          & !--- [OUT]
                                  dpgradw,         dpgradw_pl,         & !--- [OUT]
                                  I_SRC_horizontal                     ) !--- [IN]

          !--- buoyancy force
          !--- not calculated, because this term is implicit.

          !---------------------------------------------------------------------
          ! END vi_path1
          !---------------------------------------------------------------------

          !--- sum of tendency ( large step + pres-grad + div-damp + div-damp_2d )
!OCL SERIAL
          do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             drhogvx(g,k,l) = grhogvx(g,k,l) - dpgrad(g,k,l,GRD_XDIR) + ddivdvx(g,k,l) + ddivdvx_2d(g,k,l)
             drhogvy(g,k,l) = grhogvy(g,k,l) - dpgrad(g,k,l,GRD_YDIR) + ddivdvy(g,k,l) + ddivdvy_2d(g,k,l)
             drhogvz(g,k,l) = grhogvz(g,k,l) - dpgrad(g,k,l,GRD_ZDIR) + ddivdvz(g,k,l) + ddivdvz_2d(g,k,l)
             drhogw (g,k,l) = grhogw (g,k,l)                          + ddivdw (g,k,l) * NON_HYDRO_ALPHA
          enddo
          enddo
          enddo

          if ( ADM_prc_me == ADM_prc_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                drhogvx_pl(g,k,l) = grhogvx_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_XDIR) + ddivdvx_pl(g,k,l) + ddivdvx_2d_pl(g,k,l)
                drhogvy_pl(g,k,l) = grhogvy_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_YDIR) + ddivdvy_pl(g,k,l) + ddivdvy_2d_pl(g,k,l)
                drhogvz_pl(g,k,l) = grhogvz_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_ZDIR) + ddivdvz_pl(g,k,l) + ddivdvz_2d_pl(g,k,l)
                drhogw_pl (g,k,l) = grhogw_pl (g,k,l)                             + ddivdw_pl (g,k,l) * NON_HYDRO_ALPHA
             enddo
             enddo
             enddo
          endif

       else !--- NO-SPLITING
          !---------------------------------------------------------------------
          ! vi_path1 is skipped
          !---------------------------------------------------------------------

          !------ sum of tendency ( large step )
!OCL SERIAL
          do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             drhogvx(g,k,l) = grhogvx(g,k,l)
             drhogvy(g,k,l) = grhogvy(g,k,l)
             drhogvz(g,k,l) = grhogvz(g,k,l)
             drhogw (g,k,l) = grhogw (g,k,l)
          enddo
          enddo
          enddo

          if ( ADM_prc_me == ADM_prc_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                drhogvx_pl(g,k,l) = grhogvx_pl(g,k,l)
                drhogvy_pl(g,k,l) = grhogvy_pl(g,k,l)
                drhogvz_pl(g,k,l) = grhogvz_pl(g,k,l)
                drhogw_pl (g,k,l) = grhogw_pl (g,k,l)
             enddo
             enddo
             enddo
          endif

       endif ! Split/Non-split

       do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             v_split2_vh(g,k,l,1) = rhogvx_split(g,k,l) + dt * drhogvx(g,k,l)
             v_split2_vh(g,k,l,2) = rhogvy_split(g,k,l) + dt * drhogvy(g,k,l)
             v_split2_vh(g,k,l,3) = rhogvz_split(g,k,l) + dt * drhogvz(g,k,l)
          enddo
          enddo

          call BNDCND_rhovxvyvz( ADM_gall,             & !--- [IN]
                                 rhog(:,:,l),          & !--- [IN]
                                 v_split2_vh(:,:,l,1), & !--- [INOUT]
                                 v_split2_vh(:,:,l,2), & !--- [INOUT]
                                 v_split2_vh(:,:,l,3)  ) !--- [INOUT]
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                v_split2_vh_pl(g,k,l,1) = rhogvx_split_pl(g,k,l) + dt * drhogvx_pl(g,k,l)
                v_split2_vh_pl(g,k,l,2) = rhogvy_split_pl(g,k,l) + dt * drhogvy_pl(g,k,l)
                v_split2_vh_pl(g,k,l,3) = rhogvz_split_pl(g,k,l) + dt * drhogvz_pl(g,k,l)
             enddo
             enddo

             call BNDCND_rhovxvyvz( ADM_gall_pl,             & !--- [IN]
                                    rhog_pl(:,:,l),          & !--- [IN]
                                    v_split2_vh_pl(:,:,l,1), & !--- [INOUT]
                                    v_split2_vh_pl(:,:,l,2), & !--- [INOUT]
                                    v_split2_vh_pl(:,:,l,3)  ) !--- [INOUT]
          enddo
       endif

       !--- communication of horizontal momentum
       call COMM_data_transfer( v_split2_vh, v_split2_vh_pl )

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
          v_split2_vh(suf(ADM_gall_1d,1),k,l,1) = v_split2_vh(suf(ADM_gmax+1,ADM_gmin),k,l,1)
          v_split2_vh(suf(1,ADM_gall_1d),k,l,1) = v_split2_vh(suf(ADM_gmin,ADM_gmax+1),k,l,1)
          v_split2_vh(suf(ADM_gall_1d,1),k,l,2) = v_split2_vh(suf(ADM_gmax+1,ADM_gmin),k,l,2)
          v_split2_vh(suf(1,ADM_gall_1d),k,l,2) = v_split2_vh(suf(ADM_gmin,ADM_gmax+1),k,l,2)
          v_split2_vh(suf(ADM_gall_1d,1),k,l,3) = v_split2_vh(suf(ADM_gmax+1,ADM_gmin),k,l,3)
          v_split2_vh(suf(1,ADM_gall_1d),k,l,3) = v_split2_vh(suf(ADM_gmin,ADM_gmax+1),k,l,3)
       enddo
       enddo

       !--- 2nd small step : vertical implicit : rhog, rhogw, preg_prim: next step
       call vi_path2( v_split2_we(:,:,:,1), v_split2_we_pl(:,:,:,1), & !--- [OUT]
                      v_split2_we(:,:,:,2), v_split2_we_pl(:,:,:,2), & !--- [OUT]
                      v_split2_we(:,:,:,3), v_split2_we_pl(:,:,:,3), & !--- [OUT]
                      v_split2_vh(:,:,:,1), v_split2_vh_pl(:,:,:,1), & !--- [IN]
                      v_split2_vh(:,:,:,2), v_split2_vh_pl(:,:,:,2), & !--- [IN]
                      v_split2_vh(:,:,:,3), v_split2_vh_pl(:,:,:,3), & !--- [IN]
                      rhog_split,           rhog_split_pl,           & !--- [IN]
                      rhogvx_split,         rhogvx_split_pl,         & !--- [IN]
                      rhogvy_split,         rhogvy_split_pl,         & !--- [IN]
                      rhogvz_split,         rhogvz_split_pl,         & !--- [IN]
                      rhogw_split,          rhogw_split_pl,          & !--- [IN]
                      rhoge_split,          rhoge_split_pl,          & !--- [IN]
                      preg_prim_split,      preg_prim_split_pl,      & !--- [IN]
                      rhog,                 rhog_pl,                 & !--- [IN]
                      rhogvx,               rhogvx_pl,               & !--- [IN]
                      rhogvy,               rhogvy_pl,               & !--- [IN]
                      rhogvz,               rhogvz_pl,               & !--- [IN]
                      rhogw,                rhogw_pl,                & !--- [IN]
                      eth,                  eth_pl,                  & !--- [IN]
                      grhog,                grhog_pl,                & !--- [IN]
                      drhogw,               drhogw_pl,               & !--- [IN]
                      grhoge,               grhoge_pl,               & !--- [IN]
                      grhogetot0,           grhogetot0_pl,           & !--- [IN]
                      dt                                             ) !--- [IN]

       !--- communication of rhog, rhogw, and rhoge
       call COMM_data_transfer( v_split2_we, v_split2_we_pl )

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
          v_split2_we(suf(ADM_gall_1d,1),k,l,1) = v_split2_we(suf(ADM_gmax+1,ADM_gmin),k,l,1)
          v_split2_we(suf(1,ADM_gall_1d),k,l,1) = v_split2_we(suf(ADM_gmin,ADM_gmax+1),k,l,1)
          v_split2_we(suf(ADM_gall_1d,1),k,l,2) = v_split2_we(suf(ADM_gmax+1,ADM_gmin),k,l,2)
          v_split2_we(suf(1,ADM_gall_1d),k,l,2) = v_split2_we(suf(ADM_gmin,ADM_gmax+1),k,l,2)
          v_split2_we(suf(ADM_gall_1d,1),k,l,3) = v_split2_we(suf(ADM_gmax+1,ADM_gmin),k,l,3)
          v_split2_we(suf(1,ADM_gall_1d),k,l,3) = v_split2_we(suf(ADM_gmin,ADM_gmax+1),k,l,3)
       enddo
       enddo

       do l = 1, ADM_lall
          !--- update for next step
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             rhogvx_split(g,k,l) = v_split2_vh(g,k,l,1)
             rhogvy_split(g,k,l) = v_split2_vh(g,k,l,2)
             rhogvz_split(g,k,l) = v_split2_vh(g,k,l,3)
             rhog_split  (g,k,l) = v_split2_we(g,k,l,1)
             rhogw_split (g,k,l) = v_split2_we(g,k,l,2)
             rhoge_split (g,k,l) = v_split2_we(g,k,l,3)
          enddo
          enddo

          !--- calculation of mean mass flux ( for tracers )
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             v_mean_c(g,k,l,1) = v_mean_c(g,k,l,1) + rhogvx_split(g,k,l) * rweight_itr
             v_mean_c(g,k,l,2) = v_mean_c(g,k,l,2) + rhogvy_split(g,k,l) * rweight_itr
             v_mean_c(g,k,l,3) = v_mean_c(g,k,l,3) + rhogvz_split(g,k,l) * rweight_itr
             v_mean_c(g,k,l,4) = v_mean_c(g,k,l,4) + rhog_split  (g,k,l) * rweight_itr
             v_mean_c(g,k,l,5) = v_mean_c(g,k,l,5) + rhogw_split (g,k,l) * rweight_itr
          enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          !--- update for next step
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                rhogvx_split_pl(g,k,l) = v_split2_vh_pl(g,k,l,1)
                rhogvy_split_pl(g,k,l) = v_split2_vh_pl(g,k,l,2)
                rhogvz_split_pl(g,k,l) = v_split2_vh_pl(g,k,l,3)
                rhog_split_pl  (g,k,l) = v_split2_we_pl(g,k,l,1)
                rhogw_split_pl (g,k,l) = v_split2_we_pl(g,k,l,2)
                rhoge_split_pl (g,k,l) = v_split2_we_pl(g,k,l,3)
             enddo
             enddo

             !--- calculation of mean mass flux ( for tracers )
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                v_mean_c_pl(g,k,l,1) = v_mean_c_pl(g,k,l,1) + rhogvx_split_pl(g,k,l) * rweight_itr
                v_mean_c_pl(g,k,l,2) = v_mean_c_pl(g,k,l,2) + rhogvy_split_pl(g,k,l) * rweight_itr
                v_mean_c_pl(g,k,l,3) = v_mean_c_pl(g,k,l,3) + rhogvz_split_pl(g,k,l) * rweight_itr
                v_mean_c_pl(g,k,l,4) = v_mean_c_pl(g,k,l,4) + rhog_split_pl  (g,k,l) * rweight_itr
                v_mean_c_pl(g,k,l,5) = v_mean_c_pl(g,k,l,5) + rhogw_split_pl (g,k,l) * rweight_itr
             enddo
             enddo
          enddo
       endif

    enddo  ! small step end

    !--- update prognostic variables
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       rhog  (g,k,l) = rhog  (g,k,l) + rhog_split  (g,k,l)
       rhogvx(g,k,l) = rhogvx(g,k,l) + rhogvx_split(g,k,l)
       rhogvy(g,k,l) = rhogvy(g,k,l) + rhogvy_split(g,k,l)
       rhogvz(g,k,l) = rhogvz(g,k,l) + rhogvz_split(g,k,l)
       rhogw (g,k,l) = rhogw (g,k,l) + rhogw_split (g,k,l)
       rhoge (g,k,l) = rhoge (g,k,l) + rhoge_split (g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          rhog_pl  (g,k,l) = rhog_pl  (g,k,l) + rhog_split_pl  (g,k,l)
          rhogvx_pl(g,k,l) = rhogvx_pl(g,k,l) + rhogvx_split_pl(g,k,l)
          rhogvy_pl(g,k,l) = rhogvy_pl(g,k,l) + rhogvy_split_pl(g,k,l)
          rhogvz_pl(g,k,l) = rhogvz_pl(g,k,l) + rhogvz_split_pl(g,k,l)
          rhogw_pl (g,k,l) = rhogw_pl (g,k,l) + rhogw_split_pl (g,k,l)
          rhoge_pl (g,k,l) = rhoge_pl (g,k,l) + rhoge_split_pl (g,k,l)
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
       rhogw_split2,  rhogw_split2_pl,  & !--- [OUT]
       rhoge_split2,  rhoge_split2_pl,  & !--- [OUT]
       rhogvx_split2, rhogvx_split2_pl, & !--- [IN]
       rhogvy_split2, rhogvy_split2_pl, & !--- [IN]
       rhogvz_split2, rhogvz_split2_pl, & !--- [IN]
       rhog_split,    rhog_split_pl,    & !--- [IN]
       rhogvx_split,  rhogvx_split_pl,  & !--- [IN]
       rhogvy_split,  rhogvy_split_pl,  & !--- [IN]
       rhogvz_split,  rhogvz_split_pl,  & !--- [IN]
       rhogw_split,   rhogw_split_pl,   & !--- [IN]
       rhoge_split,   rhoge_split_pl,   & !--- [IN]
       preg_prim_split,   preg_prim_split_pl,   & !--- [IN]
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
       VMTR_C2WfactGz,    &
       VMTR_C2WfactGz_pl
    use mod_src, only: &
       src_flux_convergence,      &
       src_advection_convergence, &
       I_SRC_horizontal,          &
       I_SRC_default
    use mod_cnvvar, only: &
       cnvvar_rhokin_ijkl
    use mod_bsstate, only: &
       phi, phi_pl
    use mod_bndcnd, only: &
       BNDCND_rhow
    implicit none

    real(8), intent(out) :: rhog_split2       (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars (split, at n step)
    real(8), intent(out) :: rhog_split2_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogw_split2      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhogw_split2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhoge_split2      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: rhoge_split2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx_split2     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_split2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy_split2     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_split2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz_split2     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_split2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: rhog_split        (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars (split, at n step)
    real(8), intent(in)  :: rhog_split_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx_split      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_split_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy_split      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_split_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz_split      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_split_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw_split       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogw_split_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhoge_split       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhoge_split_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: preg_prim_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: preg_prim_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: rhog              (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars
    real(8), intent(in)  :: rhog_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw             (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogw_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: eth               (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: eth_pl            (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: grhog             (ADM_gall   ,ADM_kall,ADM_lall   ) ! large step tendency
    real(8), intent(in)  :: grhog_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: grhogw            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: grhogw_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: grhoge            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: grhoge_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: grhogetot         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: grhogetot_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: dt

    real(8) :: drhog_split       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term at t=n+1
    real(8) :: drhog_split_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhoge_split      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhoge_split_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhogetot_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhogetot_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: drhog1            (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term ( large step + t=n+1 )
    real(8) :: drhog1_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: drhoge1           (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: drhoge1_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dpre1             (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dpre1_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: wk_rhog           (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: wk_rhog_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: wk_rhogvx         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: wk_rhogvx_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: wk_rhogvy         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: wk_rhogvy_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: wk_rhogvz         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: wk_rhogvz_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: wk_rhogw          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: wk_rhogw_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: rhogkin0          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: rhogkin0_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogkin1          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: rhogkin1_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogkin2          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: rhogkin2_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ethtot            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: ethtot_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____vi_path2')

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
    do g = 1, ADM_gall
       wk_rhog  (g,k,l) = rhog  (g,k,l) + rhog_split  (g,k,l)
       wk_rhogvx(g,k,l) = rhogvx(g,k,l) + rhogvx_split(g,k,l)
       wk_rhogvy(g,k,l) = rhogvy(g,k,l) + rhogvy_split(g,k,l)
       wk_rhogvz(g,k,l) = rhogvz(g,k,l) + rhogvz_split(g,k,l)
       wk_rhogw (g,k,l) = rhogw (g,k,l) + rhogw_split (g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          wk_rhog_pl  (g,k,l) = rhog_pl  (g,k,l) + rhog_split_pl  (g,k,l)
          wk_rhogvx_pl(g,k,l) = rhogvx_pl(g,k,l) + rhogvx_split_pl(g,k,l)
          wk_rhogvy_pl(g,k,l) = rhogvy_pl(g,k,l) + rhogvy_split_pl(g,k,l)
          wk_rhogvz_pl(g,k,l) = rhogvz_pl(g,k,l) + rhogvz_split_pl(g,k,l)
          wk_rhogw_pl (g,k,l) = rhogw_pl (g,k,l) + rhogw_split_pl (g,k,l)
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

       !--- horizontal advection convergence
       call src_advection_convergence( rhogvx_split2, rhogvx_split2_pl, & !--- [IN]
                                       rhogvy_split2, rhogvy_split2_pl, & !--- [IN]
                                       rhogvz_split2, rhogvz_split2_pl, & !--- [IN]
                                       rhogw_split,   rhogw_split_pl,   & !--- [IN]
                                       eth,           eth_pl,           & !--- [IN]
                                       drhoge_split,  drhoge_split_pl,  & !--- [OUT]
                                       I_SRC_horizontal                 ) !--- [IN]
    else
       drhog_split    (:,:,:) = 0.D0
       drhog_split_pl (:,:,:) = 0.D0
       drhoge_split   (:,:,:) = 0.D0
       drhoge_split_pl(:,:,:) = 0.D0
    endif

    !--- update drhog, drhoge, and calc source term of pressure
    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       drhog1 (g,k,l) = grhog  (g,k,l) + drhog_split (g,k,l)
       drhoge1(g,k,l) = grhoge (g,k,l) + drhoge_split(g,k,l)
       dpre1  (g,k,l) = drhoge1(g,k,l) * CNST_RAIR / CNST_CV
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          drhog1_pl (g,k,l) = grhog_pl  (g,k,l) + drhog_split_pl (g,k,l)
          drhoge1_pl(g,k,l) = grhoge_pl (g,k,l) + drhoge_split_pl(g,k,l)
          dpre1_pl  (g,k,l) = drhoge1_pl(g,k,l) * CNST_RAIR / CNST_CV
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
       do g = 1, ADM_gall
          rhogw_split2(g,k,l) = 0.D0
       enddo
       enddo

       call BNDCND_rhow( ADM_gall,               & !--- [IN]
                         rhogvx_split2 (:,:,l),  & !--- [IN]
                         rhogvy_split2 (:,:,l),  & !--- [IN]
                         rhogvz_split2 (:,:,l),  & !--- [IN]
                         rhogw_split2  (:,:,l),  & !--- [INOUT]
                         VMTR_C2WfactGz(:,:,:,l) ) !--- [IN]
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             rhogw_split2_pl(g,k,l) = 0.D0
          enddo
          enddo

          call BNDCND_rhow( ADM_gall_pl,                & !--- [IN]
                            rhogvx_split2_pl (:,:,l),   & !--- [IN]
                            rhogvy_split2_pl (:,:,l),   & !--- [IN]
                            rhogvz_split2_pl (:,:,l),   & !--- [IN]
                            rhogw_split2_pl  (:,:,l),   & !--- [INOUT]
                            VMTR_C2WfactGz_pl(:,:,:,l)  ) !--- [IN]
       enddo
    endif

    !---< vertical implicit : solved by tridiagonal matrix
    call vi_rhow( rhogw_split2,    rhogw_split2_pl,    & !--- [INOUT]
                  rhogw_split,     rhogw_split_pl,     & !--- [IN]
                  preg_prim_split, preg_prim_split_pl, & !--- [IN]
                  rhog_split,      rhog_split_pl,      & !--- [IN]
                  drhog1,          drhog1_pl,          & !--- [IN]
                  grhogw,          grhogw_pl,          & !--- [IN]
                  dpre1,           dpre1_pl,           & !--- [IN]
                  dt                                   ) !--- [IN]

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
    ! energy correction by Etotal (Satoh,2002)
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

    call DEBUG_rapend('____vi_path2')

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
       VMTR_RGAMH,      &
       VMTR_RGAMH_pl,   &
       VMTR_RGSQRTH,    &
       VMTR_RGSQRTH_pl
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    implicit none

    real(8), intent(in) :: eth       (ADM_gall   ,ADM_kall,ADM_lall   ) ! enthalpy at the half lev
    real(8), intent(in) :: eth_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: g_tilde   (ADM_gall   ,ADM_kall,ADM_lall   ) ! effective gravitation at the half lev
    real(8), intent(in) :: g_tilde_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: dt

    real(8) :: GCVovR   ! g * Cv / R
    real(8) :: ACVovRt2 ! alfa * Cv / R / dt**2

    integer :: g, k, l
    !---------------------------------------------------------------------------
    ! Original concept
    !
    ! A_o(:,:,:) = VMTR_RGSGAM2(:,:,:)
    ! A_i(:,:,:) = VMTR_GAM2H(:,:,:) * eth(:,:,:) ! [debug] 20120727 H.Yashiro
    ! B  (:,:,:) = g_tilde(:,:,:)
    ! C_o(:,:,:) = VMTR_RGAM2H (:,:,:) * ( CNST_CV / CNST_RAIR * CNST_EGRAV )
    ! C_i(:,:,:) = 1.D0 / VMTR_RGAM2H(:,:,:)
    ! D  (:,:,:) = CNST_CV / CNST_RAIR / ( dt*dt ) / VMTR_RGSQRTH(:,:,:)
    !
    ! do k = ADM_kmin+1, ADM_kmax
    !    Mc(:,k,:) = dble(NON_HYDRO_ALPHA) *D(:,k,:)              &
    !              + GRD_rdgzh(k)                                 &
    !              * ( GRD_rdgz (k)   * A_o(:,k  ,:) * A_i(:,k,:) &
    !                + GRD_rdgz (k-1) * A_o(:,k-1,:) * A_i(:,k,:) &
    !                - 0.5D0 * ( GRD_dfac(k) - GRD_cfac(k-1) )    &
    !                * ( B(:,k,:) + C_o(:,k,:) * C_i(:,k,:) )     &
    !                )
    !    Mu(:,k,:) = - GRD_rdgzh(k) * GRD_rdgz(k) * A_o(:,k,:) * A_i(:,k+1,:) &
    !                - GRD_rdgzh(k) * 0.5D0 * GRD_cfac(k)                     &
    !                * ( B(:,k+1,:) + C_o(:,k,:) * C_i(:,k+1,:) )
    !    Ml(:,k,:) = - GRD_rdgzh(k) * GRD_rdgz(k) * A_o(:,k,:) * A_i(:,k-1,:) &
    !                + GRD_rdgzh(k) * 0.5D0 * GRD_dfac(k-1)                   &
    !                * ( B(:,k-1,:) + C_o(:,k,:) * C_i(:,k-1,:) )
    ! end do

    call DEBUG_rapstart('____vi_rhow_update_matrix')

    GCVovR   = CNST_EGRAV * CNST_CV / CNST_RAIR
    ACVovRt2 = real(NON_HYDRO_ALPHA,kind=8) * CNST_CV / CNST_RAIR / ( dt*dt )

    do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Mc(g,k,l) = ACVovRt2 / VMTR_RGSQRTH(g,k,l)                             &
                    + GRD_rdgzh(k) * ( ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k  )   &
                                       + VMTR_RGSGAM2(g,k-1,l) * GRD_rdgz(k-1) ) &
                                       * VMTR_GAM2H  (g,k  ,l) * eth(g,k,l)      &
                                     - 0.5D0 * ( GRD_dfac(k) - GRD_cfac(k-1) )   &
                                     * ( g_tilde(g,k,l) + GCVovR )               )
       enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Mu(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k)                     &
                                      * VMTR_GAM2H  (g,k+1,l) * eth(g,k+1,l)                    &
                                      + 0.5D0 * GRD_cfac(k)                                     &
                                      * ( g_tilde   (g,k+1,l)                                   &
                                        + VMTR_GAM2H(g,k+1,l)* VMTR_RGAMH(g,k,l)**2 * GCVovR  ) )
       enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Ml(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k)                     &
                                      * VMTR_GAM2H  (g,k-1,l) * eth(g,k-1,l)                    &
                                      - 0.5D0 * GRD_dfac(k-1)                                   &
                                      * ( g_tilde   (g,k-1,l)                                   &
                                        + VMTR_GAM2H(g,k-1,l) * VMTR_RGAMH(g,k,l)**2 * GCVovR ) )
       enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Mc_pl(g,k,l) = ACVovRt2 / VMTR_RGSQRTH_pl(g,k,l)                             &
                          + GRD_rdgzh(k) * ( ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k  )   &
                                             + VMTR_RGSGAM2_pl(g,k-1,l) * GRD_rdgz(k-1) ) &
                                             * VMTR_GAM2H_pl  (g,k  ,l) * eth_pl(g,k,l)   &
                                           - 0.5D0 * ( GRD_dfac(k) - GRD_cfac(k-1) )      &
                                           * ( g_tilde_pl(g,k,l) + GCVovR )               )
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Mu_pl(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k)                        &
                                            * VMTR_GAM2H_pl  (g,k+1,l) * eth_pl(g,k+1,l)                    &
                                            + 0.5D0 * GRD_cfac(k)                                           &
                                            * ( g_tilde_pl   (g,k+1,l)                                      &
                                              + VMTR_GAM2H_pl(g,k+1,l)* VMTR_RGAMH_pl(g,k,l)**2 * GCVovR  ) )
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Ml_pl(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k)                        &
                                            * VMTR_GAM2H_pl  (g,k-1,l) * eth_pl(g,k-1,l)                    &
                                            - 0.5D0 * GRD_dfac(k-1)                                         &
                                            * ( g_tilde_pl   (g,k-1,l)                                      &
                                              + VMTR_GAM2H_pl(g,k-1,l) * VMTR_RGAMH_pl(g,k,l)**2 * GCVovR ) )
          enddo
          enddo
       enddo
    endif

    call DEBUG_rapend('____vi_rhow_update_matrix')

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
       VMTR_RGAMH,       &
       VMTR_RGAMH_pl,    &
       VMTR_RGAM,        &
       VMTR_RGAM_pl,     &
       VMTR_GSGAM2H,     &
       VMTR_GSGAM2H_pl
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    implicit none

    real(8), intent(inout) :: rhogw_new   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 )
    real(8), intent(inout) :: rhogw_new_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)    :: rhogw       (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 )
    real(8), intent(in)    :: rhogw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: preg        (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(8), intent(in)    :: preg_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhog        (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho            ( G^1/2 x gam2 )
    real(8), intent(in)    :: rhog_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: Sr          (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rho  at the full level
    real(8), intent(in)    :: Sr_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: Sw          (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rhow at the half level
    real(8), intent(in)    :: Sw_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: Sp          (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for pres at the full level
    real(8), intent(in)    :: Sp_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: dt

    real(8) :: Sall    (ADM_gall,   ADM_kall)
    real(8) :: Sall_pl (ADM_gall_pl,ADM_kall)
    real(8) :: beta    (ADM_gall   )
    real(8) :: beta_pl (ADM_gall_pl)
    real(8) :: gamma   (ADM_gall,   ADM_kall)
    real(8) :: gamma_pl(ADM_gall_pl,ADM_kall)

    real(8) :: alfa
    real(8) :: CVovRt2   ! Cv / R / dt**2

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____vi_rhow')

    alfa = real(NON_HYDRO_ALPHA,kind=8)
    CVovRt2 = CNST_CV / CNST_RAIR / (dt*dt)

    do l = 1, ADM_lall
       !--- < calc Sall > ---
       do k  = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Sall(g,k) = (   ( rhogw(g,k,  l)*alfa + dt*Sw(g,k,  l) ) * VMTR_RGAMH  (g,k,  l)**2               &
                      - ( ( preg (g,k,  l)      + dt*Sp(g,k,  l) ) * VMTR_RGSGAM2(g,k,  l)                  &
                        - ( preg (g,k-1,l)      + dt*Sp(g,k-1,l) ) * VMTR_RGSGAM2(g,k-1,l)                  &
                        ) * dt * GRD_rdgzh(k)                                                               &
                      - ( ( rhog (g,k,  l)      + dt*Sr(g,k,  l) ) * VMTR_RGAM   (g,k,  l)**2 * GRD_afac(k) &
                        + ( rhog (g,k-1,l)      + dt*Sr(g,k-1,l) ) * VMTR_RGAM   (g,k-1,l)**2 * GRD_bfac(k) &
                        ) * dt * 0.5D0 * CNST_EGRAV                                                         &
                      ) * CVovRt2
       enddo
       enddo

       !--- boundary conditions
       do g = 1, ADM_gall
          rhogw_new(g,ADM_kmin,  l) = rhogw_new(g,ADM_kmin,  l) * VMTR_RGSGAM2H(g,ADM_kmin,  l)
          rhogw_new(g,ADM_kmax+1,l) = rhogw_new(g,ADM_kmax+1,l) * VMTR_RGSGAM2H(g,ADM_kmax+1,l)
          Sall     (g,ADM_kmin+1)   = Sall(g,ADM_kmin+1) - Ml(g,ADM_kmin+1,l) * rhogw_new(g,ADM_kmin,  l)
          Sall     (g,ADM_kmax  )   = Sall(g,ADM_kmax  ) - Mu(g,ADM_kmax,  l) * rhogw_new(g,ADM_kmax+1,l)
       enddo

       !--- < solve tri-daigonal matrix > ---

       ! condition at ADM_kmin+1
       k = ADM_kmin+1
       do g = 1, ADM_gall
          beta(g)          = Mc(g,k,l)

          rhogw_new(g,k,l) = Sall(g,k) / beta(g)
       enddo

       !--- forward
       do k = ADM_kmin+2, ADM_kmax
       do g = 1, ADM_gall
          gamma(g,k)       = Mu(g,k-1,l) / beta(g)

          beta(g)          = Mc(g,k,l) - Ml(g,k,l) * gamma(g,k) ! update beta

          rhogw_new(g,k,l) = ( Sall(g,k) - Ml(g,k,l) * rhogw_new(g,k-1,l) ) / beta(g)
       enddo
       enddo

       !--- backward
       do k = ADM_kmax-1, ADM_kmin+1, -1
       do g = 1, ADM_gall
          rhogw_new(g,k,l) = rhogw_new(g,k,l) - gamma(g,k+1) * rhogw_new(g,k+1,l)
       enddo
       enddo

       !--- return value ( G^1/2 x gam2 )
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          rhogw_new(g,k,l) = rhogw_new(g,k,l) * VMTR_GSGAM2H(g,k,l)
       enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          !--- < calc Sall > ---
          do k  = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Sall_pl(g,k) = (   ( rhogw_pl(g,k,  l)*alfa + dt*Sw_pl(g,k,  l) ) * VMTR_RGAMH_pl  (g,k,  l)**2               &
                            - ( ( preg_pl (g,k,  l)      + dt*Sp_pl(g,k,  l) ) * VMTR_RGSGAM2_pl(g,k,  l)                  &
                              - ( preg_pl (g,k-1,l)      + dt*Sp_pl(g,k-1,l) ) * VMTR_RGSGAM2_pl(g,k-1,l)                  &
                              ) * dt * GRD_rdgzh(k)                                                                        &
                            - ( ( rhog_pl (g,k,  l)      + dt*Sr_pl(g,k,  l) ) * VMTR_RGAM_pl   (g,k,  l)**2 * GRD_afac(k) &
                              + ( rhog_pl (g,k-1,l)      + dt*Sr_pl(g,k-1,l) ) * VMTR_RGAM_pl   (g,k-1,l)**2 * GRD_bfac(k) &
                              ) * dt * 0.5D0 * CNST_EGRAV                                                                  &
                            ) * CVovRt2
          enddo
          enddo

          !--- boundary conditions
          do g = 1, ADM_gall_pl
             rhogw_new_pl(g,ADM_kmin,  l) = rhogw_new_pl(g,ADM_kmin,  l) * VMTR_RGSGAM2H_pl(g,ADM_kmin,  l)
             rhogw_new_pl(g,ADM_kmax+1,l) = rhogw_new_pl(g,ADM_kmax+1,l) * VMTR_RGSGAM2H_pl(g,ADM_kmax+1,l)
             Sall_pl     (g,ADM_kmin+1)   = Sall_pl(g,ADM_kmin+1) - Ml_pl(g,ADM_kmin+1,l) * rhogw_new_pl(g,ADM_kmin,  l)
             Sall_pl     (g,ADM_kmax  )   = Sall_pl(g,ADM_kmax  ) - Mu_pl(g,ADM_kmax,  l) * rhogw_new_pl(g,ADM_kmax+1,l)
          enddo

          !--- < solve tri-daigonal matrix > ---

          ! condition at ADM_kmin+1
          k = ADM_kmin+1
          do g = 1, ADM_gall_pl
             beta_pl(g)          = Mc_pl(g,k,l)

             rhogw_new_pl(g,k,l) = Sall_pl(g,k) / beta_pl(g)
          enddo

          !--- forward
          do k = ADM_kmin+2, ADM_kmax
          do g = 1, ADM_gall_pl
             gamma_pl(g,k)       = Mu_pl(g,k-1,l) / beta_pl(g)

             beta_pl(g)          = Mc_pl(g,k,l) - Ml_pl(g,k,l) * gamma_pl(g,k) ! update beta

             rhogw_new_pl(g,k,l) = ( Sall_pl(g,k) - Ml_pl(g,k,l) * rhogw_new_pl(g,k-1,l) ) / beta_pl(g)
          enddo
          enddo

          !--- backward
          do k = ADM_kmax-1, ADM_kmin+1, -1
          do g = 1, ADM_gall_pl
             rhogw_new_pl(g,k,l) = rhogw_new_pl(g,k,l) - gamma_pl(g,k+1) * rhogw_new_pl(g,k+1,l)
          enddo
          enddo

          !--- return value ( G^1/2 x gam2 )
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             rhogw_new_pl(g,k,l) = rhogw_new_pl(g,k,l) * VMTR_GSGAM2H_pl(g,k,l)
          enddo
          enddo
       enddo
    endif

    call DEBUG_rapend('____vi_rhow')

    return
  end subroutine vi_rhow

end module mod_vi
!-------------------------------------------------------------------------------

