!-------------------------------------------------------------------------------
!>
!! Module artificial forcing
!!
!! @par Description
!!         This module is for the artificial forcing
!!
!! @author R.Yoshida
!!
!! @par History
!! @li      2012-10-11 (R.Yoshida) [NEW] extract from phystep
!!
!<
module mod_forcing
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: forcing

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine forcing
    use mod_adm, only :        &
         ADM_KNONE,            &
         ADM_gall_in,          &
         ADM_kall,             &
         ADM_kmin,             &
         ADM_kmax,             &
         ADM_lall,             &
         ADM_l_me,             &
         ADM_prc_run_master
    use mod_grd, only :        &
         GRD_vz,               &
         GRD_Z, GRD_ZH,        &
         GRD_ZSFC,             &
         GRD_ZSD,              &
         GRD_zs,               & ! 07/07/24 K.Suzuki add for SPRINTARS
         GRD_VEGINDX            ! 07/07/24 K.Suzuki add for SPRINTARS
    use mod_gmtr, only :       &
         GMTR_P_var,           &
         GMTR_P_IX,            &
         GMTR_P_IY,            &
         GMTR_P_IZ,            &
         GMTR_P_JX,            &
         GMTR_P_JY,            &
         GMTR_P_JZ,            &
         GMTR_lat,             &
         GMTR_lon
    use mod_vmtr, only :  &
         VMTR_GSGAM2,     &
         VMTR_GSGAM2H,    &
         VMTR_RGSGAM2,    &
         VMTR_GAM2,       &
         VMTR_GAM2H,      &
         VMTR_VOLUME, &
         VMTR_PHI
    use mod_time, only :       &
         TIME_DTL,             &
         TIME_CTIME,           & ! 07/07/24 K.Suzuki add for SPRINTARS
         TIME_CSTEP
    use mod_runconf, only :    &
         TRC_VMAX, &
         I_QV
    use mod_cnst, only :     &
         CNST_EGRAV,         &
         CNST_LH00,          &
         CNST_LH0,           &
         CNST_LHS00,         &
         CNST_LHS0,          &
         CNST_CP,            &
         CNST_CL,            &
         CNST_PRE00,         &
         CNST_KAPPA,         &
         CNST_UNDEF            ! [add] 10/11/14 A.Noda
    use mod_cnvvar, only :     &
         cnvvar_d2p,           &
         cnvvar_p2d,           &
         cnvvar_rhokin
    use mod_prgvar, only: &
         prgvar_get_in_withdiag, &
         prgvar_set_in
    use mod_sfcvar, only :    &
         sfcvar_set_in,       &
         I_PRE_SFC
    use mod_thrmdyn, only :    &
         thrmdyn_qd,           &
         thrmdyn_th,           &
         thrmdyn_eth,          &
         thrmdyn_rho,          &
         thrmdyn_tempre
    use mod_gtl, only :         &
         GTL_clip_region,       &
         GTL_clip_region_1layer,&
         GTL_clip_region_1layer_k
    use mod_diagvar, only : &
         diagvar_set_in, &
         diagvar_get_in, &
         diagvar_set_in_1layer, & ! 10/05/22 M.Satoh for Tiedtke
         diagvar_get_in_1layer, & ! 10/05/22 M.Satoh for Tiedtke
         diagvar_get_in_1layer_k  ! 11/03/02 NEC
    use mod_bndcnd, only : &
         bndcnd_thermo
    use mod_misc
    use mod_af_driver, only : &
         af_driver
    use mod_history, only: &
         history_in
    implicit none

    real(8) :: rhog  (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogvx(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogvy(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogvz(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogw (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhoge (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogq (ADM_gall_in,ADM_kall,ADM_lall,TRC_vmax)
    real(8) :: rho   (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: pre   (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: tem   (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: vx    (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: vy    (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: vz    (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: w     (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: q     (ADM_gall_in,ADM_kall,ADM_lall,TRC_vmax)

    !
    ! forcing tendency 
    !
    real(8) :: frhog(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: frhogvx(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: frhogvy(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: frhogvz(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: frhogw(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: frhoge(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: frhogetot(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: frhogq(ADM_gall_in,ADM_kall,ADM_lall,TRC_VMAX)

    real(8) :: th(ADM_gall_in,ADM_kall,ADM_lall)
    Real(8) :: qd(ADM_gall_in,ADM_kall,ADM_lall)
    !
    ! geometry, coordinate
    !
    Real(8) :: gsgam2(ADM_gall_in,ADM_kall,ADM_lall)
    Real(8) :: gsgam2h(ADM_gall_in,ADM_kall,ADM_lall)
    Real(8) :: gam2(ADM_gall_in,ADM_kall,ADM_lall)
    Real(8) :: gam2h(ADM_gall_in,ADM_kall,ADM_lall)
    Real(8) :: z(ADM_gall_in,ADM_kall,ADM_lall)
    Real(8) :: zh(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: vol(ADM_gall_in,ADM_kall,ADM_lall)
    !
    Real(8) :: zs(ADM_gall_in,ADM_lall)
    Real(8) :: zsd(ADM_gall_in,ADM_lall)
    real(8) :: lat(ADM_gall_in,ADM_lall)
    real(8) :: lon(ADM_gall_in,ADM_lall)
    real(8) :: phi(ADM_gall_in,ADM_kall,ADM_lall)
    !
    real(8) :: ix(ADM_gall_in,ADM_lall)
    real(8) :: iy(ADM_gall_in,ADM_lall)
    real(8) :: iz(ADM_gall_in,ADM_lall)
    real(8) :: jx(ADM_gall_in,ADM_lall)
    real(8) :: jy(ADM_gall_in,ADM_lall)
    real(8) :: jz(ADM_gall_in,ADM_lall)
    !
    !--- surface value
    !
    real(8) :: rho_sfc(ADM_gall_in,ADM_KNONE,ADM_lall)
    real(8) :: pre_sfc(ADM_gall_in,ADM_KNONE,ADM_lall)
    real(8) :: rho_sfc2(ADM_gall_in,ADM_KNONE,ADM_lall)
    real(8) :: pre_sfc2(ADM_gall_in,ADM_KNONE,ADM_lall)
    real(8) :: tem_sfc (ADM_gall_in,ADM_KNONE,ADM_lall)
    real(8) :: rho_sfc3(ADM_gall_in,ADM_KNONE,ADM_lall)
    real(8) :: pre_sfc3(ADM_gall_in,ADM_KNONE,ADM_lall)

    integer :: l
    !---------------------------------------------------------------------------

    call GTL_clip_region(VMTR_GSGAM2 (:,:,:),gsgam2, 1,ADM_kall)
    call GTL_clip_region(VMTR_GSGAM2H(:,:,:),gsgam2h,1,ADM_kall)
    call GTL_clip_region(VMTR_PHI    (:,:,:),phi,    1,ADM_kall)
    call GTL_clip_region(GRD_vz(:,:,:,GRD_Z),z,      1,ADM_kall)

    call GTL_clip_region_1layer(GMTR_lat(:,:),lat)

    !--- get the prognostic and diagnostic variables
    call prgvar_get_in_withdiag( rhog,   & ! [IN]
                                 rhogvx, & ! [IN]
                                 rhogvy, & ! [IN]
                                 rhogvz, & ! [IN]
                                 rhogw,  & ! [IN]
                                 rhoge,  & ! [IN]
                                 rhogq,  & ! [IN]
                                 rho,    & ! [IN]
                                 pre,    & ! [IN]
                                 tem,    & ! [IN]
                                 vx,     & ! [IN]
                                 vy,     & ! [IN]
                                 vz,     & ! [IN]
                                 w,      & ! [IN]
                                 q       ) ! [IN]

    !--- get the diagnostic values
    do l = 1, ADM_lall
       call bndcnd_thermo( ADM_gall_in, & ! [IN]
                           tem(:,:,l),  & ! [INOUT]
                           rho(:,:,l),  & ! [INOUT]
                           pre(:,:,l),  & ! [INOUT]
                           phi(:,:,l)   ) ! [IN]

       vx(:,ADM_kmax+1,l) = vx(:,ADM_kmax,l)
       vy(:,ADM_kmax+1,l) = vy(:,ADM_kmax,l)
       vz(:,ADM_kmax+1,l) = vz(:,ADM_kmax,l)
       vx(:,ADM_kmin-1,l) = vx(:,ADM_kmin,l)
       vy(:,ADM_kmin-1,l) = vy(:,ADM_kmin,l)
       vz(:,ADM_kmin-1,l) = vz(:,ADM_kmin,l)

       q(:,ADM_kmax+1,l,:) = 0.D0
       q(:,ADM_kmin-1,l,:) = 0.D0
    enddo

    do l = 1, ADM_lall
       call af_driver( ADM_gall_in,        & ! [IN]
                       rho      (:,:,l),   & ! [IN]
                       pre      (:,:,l),   & ! [IN]
                       tem      (:,:,l),   & ! [IN]
                       vx       (:,:,l),   & ! [IN]
                       vy       (:,:,l),   & ! [IN]
                       vz       (:,:,l),   & ! [IN]
                       w        (:,:,l),   & ! [IN]
                       lat      (:,l),     & ! [IN]
                       z        (:,:,l),   & ! [IN]
                       frhog    (:,:,l),   & ! [INOUT]
                       frhogvx  (:,:,l),   & ! [INOUT]
                       frhogvy  (:,:,l),   & ! [INOUT]
                       frhogvz  (:,:,l),   & ! [INOUT]
                       frhogw   (:,:,l),   & ! [INOUT]
                       frhoge   (:,:,l),   & ! [INOUT]
                       frhogetot(:,:,l),   & ! [INOUT]
                       frhogq   (:,:,l,:), & ! [INOUT]
                       GSGAM2   (:,:,l),   & ! [IN]
                       GSGAM2H  (:,:,l)    ) ! [IN]
    enddo

    rhog  (:,:,:) = rhog  (:,:,:) + TIME_DTL * frhog  (:,:,:)
    rhogvx(:,:,:) = rhogvx(:,:,:) + TIME_DTL * frhogvx(:,:,:)
    rhogvy(:,:,:) = rhogvy(:,:,:) + TIME_DTL * frhogvy(:,:,:)
    rhogvz(:,:,:) = rhogvz(:,:,:) + TIME_DTL * frhogvz(:,:,:)
    rhogw (:,:,:) = rhogw (:,:,:) + TIME_DTL * frhogw (:,:,:)
    rhoge (:,:,:) = rhoge (:,:,:) + TIME_DTL * frhoge (:,:,:)

    rhogq(:,:,:,:) = rhogq(:,:,:,:) + TIME_DTL * frhogq(:,:,:,:)

    call prgvar_set_in( rhog,   & ! [IN]
                        rhogvx, & ! [IN]
                        rhogvy, & ! [IN]
                        rhogvz, & ! [IN]
                        rhogw,  & ! [IN]
                        rhoge,  & ! [IN]
                        rhogq   ) ! [IN]

    return
  end subroutine forcing

end module mod_forcing
!-------------------------------------------------------------------------------
