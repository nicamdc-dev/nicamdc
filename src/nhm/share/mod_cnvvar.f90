!-------------------------------------------------------------------------------
!> Module variable conversion
!!
!! @par Description
!!         Conversion tools for prognostic variables
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_cnvvar
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof

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
  !++ Public procedures
  !
  public :: cnvvar_prg2diag
  public :: cnvvar_diag2prg
  public :: cnvvar_prg2diag_in
  public :: cnvvar_diag2prg_in
  public :: cnvvar_rhogkin
  public :: cnvvar_rhogkin_in
  public :: cnvvar_uv2vh
  public :: cnvvar_vh2uv
  public :: cnvvar_vh2uv_in
  public :: cnvvar_vh2uv_2D

  interface cnvvar_prg2diag_in
     module procedure cnvvar_prg2diag_in_SP
     module procedure cnvvar_prg2diag_in_DP
  end interface cnvvar_prg2diag_in

  interface cnvvar_diag2prg_in
     module procedure cnvvar_diag2prg_in_SP
     module procedure cnvvar_diag2prg_in_DP
  end interface cnvvar_diag2prg_in

  interface cnvvar_rhogkin_in
     module procedure cnvvar_rhogkin_in_SP
     module procedure cnvvar_rhogkin_in_DP
  end interface cnvvar_rhogkin_in

  interface cnvvar_vh2uv
     module procedure cnvvar_vh2uv_SP
     module procedure cnvvar_vh2uv_DP
  end interface cnvvar_vh2uv

  interface cnvvar_vh2uv_in
     module procedure cnvvar_vh2uv_in_SP
     module procedure cnvvar_vh2uv_in_DP
  end interface cnvvar_vh2uv_in

  interface cnvvar_vh2uv_2D
     module procedure cnvvar_vh2uv_2D_SP
     module procedure cnvvar_vh2uv_2D_DP
  end interface cnvvar_vh2uv_2D

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
  subroutine cnvvar_prg2diag(&
       prg,  prg_pl, &
       diag, diag_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_vmtr, only: &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    use mod_runconf, only: &
       PRG_vmax,  &
       DIAG_vmax, &
       TRC_vmax
    use mod_thrmdyn, only: &
       THRMDYN_tempre
    implicit none

    real(RP), intent(in)  :: prg    (ADM_gall   ,ADM_kall,ADM_lall   ,PRG_vmax )
    real(RP), intent(in)  :: prg_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax )
    real(RP), intent(out) :: diag   (ADM_gall   ,ADM_kall,ADM_lall   ,DIAG_vmax)
    real(RP), intent(out) :: diag_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,DIAG_vmax)

    real(RP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ein      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ein_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhog_h   (ADM_gall,   ADM_kall)
    real(RP) :: rhog_h_pl(ADM_gall_pl,ADM_kall)

    integer  :: g, k, l, iq
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       rho (g,k,l)      = prg(g,k,l,I_RHOG  ) * VMTR_RGSGAM2(g,k,l)
       diag(g,k,l,I_vx) = prg(g,k,l,I_RHOGVX) / prg(g,k,l,I_RHOG)
       diag(g,k,l,I_vy) = prg(g,k,l,I_RHOGVY) / prg(g,k,l,I_RHOG)
       diag(g,k,l,I_vz) = prg(g,k,l,I_RHOGVZ) / prg(g,k,l,I_RHOG)
       ein (g,k,l)      = prg(g,k,l,I_RHOGE ) / prg(g,k,l,I_RHOG)
    enddo
    enddo
    enddo

    do iq = 1, TRC_vmax
    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do g  = 1, ADM_gall
       diag(g,k,l,DIAG_vmax0+iq) = prg(g,k,l,PRG_vmax0+iq) / prg(g,k,l,I_RHOG)
    enddo
    enddo
    enddo
    enddo

    call THRMDYN_tempre( ADM_gall,                  & ! [IN]
                         ADM_kall,                  & ! [IN]
                         ADM_lall,                  & ! [IN]
                         ein (:,:,:),               & ! [IN]
                         rho (:,:,:),               & ! [IN]
                         diag(:,:,:,I_qstr:I_qend), & ! [IN]
                         diag(:,:,:,I_tem),         & ! [OUT]
                         diag(:,:,:,I_pre)          ) ! [OUT]

    do l = 1, ADM_lall
       !------ interpolation of rhog_h
       do k = 2, ADM_kall
       do g = 1, ADM_gall
          rhog_h(g,k) = ( VMTR_C2Wfact(g,k,1,l) * prg(g,k  ,l,I_RHOG) &
                        + VMTR_C2Wfact(g,k,2,l) * prg(g,k-1,l,I_RHOG) )
       enddo
       enddo
       do g = 1, ADM_gall
          rhog_h(g,1) = rhog_h(g,2)
       enddo

       do k = 1, ADM_kall
       do g = 1, ADM_gall
          diag(g,k,l,I_w) = prg(g,k,l,I_RHOGW) / rhog_h(g,k)
       enddo
       enddo
    enddo

    if ( ADM_have_pl ) then

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          rho_pl (g,k,l)      = prg_pl(g,k,l,I_RHOG  ) * VMTR_RGSGAM2_pl(g,k,l)
          diag_pl(g,k,l,I_vx) = prg_pl(g,k,l,I_RHOGVX) / prg_pl(g,k,l,I_RHOG)
          diag_pl(g,k,l,I_vy) = prg_pl(g,k,l,I_RHOGVY) / prg_pl(g,k,l,I_RHOG)
          diag_pl(g,k,l,I_vz) = prg_pl(g,k,l,I_RHOGVZ) / prg_pl(g,k,l,I_RHOG)
          ein_pl (g,k,l)      = prg_pl(g,k,l,I_RHOGE ) / prg_pl(g,k,l,I_RHOG)
       enddo
       enddo
       enddo

       do iq = 1, TRC_vmax
       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do g  = 1, ADM_gall_pl
          diag_pl(g,k,l,DIAG_vmax0+iq) = prg_pl(g,k,l,PRG_vmax0+iq) / prg_pl(g,k,l,I_RHOG)
       enddo
       enddo
       enddo
       enddo

       call THRMDYN_tempre( ADM_gall_pl,                  & ! [IN]
                            ADM_kall,                     & ! [IN]
                            ADM_lall_pl,                  & ! [IN]
                            ein_pl (:,:,:),               & ! [IN]
                            rho_pl (:,:,:),               & ! [IN]
                            diag_pl(:,:,:,I_qstr:I_qend), & ! [IN]
                            diag_pl(:,:,:,I_tem),         & ! [OUT]
                            diag_pl(:,:,:,I_pre)          ) ! [OUT]

       do l = 1, ADM_lall_pl
          !------ interpolation of rhog_h
          do k = 2, ADM_kall
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,k) = ( VMTR_C2Wfact_pl(g,k,1,l) * prg_pl(g,k  ,l,I_RHOG) &
                              + VMTR_C2Wfact_pl(g,k,2,l) * prg_pl(g,k-1,l,I_RHOG) )
          enddo
          enddo
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,1) = rhog_h_pl(g,2)
          enddo

          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             diag_pl(g,k,l,I_w) = prg_pl(g,k,l,I_RHOGW) / rhog_h_pl(g,k)
          enddo
          enddo
       enddo

    endif

    return
  end subroutine cnvvar_prg2diag

  !-----------------------------------------------------------------------------
  subroutine cnvvar_diag2prg( &
       prg,  prg_pl, &
       diag, diag_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_vmtr, only: &
       VMTR_GSGAM2,     &
       VMTR_GSGAM2_pl,  &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    use mod_runconf, only: &
       PRG_vmax,  &
       DIAG_vmax, &
       TRC_vmax
    use mod_thrmdyn, only: &
       THRMDYN_rhoein
    implicit none

    real(RP), intent(out) :: prg    (ADM_gall   ,ADM_kall,ADM_lall   ,PRG_vmax )
    real(RP), intent(out) :: prg_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax )
    real(RP), intent(in)  :: diag   (ADM_gall   ,ADM_kall,ADM_lall   ,DIAG_vmax)
    real(RP), intent(in)  :: diag_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,DIAG_vmax)

    real(RP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ein      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ein_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhog_h   (ADM_gall,   ADM_kall)
    real(RP) :: rhog_h_pl(ADM_gall_pl,ADM_kall)

    integer  :: g, k, l, iq
    !---------------------------------------------------------------------------

    call THRMDYN_rhoein( ADM_gall,                  & ! [IN]
                         ADM_kall,                  & ! [IN]
                         ADM_lall,                  & ! [IN]
                         diag(:,:,:,I_tem),         & ! [IN]
                         diag(:,:,:,I_pre),         & ! [IN]
                         diag(:,:,:,I_qstr:I_qend), & ! [IN]
                         rho (:,:,:),               & ! [OUT]
                         ein (:,:,:)                ) ! [OUT]

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       prg(g,k,l,I_RHOG  ) = rho(g,k,l) * VMTR_GSGAM2(g,k,l)
       prg(g,k,l,I_RHOGVX) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vx)
       prg(g,k,l,I_RHOGVY) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vy)
       prg(g,k,l,I_RHOGVZ) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vz)
       prg(g,k,l,I_RHOGE ) = prg(g,k,l,I_RHOG) * ein (g,k,l)
    enddo
    enddo
    enddo

    do iq = 1, TRC_vmax
    do l  = 1, ADM_lall
    do k  = 1, ADM_kall
    do g  = 1, ADM_gall
       prg(g,k,l,PRG_vmax0+iq) = prg(g,k,l,I_RHOG) * diag(g,k,l,DIAG_vmax0+iq)
    enddo
    enddo
    enddo
    enddo

    do l = 1, ADM_lall
       !------ interpolation of rhog_h
       do k = 2, ADM_kall
       do g = 1, ADM_gall
          rhog_h(g,k) = ( VMTR_C2Wfact(g,k,1,l) * prg(g,k  ,l,I_RHOG) &
                        + VMTR_C2Wfact(g,k,2,l) * prg(g,k-1,l,I_RHOG) )
       enddo
       enddo
       do g = 1, ADM_gall
          rhog_h(g,1) = rhog_h(g,2)
       enddo

       do k = 1, ADM_kall
       do g = 1, ADM_gall
          prg(g,k,l,I_RHOGW) = rhog_h(g,k) * diag(g,k,l,I_w)
       enddo
       enddo
    enddo

    if ( ADM_have_pl ) then

       call THRMDYN_rhoein( ADM_gall_pl,                  & ! [IN]
                            ADM_kall,                     & ! [IN]
                            ADM_lall_pl,                  & ! [IN]
                            diag_pl(:,:,:,I_tem),         & ! [IN]
                            diag_pl(:,:,:,I_pre),         & ! [IN]
                            diag_pl(:,:,:,I_qstr:I_qend), & ! [IN]
                            rho_pl (:,:,:),               & ! [OUT]
                            ein_pl (:,:,:)                ) ! [OUT]

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          prg_pl(g,k,l,I_RHOG  ) = rho_pl(g,k,l) * VMTR_GSGAM2_pl(g,k,l)
          prg_pl(g,k,l,I_RHOGVX) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vx)
          prg_pl(g,k,l,I_RHOGVY) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vy)
          prg_pl(g,k,l,I_RHOGVZ) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vz)
          prg_pl(g,k,l,I_RHOGE ) = prg_pl(g,k,l,I_RHOG) * ein_pl (g,k,l)
       enddo
       enddo
       enddo

       do iq = 1, TRC_vmax
       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do g  = 1, ADM_gall_pl
          prg_pl(g,k,l,PRG_vmax0+iq) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,DIAG_vmax0+iq)
       enddo
       enddo
       enddo
       enddo

       do l = 1, ADM_lall_pl
          !------ interpolation of rhog_h
          do k = 2, ADM_kall
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,k) = ( VMTR_C2Wfact_pl(g,k,1,l) * prg_pl(g,k  ,l,I_RHOG) &
                              + VMTR_C2Wfact_pl(g,k,2,l) * prg_pl(g,k-1,l,I_RHOG) )
          enddo
          enddo
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,1) = rhog_h_pl(g,2)
          enddo

          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             prg_pl(g,k,l,I_RHOGW) = rhog_h_pl(g,k) * diag_pl(g,k,l,I_w)
          enddo
          enddo
       enddo

    endif

    return
  end subroutine cnvvar_diag2prg

  !-----------------------------------------------------------------------------
  subroutine cnvvar_prg2diag_in_SP(&
       ijdim,       &
       rho,         &
       pre,         &
       tem,         &
       vx,          &
       vy,          &
       vz,          &
       w,           &
       q,           &
       rhog,        &
       rhogvx,      &
       rhogvy,      &
       rhogvz,      &
       rhogw,       &
       rhoge,       &
       rhogq,       &
       phi,         &
       gsqgam2,     &
       gsqgam2h,    &
       set_boundary )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       GRD_afact, &
       GRD_bfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_thrmdyn, only: &
       THRMDYN_tempre
    use mod_bndcnd, only: &
       BNDCND_thermo
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(out) :: rho     (ijdim,kdim)
    real(SP), intent(out) :: pre     (ijdim,kdim)
    real(SP), intent(out) :: tem     (ijdim,kdim)
    real(SP), intent(out) :: vx      (ijdim,kdim)
    real(SP), intent(out) :: vy      (ijdim,kdim)
    real(SP), intent(out) :: vz      (ijdim,kdim)
    real(SP), intent(out) :: w       (ijdim,kdim)
    real(SP), intent(out) :: q       (ijdim,kdim,nqmax)
    real(SP), intent(in)  :: rhog    (ijdim,kdim)
    real(SP), intent(in)  :: rhogvx  (ijdim,kdim)
    real(SP), intent(in)  :: rhogvy  (ijdim,kdim)
    real(SP), intent(in)  :: rhogvz  (ijdim,kdim)
    real(SP), intent(in)  :: rhogw   (ijdim,kdim)
    real(SP), intent(in)  :: rhoge   (ijdim,kdim)
    real(SP), intent(in)  :: rhogq   (ijdim,kdim,nqmax)
    real(SP), intent(in)  :: phi     (ijdim,kdim)
    real(SP), intent(in)  :: gsqgam2 (ijdim,kdim)
    real(SP), intent(in)  :: gsqgam2h(ijdim,kdim)
    logical,  intent(in)  :: set_boundary

    real(SP) :: ein(ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = rhog  (ij,k) / gsqgam2(ij,k)
       vx (ij,k) = rhogvx(ij,k) / rhog(ij,k)
       vy (ij,k) = rhogvy(ij,k) / rhog(ij,k)
       vz (ij,k) = rhogvz(ij,k) / rhog(ij,k)
       ein(ij,k) = rhoge (ij,k) / rhog(ij,k)
    enddo
    enddo

!OCL XFILL
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo

    call THRMDYN_tempre( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         ein(:,:),   & ! [IN]
                         rho(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         tem(:,:),   & ! [OUT]
                         pre(:,:)    ) ! [OUT]

    k = kmin-1
    do ij = 1, ijdim
       w(ij,k) = rhogw(ij,k) / ( ( GRD_afact(kmin) * rho(ij,k+1) &
                                 + GRD_bfact(kmin) * rho(ij,k  ) ) * gsqgam2h(ij,k+1) )
    enddo

    do k  = kmin, kmax+1
    do ij = 1, ijdim
       w(ij,k) = rhogw(ij,k) / ( ( GRD_afact(k) * rho(ij,k  ) &
                                 + GRD_bfact(k) * rho(ij,k-1) ) * gsqgam2h(ij,k) )
    enddo
    enddo

    if ( set_boundary ) then
       call BNDCND_thermo( ijdim, & ! [IN]
                           tem,   & ! [INOUT]
                           rho,   & ! [INOUT]
                           pre,   & ! [INOUT]
                           phi    ) ! [IN]

       vx(:,kmax+1) = vx(:,kmax)
       vy(:,kmax+1) = vy(:,kmax)
       vz(:,kmax+1) = vz(:,kmax)
       vx(:,kmin-1) = vx(:,kmin)
       vy(:,kmin-1) = vy(:,kmin)
       vz(:,kmin-1) = vz(:,kmin)

       do nq = 1, nqmax
!OCL XFILL
          q(:,kmax+1,nq) = 0.0_SP
!OCL XFILL
          q(:,kmin-1,nq) = 0.0_SP
       enddo
    endif

    return
  end subroutine cnvvar_prg2diag_in_SP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_prg2diag_in_DP(&
       ijdim,       &
       rho,         &
       pre,         &
       tem,         &
       vx,          &
       vy,          &
       vz,          &
       w,           &
       q,           &
       rhog,        &
       rhogvx,      &
       rhogvy,      &
       rhogvz,      &
       rhogw,       &
       rhoge,       &
       rhogq,       &
       phi,         &
       gsqgam2,     &
       gsqgam2h,    &
       set_boundary )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       GRD_afact, &
       GRD_bfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_thrmdyn, only: &
       THRMDYN_tempre
    use mod_bndcnd, only: &
       BNDCND_thermo
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(out) :: rho     (ijdim,kdim)
    real(DP), intent(out) :: pre     (ijdim,kdim)
    real(DP), intent(out) :: tem     (ijdim,kdim)
    real(DP), intent(out) :: vx      (ijdim,kdim)
    real(DP), intent(out) :: vy      (ijdim,kdim)
    real(DP), intent(out) :: vz      (ijdim,kdim)
    real(DP), intent(out) :: w       (ijdim,kdim)
    real(DP), intent(out) :: q       (ijdim,kdim,nqmax)
    real(DP), intent(in)  :: rhog    (ijdim,kdim)
    real(DP), intent(in)  :: rhogvx  (ijdim,kdim)
    real(DP), intent(in)  :: rhogvy  (ijdim,kdim)
    real(DP), intent(in)  :: rhogvz  (ijdim,kdim)
    real(DP), intent(in)  :: rhogw   (ijdim,kdim)
    real(DP), intent(in)  :: rhoge   (ijdim,kdim)
    real(DP), intent(in)  :: rhogq   (ijdim,kdim,nqmax)
    real(DP), intent(in)  :: phi     (ijdim,kdim)
    real(DP), intent(in)  :: gsqgam2 (ijdim,kdim)
    real(DP), intent(in)  :: gsqgam2h(ijdim,kdim)
    logical,  intent(in)  :: set_boundary

    real(DP) :: ein(ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = rhog  (ij,k) / gsqgam2(ij,k)
       vx (ij,k) = rhogvx(ij,k) / rhog(ij,k)
       vy (ij,k) = rhogvy(ij,k) / rhog(ij,k)
       vz (ij,k) = rhogvz(ij,k) / rhog(ij,k)
       ein(ij,k) = rhoge (ij,k) / rhog(ij,k)
    enddo
    enddo

!OCL XFILL
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo

    call THRMDYN_tempre( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         ein(:,:),   & ! [IN]
                         rho(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         tem(:,:),   & ! [OUT]
                         pre(:,:)    ) ! [OUT]

    k = kmin-1
    do ij = 1, ijdim
       w(ij,k) = rhogw(ij,k) / ( ( GRD_afact(kmin) * rho(ij,k+1) &
                                 + GRD_bfact(kmin) * rho(ij,k  ) ) * gsqgam2h(ij,k+1) )
    enddo

    do k  = kmin, kmax+1
    do ij = 1, ijdim
       w(ij,k) = rhogw(ij,k) / ( ( GRD_afact(k) * rho(ij,k  ) &
                                 + GRD_bfact(k) * rho(ij,k-1) ) * gsqgam2h(ij,k) )
    enddo
    enddo

    if ( set_boundary ) then
       call BNDCND_thermo( ijdim, & ! [IN]
                           tem,   & ! [INOUT]
                           rho,   & ! [INOUT]
                           pre,   & ! [INOUT]
                           phi    ) ! [IN]

       vx(:,kmax+1) = vx(:,kmax)
       vy(:,kmax+1) = vy(:,kmax)
       vz(:,kmax+1) = vz(:,kmax)
       vx(:,kmin-1) = vx(:,kmin)
       vy(:,kmin-1) = vy(:,kmin)
       vz(:,kmin-1) = vz(:,kmin)

       do nq = 1, nqmax
!OCL XFILL
          q(:,kmax+1,nq) = 0.0_DP
!OCL XFILL
          q(:,kmin-1,nq) = 0.0_DP
       enddo
    endif

    return
  end subroutine cnvvar_prg2diag_in_DP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_diag2prg_in_SP(&
       ijdim,       &
       pre,         &
       tem,         &
       vx,          &
       vy,          &
       vz,          &
       w,           &
       q,           &
       rhog,        &
       rhogvx,      &
       rhogvy,      &
       rhogvz,      &
       rhogw,       &
       rhoge,       &
       rhogq,       &
       gsqgam2,     &
       gsqgam2h     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       GRD_afact, &
       GRD_bfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_thrmdyn, only: &
       THRMDYN_rhoein
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: pre     (ijdim,kdim)
    real(SP), intent(in)  :: tem     (ijdim,kdim)
    real(SP), intent(in)  :: vx      (ijdim,kdim)
    real(SP), intent(in)  :: vy      (ijdim,kdim)
    real(SP), intent(in)  :: vz      (ijdim,kdim)
    real(SP), intent(in)  :: w       (ijdim,kdim)
    real(SP), intent(in)  :: q       (ijdim,kdim,nqmax)
    real(SP), intent(out) :: rhog    (ijdim,kdim)
    real(SP), intent(out) :: rhogvx  (ijdim,kdim)
    real(SP), intent(out) :: rhogvy  (ijdim,kdim)
    real(SP), intent(out) :: rhogvz  (ijdim,kdim)
    real(SP), intent(out) :: rhogw   (ijdim,kdim)
    real(SP), intent(out) :: rhoge   (ijdim,kdim)
    real(SP), intent(out) :: rhogq   (ijdim,kdim,nqmax)
    real(SP), intent(in)  :: gsqgam2 (ijdim,kdim)
    real(SP), intent(in)  :: gsqgam2h(ijdim,kdim)

    real(SP) :: rho(ijdim,kdim)
    real(SP) :: ein(ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    call THRMDYN_rhoein( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         tem(:,:),   & ! [IN]
                         pre(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         rho(:,:),   & ! [OUT]
                         ein(:,:)    ) ! [OUT]

    do k  = 1, kdim
    do ij = 1, ijdim
       rhog  (ij,k) = rho(ij,k) * gsqgam2(ij,k)
       rhogvx(ij,k) = vx (ij,k) * rhog(ij,k)
       rhogvy(ij,k) = vy (ij,k) * rhog(ij,k)
       rhogvz(ij,k) = vz (ij,k) * rhog(ij,k)
       rhoge (ij,k) = ein(ij,k) * rhog(ij,k)
    enddo
    enddo

!OCL XFILL
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       rhogq(ij,k,nq) = q(ij,k,nq) * rhog(ij,k)
    enddo
    enddo
    enddo

    k = kmin-1
    do ij = 1, ijdim
       rhogw(ij,k) = w(ij,k) * ( ( GRD_afact(kmin) * rho(ij,k+1) &
                                 + GRD_bfact(kmin) * rho(ij,k  ) ) * gsqgam2h(ij,k+1) )
    enddo

    do k  = kmin, kmax+1
    do ij = 1, ijdim
       rhogw(ij,k) = w(ij,k) * ( ( GRD_afact(k) * rho(ij,k  ) &
                                 + GRD_bfact(k) * rho(ij,k-1) ) * gsqgam2h(ij,k) )
    enddo
    enddo

    return
  end subroutine cnvvar_diag2prg_in_SP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_diag2prg_in_DP(&
       ijdim,       &
       pre,         &
       tem,         &
       vx,          &
       vy,          &
       vz,          &
       w,           &
       q,           &
       rhog,        &
       rhogvx,      &
       rhogvy,      &
       rhogvz,      &
       rhogw,       &
       rhoge,       &
       rhogq,       &
       gsqgam2,     &
       gsqgam2h     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       GRD_afact, &
       GRD_bfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_thrmdyn, only: &
       THRMDYN_rhoein
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: pre     (ijdim,kdim)
    real(DP), intent(in)  :: tem     (ijdim,kdim)
    real(DP), intent(in)  :: vx      (ijdim,kdim)
    real(DP), intent(in)  :: vy      (ijdim,kdim)
    real(DP), intent(in)  :: vz      (ijdim,kdim)
    real(DP), intent(in)  :: w       (ijdim,kdim)
    real(DP), intent(in)  :: q       (ijdim,kdim,nqmax)
    real(DP), intent(out) :: rhog    (ijdim,kdim)
    real(DP), intent(out) :: rhogvx  (ijdim,kdim)
    real(DP), intent(out) :: rhogvy  (ijdim,kdim)
    real(DP), intent(out) :: rhogvz  (ijdim,kdim)
    real(DP), intent(out) :: rhogw   (ijdim,kdim)
    real(DP), intent(out) :: rhoge   (ijdim,kdim)
    real(DP), intent(out) :: rhogq   (ijdim,kdim,nqmax)
    real(DP), intent(in)  :: gsqgam2 (ijdim,kdim)
    real(DP), intent(in)  :: gsqgam2h(ijdim,kdim)

    real(DP) :: rho(ijdim,kdim)
    real(DP) :: ein(ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    call THRMDYN_rhoein( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         tem(:,:),   & ! [IN]
                         pre(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         rho(:,:),   & ! [OUT]
                         ein(:,:)    ) ! [OUT]

    do k  = 1, kdim
    do ij = 1, ijdim
       rhog  (ij,k) = rho(ij,k) * gsqgam2(ij,k)
       rhogvx(ij,k) = vx (ij,k) * rhog(ij,k)
       rhogvy(ij,k) = vy (ij,k) * rhog(ij,k)
       rhogvz(ij,k) = vz (ij,k) * rhog(ij,k)
       rhoge (ij,k) = ein(ij,k) * rhog(ij,k)
    enddo
    enddo

!OCL XFILL
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       rhogq(ij,k,nq) = q(ij,k,nq) * rhog(ij,k)
    enddo
    enddo
    enddo

    k = kmin-1
    do ij = 1, ijdim
       rhogw(ij,k) = w(ij,k) * ( ( GRD_afact(kmin) * rho(ij,k+1) &
                                 + GRD_bfact(kmin) * rho(ij,k  ) ) * gsqgam2h(ij,k+1) )
    enddo

    do k  = kmin, kmax+1
    do ij = 1, ijdim
       rhogw(ij,k) = w(ij,k) * ( ( GRD_afact(k) * rho(ij,k  ) &
                                 + GRD_bfact(k) * rho(ij,k-1) ) * gsqgam2h(ij,k) )
    enddo
    enddo

    return
  end subroutine cnvvar_diag2prg_in_DP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_rhogkin( &
       rhog,    rhog_pl,   &
       rhogvx,  rhogvx_pl, &
       rhogvy,  rhogvy_pl, &
       rhogvz,  rhogvz_pl, &
       rhogw,   rhogw_pl,  &
       rhogkin, rhogkin_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl, &
       VMTR_W2Cfact,    &
       VMTR_W2Cfact_pl
    implicit none

    real(RP), intent(in)  :: rhog      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 )
    real(RP), intent(in)  :: rhog_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vx
    real(RP), intent(in)  :: rhogvx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vy
    real(RP), intent(in)  :: rhogvy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vz
    real(RP), intent(in)  :: rhogvz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw     (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X w
    real(RP), intent(in)  :: rhogw_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogkin   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X kin
    real(RP), intent(out) :: rhogkin_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhogkin_h   (ADM_gall,   ADM_kall) ! rho X ( G^1/2 X gamma2 ) X kin (horizontal)
    real(RP) :: rhogkin_h_pl(ADM_gall_pl,ADM_kall)
    real(RP) :: rhogkin_v   (ADM_gall,   ADM_kall) ! rho X ( G^1/2 X gamma2 ) X kin (vertical)
    real(RP) :: rhogkin_v_pl(ADM_gall_pl,ADM_kall)

    integer  :: gall, kmin, kmax, lall

    integer  :: g, k, l
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(rhogkin_h,rhogkin_v) &
    !$acc pcopy(rhogkin) &
    !$acc pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw) &
    !$acc pcopyin(VMTR_C2Wfact,VMTR_W2Cfact)

    call PROF_rapstart('CNV_rhogkin',2)

    gall = ADM_gall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    !$acc kernels pcopy(rhogkin) present(rhogkin_h,rhogkin_v) &
    !$acc pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,VMTR_C2Wfact,VMTR_W2Cfact)
    !$omp parallel default(none),private(g,k,l), &
    !$omp shared(gall,kmin,kmax,lall,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
    !$omp        rhogkin_h,rhogkin_v,VMTR_C2Wfact,VMTR_W2Cfact)
    do l = 1, lall
       !--- horizontal kinetic energy
!OCL XFILL
       !$omp do
       do k = kmin, kmax
       do g = 1, gall
          rhogkin_h(g,k) = 0.5_RP * ( rhogvx(g,k,l) * rhogvx(g,k,l) &
                                    + rhogvy(g,k,l) * rhogvy(g,k,l) &
                                    + rhogvz(g,k,l) * rhogvz(g,k,l) ) / rhog(g,k,l)
       enddo
       enddo
       !$omp end do nowait

       !--- vertical kinetic energy
!OCL XFILL
       !$omp do
       do k = kmin+1, kmax
       do g = 1, gall
          rhogkin_v(g,k) = 0.5_RP * ( rhogw(g,k,l) * rhogw(g,k,l) ) &
                         / ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k  ,l) &
                           + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )
       enddo
       enddo
       !$omp end do nowait

!OCL XFILL
       !$omp do
       do g = 1, gall
          rhogkin_v(g,kmin  ) = 0.0_RP
          rhogkin_v(g,kmax+1) = 0.0_RP
       enddo
       !$omp end do

       !--- total kinetic energy
!OCL XFILL
       !$omp do
       do k = kmin, kmax
       do g = 1, gall
          rhogkin(g,k,l) = rhogkin_h(g,k)                             & ! horizontal
                         + ( VMTR_W2Cfact(g,k,1,l) * rhogkin_v(g,k+1) & ! vertical
                           + VMTR_W2Cfact(g,k,2,l) * rhogkin_v(g,k  ) )
       enddo
       enddo
       !$omp end do nowait

!OCL XFILL
       !$omp do
       do g = 1, gall
          rhogkin(g,kmin-1,l) = 0.0_RP
          rhogkin(g,kmax+1,l) = 0.0_RP
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          !--- horizontal kinetic energy
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogkin_h_pl(g,k) = 0.5_RP * ( rhogvx_pl(g,k,l) * rhogvx_pl(g,k,l) &
                                          + rhogvy_pl(g,k,l) * rhogvy_pl(g,k,l) &
                                          + rhogvz_pl(g,k,l) * rhogvz_pl(g,k,l) ) / rhog_pl(g,k,l)
          enddo
          enddo

          !--- vertical kinetic energy
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogkin_v_pl(g,k) = 0.5_RP * ( rhogw_pl(g,k,l) * rhogw_pl(g,k,l) ) &
                               / ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k  ,l) &
                                 + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )
          enddo
          enddo
          rhogkin_v_pl(:,ADM_kmin  ) = 0.0_RP
          rhogkin_v_pl(:,ADM_kmax+1) = 0.0_RP

          !--- total kinetic energy
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogkin_pl(g,k,l) = rhogkin_h_pl(g,k)                                & ! horizontal
                               + ( VMTR_W2Cfact_pl(g,k,1,l) * rhogkin_v_pl(g,k+1) & ! vertical
                                 + VMTR_W2Cfact_pl(g,k,2,l) * rhogkin_v_pl(g,k  ) )
          enddo
          enddo
          rhogkin_pl(:,ADM_kmin-1,l) = 0.0_RP
          rhogkin_pl(:,ADM_kmax+1,l) = 0.0_RP
       enddo
    endif

    call PROF_rapend('CNV_rhogkin',2)

    !$acc end data

    return
  end subroutine cnvvar_rhogkin

  !-----------------------------------------------------------------------------
  subroutine cnvvar_rhogkin_in_SP( &
       ijdim,     &
       kdim,      &
       rhog,      &
       rhogvx,    &
       rhogvy,    &
       rhogvz,    &
       rhogw,     &
       C2Wfact,   &
       W2Cfact,   &
       rhogkin,   &
       rhogkin_h, &
       rhogkin_v  )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: rhog     (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 )
    real(SP), intent(in)  :: rhogvx   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vx
    real(SP), intent(in)  :: rhogvy   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vy
    real(SP), intent(in)  :: rhogvz   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vz
    real(SP), intent(in)  :: rhogw    (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X w
    real(SP), intent(in)  :: C2Wfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(SP), intent(in)  :: W2Cfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(SP), intent(out) :: rhogkin  (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin
    real(SP), intent(out) :: rhogkin_h(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (horizontal)
    real(SP), intent(out) :: rhogkin_v(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (vertical)

    integer  :: gall, kmin, kmax

    integer  :: g, k
    !---------------------------------------------------------------------------

    call PROF_rapstart('CNV_rhogkin',2)

    gall = ijdim
    kmin = 2
    kmax = kdim-1

    !$omp parallel default(none),private(g,k),                           &
    !$omp shared(gall,kmin,kmax,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
    !$omp        rhogkin_h,rhogkin_v,C2Wfact,W2Cfact)

    !--- horizontal kinetic energy
!OCL XFILL
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin_h(g,k) = 0.5_SP * ( rhogvx(g,k) * rhogvx(g,k) &
                                 + rhogvy(g,k) * rhogvy(g,k) &
                                 + rhogvz(g,k) * rhogvz(g,k) ) / rhog(g,k)
    enddo
    enddo
    !$omp end do nowait

!OCL XFILL
    !$omp do
    do g = 1, gall
       rhogkin_h(g,kmin-1) = 0.0_SP
       rhogkin_h(g,kmax+1) = 0.0_SP
    enddo
    !$omp end do nowait

    !--- vertical kinetic energy
    !$omp do
    do k = kmin+1, kmax
    do g = 1, gall
       rhogkin_v(g,k) = 0.5_SP * ( rhogw(g,k) * rhogw(g,k) ) &
                      / ( C2Wfact(g,k,1) * rhog(g,k  ) &
                        + C2Wfact(g,k,2) * rhog(g,k-1) )
    enddo
    enddo
    !$omp end do nowait

!OCL XFILL
    !$omp do
    do g = 1, gall
       rhogkin_v(g,kmin-1) = 0.0_SP
       rhogkin_v(g,kmin  ) = 0.0_SP
       rhogkin_v(g,kmax+1) = 0.0_SP
    enddo
    !$omp end do

    !--- total kinetic energy
!OCL XFILL
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin(g,k) = rhogkin_h(g,k)                      & ! horizontal
                    + ( W2Cfact(g,k,1) * rhogkin_v(g,k+1) & ! vertical
                      + W2Cfact(g,k,2) * rhogkin_v(g,k  ) )
    enddo
    enddo
    !$omp end do nowait

!OCL XFILL
    !$omp do
    do g = 1, gall
       rhogkin(g,kmin-1) = 0.0_SP
       rhogkin(g,kmax+1) = 0.0_SP
    enddo
    !$omp end do

    !$omp end parallel

    call PROF_rapend('CNV_rhogkin',2)

    return
  end subroutine cnvvar_rhogkin_in_SP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_rhogkin_in_DP( &
       ijdim,     &
       kdim,      &
       rhog,      &
       rhogvx,    &
       rhogvy,    &
       rhogvz,    &
       rhogw,     &
       C2Wfact,   &
       W2Cfact,   &
       rhogkin,   &
       rhogkin_h, &
       rhogkin_v  )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: rhog     (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 )
    real(DP), intent(in)  :: rhogvx   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vx
    real(DP), intent(in)  :: rhogvy   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vy
    real(DP), intent(in)  :: rhogvz   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vz
    real(DP), intent(in)  :: rhogw    (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X w
    real(DP), intent(in)  :: C2Wfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(DP), intent(in)  :: W2Cfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(DP), intent(out) :: rhogkin  (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin
    real(DP), intent(out) :: rhogkin_h(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (horizontal)
    real(DP), intent(out) :: rhogkin_v(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (vertical)

    integer  :: gall, kmin, kmax

    integer  :: g, k
    !---------------------------------------------------------------------------

    call PROF_rapstart('CNV_rhogkin',2)

    gall = ijdim
    kmin = 2
    kmax = kdim-1

    !$omp parallel default(none),private(g,k),                           &
    !$omp shared(gall,kmin,kmax,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
    !$omp        rhogkin_h,rhogkin_v,C2Wfact,W2Cfact)

    !--- horizontal kinetic energy
!OCL XFILL
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin_h(g,k) = 0.5_DP * ( rhogvx(g,k) * rhogvx(g,k) &
                                 + rhogvy(g,k) * rhogvy(g,k) &
                                 + rhogvz(g,k) * rhogvz(g,k) ) / rhog(g,k)
    enddo
    enddo
    !$omp end do nowait

!OCL XFILL
    !$omp do
    do g = 1, gall
       rhogkin_h(g,kmin-1) = 0.0_DP
       rhogkin_h(g,kmax+1) = 0.0_DP
    enddo
    !$omp end do nowait

    !--- vertical kinetic energy
    !$omp do
    do k = kmin+1, kmax
    do g = 1, gall
       rhogkin_v(g,k) = 0.5_DP * ( rhogw(g,k) * rhogw(g,k) ) &
                      / ( C2Wfact(g,k,1) * rhog(g,k  ) &
                        + C2Wfact(g,k,2) * rhog(g,k-1) )
    enddo
    enddo
    !$omp end do nowait

!OCL XFILL
    !$omp do
    do g = 1, gall
       rhogkin_v(g,kmin-1) = 0.0_DP
       rhogkin_v(g,kmin  ) = 0.0_DP
       rhogkin_v(g,kmax+1) = 0.0_DP
    enddo
    !$omp end do

    !--- total kinetic energy
!OCL XFILL
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin(g,k) = rhogkin_h(g,k)                      & ! horizontal
                    + ( W2Cfact(g,k,1) * rhogkin_v(g,k+1) & ! vertical
                      + W2Cfact(g,k,2) * rhogkin_v(g,k  ) )
    enddo
    enddo
    !$omp end do nowait

!OCL XFILL
    !$omp do
    do g = 1, gall
       rhogkin(g,kmin-1) = 0.0_DP
       rhogkin(g,kmax+1) = 0.0_DP
    enddo
    !$omp end do

    !$omp end parallel

    call PROF_rapend('CNV_rhogkin',2)

    return
  end subroutine cnvvar_rhogkin_in_DP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_uv2vh( &
       ucos, ucos_pl, &
       vcos, vcos_pl, &
       vx,   vx_pl,   &
       vy,   vy_pl,   &
       vz,   vz_pl    )
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_grd, only: &
       GRD_LAT,   &
       GRD_LAT_pl
    use mod_gmtr, only: &
       GMTR_p_IX, &
       GMTR_p_IY, &
       GMTR_p_IZ, &
       GMTR_p_JX, &
       GMTR_p_JY, &
       GMTR_p_JZ, &
       GMTR_p,    &
       GMTR_p_pl
    implicit none

    real(RP), intent(in)  :: ucos   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: ucos_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vcos   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vcos_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vx     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vy     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vz     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: u, v, coslat, sw

    integer  :: g, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       coslat = cos(GRD_LAT(g,l))

       sw = 0.5_RP + sign(0.5_RP,-abs(coslat)) ! if (coslat == 0), u=v=0

       u = ucos(g,k,l) * ( 1.0_RP - sw ) / ( coslat - sw )
       v = vcos(g,k,l) * ( 1.0_RP - sw ) / ( coslat - sw )

       vx(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IX) &
                 + v * GMTR_p(g,k0,l,GMTR_p_JX)
       vy(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IY) &
                 + v * GMTR_p(g,k0,l,GMTR_p_JY)
       vz(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IZ) &
                 + v * GMTR_p(g,k0,l,GMTR_p_JZ)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          coslat = cos(GRD_LAT_pl(g,l))

          sw = 0.5_RP + sign(0.5_RP,-abs(coslat)) ! if (coslat == 0), u=v=0

          u = ucos_pl(g,k,l) * ( 1.0_RP - sw ) / ( coslat - sw )
          v = vcos_pl(g,k,l) * ( 1.0_RP - sw ) / ( coslat - sw )

          vx_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JX)
          vy_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JY)
          vz_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IZ) &
                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JZ)
       enddo
       enddo
       enddo
    endif

    return
  end subroutine cnvvar_uv2vh

  !-----------------------------------------------------------------------------
  subroutine cnvvar_vh2uv_SP( &
       u,  u_pl,  &
       v,  v_pl,  &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl, &
       withcos    )
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_grd, only: &
       GRD_LAT,   &
       GRD_LAT_pl
    use mod_gmtr, only: &
       GMTR_p_IX, &
       GMTR_p_IY, &
       GMTR_p_IZ, &
       GMTR_p_JX, &
       GMTR_p_JY, &
       GMTR_p_JZ, &
       GMTR_p,    &
       GMTR_p_pl
    implicit none

    real(SP), intent(out) :: u    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(SP), intent(out) :: u_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(SP), intent(out) :: v    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(SP), intent(out) :: v_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(SP), intent(in)  :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(SP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(SP), intent(in)  :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(SP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(SP), intent(in)  :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(SP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    logical,  intent(in), optional :: withcos

    real(SP) :: coslat   (ADM_gall   ,ADM_lall   )
    real(SP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)

    integer  :: g, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    coslat   (:,:) = 1.0_SP
    coslat_pl(:,:) = 1.0_SP

    if ( present(withcos) ) then
       if ( withcos ) then
          coslat(:,:) = cos(GRD_LAT(:,:))
          if ( ADM_have_pl ) then
             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
          endif
       endif
    endif

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       u(g,k,l) = ( vx(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IX) &
                  + vy(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IY) &
                  + vz(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IZ) ) * coslat(g,l)
       v(g,k,l) = ( vx(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JX) &
                  + vy(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JY) &
                  + vz(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JZ) ) * coslat(g,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          u_pl(g,k,l) = ( vx_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
                        + vy_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
                        + vz_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IZ) ) * coslat_pl(g,l)
          v_pl(g,k,l) = ( vx_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JX) &
                        + vy_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JY) &
                        + vz_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JZ) ) * coslat_pl(g,l)
       enddo
       enddo
       enddo
    endif

    return
  end subroutine cnvvar_vh2uv_SP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_vh2uv_DP( &
       u,  u_pl,  &
       v,  v_pl,  &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl, &
       withcos    )
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_grd, only: &
       GRD_LAT,   &
       GRD_LAT_pl
    use mod_gmtr, only: &
       GMTR_p_IX, &
       GMTR_p_IY, &
       GMTR_p_IZ, &
       GMTR_p_JX, &
       GMTR_p_JY, &
       GMTR_p_JZ, &
       GMTR_p,    &
       GMTR_p_pl
    implicit none

    real(DP), intent(out) :: u    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(DP), intent(out) :: u_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(DP), intent(out) :: v    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(DP), intent(out) :: v_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(DP), intent(in)  :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(DP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(DP), intent(in)  :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(DP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(DP), intent(in)  :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(DP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    logical,  intent(in), optional :: withcos

    real(DP) :: coslat   (ADM_gall   ,ADM_lall   )
    real(DP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)

    integer  :: g, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    coslat   (:,:) = 1.0_DP
    coslat_pl(:,:) = 1.0_DP

    if ( present(withcos) ) then
       if ( withcos ) then
          coslat(:,:) = cos(GRD_LAT(:,:))
          if ( ADM_have_pl ) then
             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
          endif
       endif
    endif

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       u(g,k,l) = ( vx(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IX) &
                  + vy(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IY) &
                  + vz(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IZ) ) * coslat(g,l)
       v(g,k,l) = ( vx(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JX) &
                  + vy(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JY) &
                  + vz(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JZ) ) * coslat(g,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          u_pl(g,k,l) = ( vx_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
                        + vy_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
                        + vz_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IZ) ) * coslat_pl(g,l)
          v_pl(g,k,l) = ( vx_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JX) &
                        + vy_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JY) &
                        + vz_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JZ) ) * coslat_pl(g,l)
       enddo
       enddo
       enddo
    endif

    return
  end subroutine cnvvar_vh2uv_DP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_vh2uv_in_SP( &
       ijdim, &
       u,     &
       v,     &
       vx,    &
       vy,    &
       vz,    &
       ix,    &
       iy,    &
       iz,    &
       jx,    &
       jy,    &
       jz     )
    use mod_adm, only: &
       kdim => ADM_kall
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(out) :: u (ijdim,kdim)
    real(SP), intent(out) :: v (ijdim,kdim)
    real(SP), intent(in)  :: vx(ijdim,kdim)
    real(SP), intent(in)  :: vy(ijdim,kdim)
    real(SP), intent(in)  :: vz(ijdim,kdim)
    real(SP), intent(in)  :: ix(ijdim)
    real(SP), intent(in)  :: iy(ijdim)
    real(SP), intent(in)  :: iz(ijdim)
    real(SP), intent(in)  :: jx(ijdim)
    real(SP), intent(in)  :: jy(ijdim)
    real(SP), intent(in)  :: jz(ijdim)

    integer  :: ij, k
    !---------------------------------------------------------------------------

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,u,v,vx,vy,vz,ix,iy,iz,jx,jy,jz)
    do k  = 1, kdim
    do ij = 1, ijdim
       u(ij,k) = vx(ij,k) * ix(ij) &
               + vy(ij,k) * iy(ij) &
               + vz(ij,k) * iz(ij)
       v(ij,k) = vx(ij,k) * jx(ij) &
               + vy(ij,k) * jy(ij) &
               + vz(ij,k) * jz(ij)
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine cnvvar_vh2uv_in_SP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_vh2uv_in_DP( &
       ijdim, &
       u,     &
       v,     &
       vx,    &
       vy,    &
       vz,    &
       ix,    &
       iy,    &
       iz,    &
       jx,    &
       jy,    &
       jz     )
    use mod_adm, only: &
       kdim => ADM_kall
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(out) :: u (ijdim,kdim)
    real(DP), intent(out) :: v (ijdim,kdim)
    real(DP), intent(in)  :: vx(ijdim,kdim)
    real(DP), intent(in)  :: vy(ijdim,kdim)
    real(DP), intent(in)  :: vz(ijdim,kdim)
    real(DP), intent(in)  :: ix(ijdim)
    real(DP), intent(in)  :: iy(ijdim)
    real(DP), intent(in)  :: iz(ijdim)
    real(DP), intent(in)  :: jx(ijdim)
    real(DP), intent(in)  :: jy(ijdim)
    real(DP), intent(in)  :: jz(ijdim)

    integer  :: ij, k
    !---------------------------------------------------------------------------

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,u,v,vx,vy,vz,ix,iy,iz,jx,jy,jz)
    do k  = 1, kdim
    do ij = 1, ijdim
       u(ij,k) = vx(ij,k) * ix(ij) &
               + vy(ij,k) * iy(ij) &
               + vz(ij,k) * iz(ij)
       v(ij,k) = vx(ij,k) * jx(ij) &
               + vy(ij,k) * jy(ij) &
               + vz(ij,k) * jz(ij)
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine cnvvar_vh2uv_in_DP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_vh2uv_2D_SP( &
       u,  u_pl,  &
       v,  v_pl,  &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl, &
       withcos    )
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl
    use mod_grd, only: &
       GRD_LAT,   &
       GRD_LAT_pl
    use mod_gmtr, only: &
       GMTR_p_IX, &
       GMTR_p_IY, &
       GMTR_p_IZ, &
       GMTR_p_JX, &
       GMTR_p_JY, &
       GMTR_p_JZ, &
       GMTR_p,    &
       GMTR_p_pl
    implicit none

    real(SP), intent(out) :: u    (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(SP), intent(out) :: u_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(SP), intent(out) :: v    (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(SP), intent(out) :: v_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(SP), intent(in)  :: vx   (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(SP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(SP), intent(in)  :: vy   (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(SP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(SP), intent(in)  :: vz   (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(SP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)

    logical,  intent(in), optional :: withcos

    real(SP) :: coslat   (ADM_gall   ,ADM_lall   )
    real(SP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)

    integer  :: g, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    coslat   (:,:) = 1.0_SP
    coslat_pl(:,:) = 1.0_SP

    if ( present(withcos) ) then
       if ( withcos ) then
          coslat(:,:) = cos(GRD_LAT(:,:))
          if ( ADM_have_pl ) then
             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
          endif
       endif
    endif

    do l = 1, ADM_lall
    do g = 1, ADM_gall
       u(g,k0,l) = ( vx(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IX) &
                   + vy(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IY) &
                   + vz(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IZ) ) * coslat(g,l)
       v(g,k0,l) = ( vx(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JX) &
                   + vy(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JY) &
                   + vz(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JZ) ) * coslat(g,l)
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl
          u_pl(g,k0,l) = ( vx_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
                         + vy_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
                         + vz_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IZ) ) * coslat_pl(g,l)
          v_pl(g,k0,l) = ( vx_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JX) &
                         + vy_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JY) &
                         + vz_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JZ) ) * coslat_pl(g,l)
       enddo
       enddo
    endif

    return
  end subroutine cnvvar_vh2uv_2D_SP

  !-----------------------------------------------------------------------------
  subroutine cnvvar_vh2uv_2D_DP( &
       u,  u_pl,  &
       v,  v_pl,  &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl, &
       withcos    )
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl
    use mod_grd, only: &
       GRD_LAT,   &
       GRD_LAT_pl
    use mod_gmtr, only: &
       GMTR_p_IX, &
       GMTR_p_IY, &
       GMTR_p_IZ, &
       GMTR_p_JX, &
       GMTR_p_JY, &
       GMTR_p_JZ, &
       GMTR_p,    &
       GMTR_p_pl
    implicit none

    real(DP), intent(out) :: u    (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(DP), intent(out) :: u_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(DP), intent(out) :: v    (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(DP), intent(out) :: v_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(DP), intent(in)  :: vx   (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(DP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(DP), intent(in)  :: vy   (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(DP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(DP), intent(in)  :: vz   (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(DP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)

    logical,  intent(in), optional :: withcos

    real(DP) :: coslat   (ADM_gall   ,ADM_lall   )
    real(DP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)

    integer  :: g, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    coslat   (:,:) = 1.0_DP
    coslat_pl(:,:) = 1.0_DP

    if ( present(withcos) ) then
       if ( withcos ) then
          coslat(:,:) = cos(GRD_LAT(:,:))
          if ( ADM_have_pl ) then
             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
          endif
       endif
    endif

    do l = 1, ADM_lall
    do g = 1, ADM_gall
       u(g,k0,l) = ( vx(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IX) &
                   + vy(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IY) &
                   + vz(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IZ) ) * coslat(g,l)
       v(g,k0,l) = ( vx(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JX) &
                   + vy(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JY) &
                   + vz(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JZ) ) * coslat(g,l)
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl
          u_pl(g,k0,l) = ( vx_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
                         + vy_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
                         + vz_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IZ) ) * coslat_pl(g,l)
          v_pl(g,k0,l) = ( vx_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JX) &
                         + vy_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JY) &
                         + vz_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JZ) ) * coslat_pl(g,l)
       enddo
       enddo
    endif

    return
  end subroutine cnvvar_vh2uv_2D_DP

end module mod_cnvvar
