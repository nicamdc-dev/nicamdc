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
  public :: cnvvar_rhogkin
  public :: cnvvar_uv2vh
  public :: cnvvar_vh2uv

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
    use mod_vmtr, only : &
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

    integer :: g, k, l, iq
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

    integer :: g, k, l, iq
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

    integer  :: gall, kmin, kmax

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call PROF_rapstart('CNV_rhogkin',2)

    gall = ADM_gall
    kmin = ADM_kmin
    kmax = ADM_kmax

    do l = 1, ADM_lall
       !$omp parallel default(none),private(g,k),                             &
       !$omp shared(l,gall,kmin,kmax,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
       !$omp        rhogkin_h,rhogkin_v,VMTR_C2Wfact,VMTR_W2Cfact)

       !--- horizontal kinetic energy
       !$omp do
       do k = kmin, kmax
       do g = 1, gall
          rhogkin_h(g,k) = 0.5_RP * ( rhogvx(g,k,l) * rhogvx(g,k,l) &
                                    + rhogvy(g,k,l) * rhogvy(g,k,l) &
                                    + rhogvz(g,k,l) * rhogvz(g,k,l) ) / rhog(g,k,l)
       enddo
       enddo
       !$omp end do

       !--- vertical kinetic energy
       !$omp do
       do k = kmin+1, kmax
       do g = 1, gall
          rhogkin_v(g,k) = 0.5_RP * ( rhogw(g,k,l) * rhogw(g,k,l) ) &
                         / ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k  ,l) &
                           + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )
       enddo
       enddo
       !$omp end do

       !$omp workshare
       rhogkin_v(:,kmin  ) = 0.0_RP
       rhogkin_v(:,kmax+1) = 0.0_RP
       !$omp end workshare

       !--- total kinetic energy
       !$omp do
       do k = kmin, kmax
       do g = 1, gall
          rhogkin(g,k,l) = rhogkin_h(g,k)                             & ! horizontal
                         + ( VMTR_W2Cfact(g,k,1,l) * rhogkin_v(g,k+1) & ! vertical
                           + VMTR_W2Cfact(g,k,2,l) * rhogkin_v(g,k  ) )
       enddo
       enddo
       !$omp end do

       !$omp workshare
       rhogkin(:,kmin-1,l) = 0.0_RP
       rhogkin(:,kmax+1,l) = 0.0_RP
       !$omp end workshare

       !$omp end parallel
    enddo

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

    return
  end subroutine cnvvar_rhogkin

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
       GRD_LAT, &
       GRD_s,   &
       GRD_s_pl
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
       coslat = cos(GRD_s(g,k0,l,GRD_LAT))

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
          coslat = cos(GRD_s_pl(g,k0,l,GRD_LAT))

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
  subroutine cnvvar_vh2uv( &
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
       GRD_LAT, &
       GRD_s,   &
       GRD_s_pl
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

    real(RP), intent(out) :: u    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: u_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: v    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: v_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    logical,  intent(in), optional :: withcos

    real(RP) :: coslat   (ADM_gall   ,ADM_lall   )
    real(RP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)

    integer  :: g, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    coslat   (:,:) = 1.0_RP
    coslat_pl(:,:) = 1.0_RP

    if ( present(withcos) ) then
       if ( withcos ) then
          coslat(:,:) = cos(GRD_s(:,k0,:,GRD_LAT))
          if ( ADM_have_pl ) then
             coslat_pl(:,:) = cos(GRD_s_pl(:,k0,:,GRD_LAT))
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
  end subroutine cnvvar_vh2uv

end module mod_cnvvar
