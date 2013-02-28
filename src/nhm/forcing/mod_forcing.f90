!-------------------------------------------------------------------------------
module mod_forcing
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module is for the artificial forcing
  !       
  ! 
  !++ Current Corresponding Author : R. Yoshida
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      12-10-11  R.Yoshida: convert from phystep for dry dyn-core experiments
  !
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules 
  !
  !-----------------------------------------------------------------------------
  !
  implicit none
  public :: forcing
  !
contains
  !
  subroutine forcing
    !
    use mod_adm, only :        &
         ADM_KNONE,            &
         ADM_gall_in,          &
         ADM_kall,             &
         ADM_kmin,             &
         ADM_kmax,             &
         ADM_lall,             &
         ADM_l_me,             & ! 2010.5.11. M.Satoh
         ADM_prc_run_master,   &
         ADM_LOG_FID
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
         VMTR_VOLUME
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
    use mod_prgvar, only :     &
         prgvar_get_in,        &
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
    use mod_bsstate, only :   &
         phi
    use mod_bndcnd, only : &
         bndcnd_thermo
    use mod_misc
    use mod_af_driver, only : &
         af_driver
    use mod_history, only: &
         history_in
!A.F.    use mod_limiter, only : & ! 2010.5.11 M.Satoh
!A.F.         trc_limiter
    !
    implicit none
    !
    ! Prognostic variables
    !
    !--- rho X ( G^{1/2} X gamma2 )
    Real(8) :: rhog(ADM_gall_in,ADM_kall,ADM_lall)
    !--- rho X ( G^{1/2} X gamma2 ) X vx
    Real(8) :: rhogvx(ADM_gall_in,ADM_kall,ADM_lall)
    !--- rho X ( G^{1/2} X gamma2 ) X vy
    Real(8) :: rhogvy(ADM_gall_in,ADM_kall,ADM_lall)
    !--- rho X ( G^{1/2} X gamma2 ) X vz
    Real(8) :: rhogvz(ADM_gall_in,ADM_kall,ADM_lall)
    !--- rho X ( G^{1/2} X gamma2 ) X w
    Real(8) :: rhogw(ADM_gall_in,ADM_kall,ADM_lall)
    !--- rho X ( G^{1/2} X gamma2 ) X ein
    Real(8) :: rhoge(ADM_gall_in,ADM_kall,ADM_lall)
    !--- rho X ( G^{1/2} X gamma2 ) X q
    Real(8) :: rhogq(ADM_gall_in,ADM_kall,ADM_lall,TRC_VMAX)
    !--- 2011/03/02 NEC add
    Real(8) :: rhogq_l(ADM_gall_in,ADM_kall,TRC_VMAX)
    !
    ! forcing tendency 
    !
    !--- forcing tendensy of rhog  ( G^{1/2} X gamma2 )
    real(8) :: frhog(ADM_gall_in,ADM_kall,ADM_lall)
    !--- forcing tendensy of rhogvx  ( G^{1/2} X gamma2 )
    real(8) :: frhogvx(ADM_gall_in,ADM_kall,ADM_lall)
    !--- forcing tendensy of rhogvy  ( G^{1/2} X gamma2 )
    real(8) :: frhogvy(ADM_gall_in,ADM_kall,ADM_lall)
    !--- forcing tendensy of rhogvz  ( G^{1/2} X gamma2 )
    real(8) :: frhogvz(ADM_gall_in,ADM_kall,ADM_lall)
    !--- forcing tendensy of rhogw  ( G^{1/2} X gamma2 )
    real(8) :: frhogw(ADM_gall_in,ADM_kall,ADM_lall)
    !--- forcing tendensy of rhoge  ( G^{1/2} X gamma2 )
    real(8) :: frhoge(ADM_gall_in,ADM_kall,ADM_lall)
    !--- tendensy of rhogetot  ( G^{1/2} X gamma2 )
    real(8) :: frhogetot(ADM_gall_in,ADM_kall,ADM_lall)
    !--- forcing tendensy of rhogq  ( G^{1/2} X gamma2 )
    real(8) :: frhogq(ADM_gall_in,ADM_kall,ADM_lall,TRC_VMAX)
    !--- 2011/03/02 NEC add
    real(8) :: frhogq_l(ADM_gall_in,ADM_kall,TRC_VMAX)
    !
    ! dynamics
    !
    !--- density ( physical )
    Real(8) :: rho(ADM_gall_in,ADM_kall,ADM_lall)
    !--- temperature ( physical )
    Real(8) :: tem(ADM_gall_in,ADM_kall,ADM_lall)
    !--- pressure ( physical )
    Real(8) :: pre(ADM_gall_in,ADM_kall,ADM_lall)
    !--- potential temperature ( physical )
    real(8) :: th(ADM_gall_in,ADM_kall,ADM_lall)
    !--- horizontal velocity_x  ( physical )
    Real(8) :: vx(ADM_gall_in,ADM_kall,ADM_lall)
    !--- horizontal velocity_y  ( physical )
    Real(8) :: vy(ADM_gall_in,ADM_kall,ADM_lall)
    !--- horizontal velocity_z  ( physical )
    Real(8) :: vz(ADM_gall_in,ADM_kall,ADM_lall)
    !--- vertical velocity ( physical )
    Real(8) :: w(ADM_gall_in,ADM_kall,ADM_lall)
    !
    ! mixing ratio of water substance ( physical )
    !
    Real(8) :: q(ADM_gall_in,ADM_kall,ADM_lall,TRC_VMAX)
    !--- 2011/03/02 NEC add
    Real(8) :: q_l(ADM_gall_in,ADM_kall,TRC_VMAX)

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
    real(8) :: phi_in(ADM_gall_in,ADM_kall,ADM_lall)
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
    !
    real(8) :: q_flux_sfc(ADM_gall_in,ADM_KNONE,ADM_lall,TRC_VMAX) ! tracer flux from surface
    real(8) :: gindex_sfc(ADM_gall_in,ADM_KNONE,ADM_lall)  ! vegetation index    
    integer :: index_sfc(ADM_gall_in,ADM_KNONE,ADM_lall)   ! 07/10/24 T.Mitsui
    integer :: l
    !
    !
    !============================================================================
    !
    !--- structure   
    !
    !============================================================================
    !
    ! getting physical prognostic and diagnostic variables
    !
    ! non-negative fixer of tracers
    ! boundary condition of atmosphere ( k=1 and k=kall )
    ! af_driver  : additional forcing for specific experiment
    !
    ! setting physical prognostic and diagnostic variables
    !
    !============================================================================
    !--- intitialization of forcing
    frhog     = 0.0D0
    frhogvx   = 0.0D0
    frhogvy   = 0.0D0
    frhogvz   = 0.0D0
    frhogw    = 0.0D0
    frhoge    = 0.0D0
    frhogetot = 0.0D0
    frhogq    = 0.0D0
    frhogq_l  = 0.0d0  ![add] 2012/02/01 T.Seiki
    !
    q_flux_sfc(:,:,:,:)=0.d0     ! 09/04/14 [Add] T.Mitsui 
    !--- clipping
    call GTL_clip_region(VMTR_GSGAM2,gsgam2,1,ADM_kall)
    call GTL_clip_region(VMTR_GSGAM2H,gsgam2h,1,ADM_kall)
    call GTL_clip_region(VMTR_GAM2,gam2,1,ADM_kall)
    call GTL_clip_region(VMTR_GAM2H,gam2h,1,ADM_kall)
    call GTL_clip_region(VMTR_VOLUME,vol,1,ADM_kall)
    call GTL_clip_region(GRD_vz(:,:,:,GRD_Z),z,1,ADM_kall)
    call GTL_clip_region(GRD_vz(:,:,:,GRD_ZH),zh,1,ADM_kall)
    call GTL_clip_region(phi,phi_in,1,ADM_kall)
    !
    call GTL_clip_region_1layer_k(GMTR_P_var(:,:,:,GMTR_P_IX),IX,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer_k(GMTR_P_var(:,:,:,GMTR_P_IY),IY,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer_k(GMTR_P_var(:,:,:,GMTR_P_IZ),IZ,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer_k(GMTR_P_var(:,:,:,GMTR_P_JX),JX,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer_k(GMTR_P_var(:,:,:,GMTR_P_JY),JY,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer_k(GMTR_P_var(:,:,:,GMTR_P_JZ),JZ,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer_k(GRD_zs(:,:,:,GRD_ZSFC),zs,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer_k(GRD_zs(:,:,:,GRD_ZSD),zsd,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer_k(GRD_zs(:,:,:,GRD_VEGINDX),gindex_sfc,ADM_KNONE,ADM_KNONE)
    call GTL_clip_region_1layer(GMTR_lat(:,:),lat)
    call GTL_clip_region_1layer(GMTR_lon(:,:),lon)
    index_sfc(:,:,:) = nint ( gindex_sfc(:,:,:) )
    !
    !--- get the prognostic variables
    call prgvar_get_in( &  
         rhog,          &
         rhogvx,        &
         rhogvy,        &
         rhogvz,        &
         rhogw,         &
         rhoge,         &
         rhogq)
    !==========
    !------ Generation of diagnostic values
    !------ and setting the boudary conditions
    do l=1,ADM_lall
       ADM_l_me = l ! 2010.5.11. M.Satoh
       ! 2011/03/02 NEC
       rhogq_l(:,:,:)=rhogq(:,:,l,:)

       call cnvvar_p2d(    &
            ADM_gall_in,   &  !--- in
            rho(:,:,l),    &  !--- OUT
            pre(:,:,l),    &  !--- OUT
            tem(:,:,l),    &  !--- OUT
            vx(:,:,l),     &  !--- OUT
            vy(:,:,l),     &  !--- OUT
            vz(:,:,l),     &  !--- OUT
            w(:,:,l),      &  !--- OUT
            !q(:,:,l,:),    &  !--- OUT
            q_l(:,:,:),    &  !--- OUT ! 2011/03/02 NEC
            rhog(:,:,l),   &  !--- IN
            rhogvx(:,:,l), &  !--- IN
            rhogvy(:,:,l), &  !--- IN
            rhogvz(:,:,l), &  !--- IN
            rhogw(:,:,l),  &  !--- IN
            rhoge(:,:,l),  &  !--- IN
            !rhogq(:,:,l,:),&  !--- IN
            rhogq_l(:,:,:),&  !--- IN ! 2011/03/02 NEC
            gsgam2(:,:,l), &  !--- IN
            gsgam2h(:,:,l) &  !--- IN
            )

       call bndcnd_thermo( &
            ADM_gall_in,   &  !--- IN : number of horizontal grid
            tem(:,:,l),    &  !--- INOUT : temperature  
            rho(:,:,l),    &  !--- INOUT : density 
            pre(:,:,l),    &  !--- INOUT : pressure
            phi_in(:,:,l) )   !--- IN    : geopotential
       vx(:,ADM_kmax+1,l) = vx(:,ADM_kmax,l)
       vy(:,ADM_kmax+1,l) = vy(:,ADM_kmax,l)
       vz(:,ADM_kmax+1,l) = vz(:,ADM_kmax,l)
       vx(:,ADM_kmin-1,l) = vx(:,ADM_kmin,l)
       vy(:,ADM_kmin-1,l) = vy(:,ADM_kmin,l)
       vz(:,ADM_kmin-1,l) = vz(:,ADM_kmin,l)
       ! 2011/03/02 NEC [add]
       q(:,ADM_kmin:ADM_kmax,l,:)=q_l(:,ADM_kmin:ADM_kmax,:)
       !
       q(:,ADM_kmax+1,l,:) = 0.0D0
       q(:,ADM_kmin-1,l,:) = 0.0D0

       call thrmdyn_th( &
            ADM_gall_in,&  !--- in
            th(:,:,l),  &  !--- OUT
            tem(:,:,l), &  !--- IN
            pre(:,:,l) )   !--- IN

       call thrmdyn_qd(  &
            ADM_gall_in, &
            qd(:,:,l),   &
            q(:,:,l,:)   )

       ! [add] R.Yoshida 20120731 (for PS output)
       call sv_pre_sfc(     & 
           ADM_gall_in,                    & !--- IN
           rho(:,:,l),                     & !--- IN
           pre(:,:,l),                     & !--- IN
           z(:,:,l),                       & !--- IN
           zs(:,l),                        & !--- IN
           rho_sfc(:,ADM_KNONE,l),         & !--- OUT 2011/08/16b M.Satoh [add]
           pre_sfc(:,ADM_KNONE,l),          & !--- OUT
           tem(:,:,l),                     & !--- IN
           qd(:,:,l),                     & !--- IN
           q(:,:,l,I_QV),                     & !--- IN
           !rho_sfc2(:,ADM_KNONE,l),         & !--- OUT 2011/08/16b M.Satoh [add]
           !pre_sfc2(:,ADM_KNONE,l),          & !--- OUT
           tem_sfc (:,ADM_KNONE,l)          & !--- OUT
           !rho_sfc3(:,ADM_KNONE,l),          & !--- OUT
           !pre_sfc3(:,ADM_KNONE,l)          & !--- OUT
          )

    end do ! l-loop

    ADM_l_me = 0 ! reset: 2010.5.11. M.Satoh
    !
    frhog = 0.0D0
    frhogvx = 0.0D0
    frhogvy = 0.0D0
    frhogvz = 0.0D0
    frhogw = 0.0D0
    frhoge = 0.0D0
    frhogetot = 0.0D0
    frhogq = 0.0D0
    !
    do l=1,ADM_lall
       ADM_l_me = l ! 2010.5.11. M.Satoh
       ! 2011/03/02 NEC
       q_l(:,:,:)=q(:,:,l,:)
       !
       !--- tendency of radiation
       !--- additional forcing ( default : NONE )
       !--- a part of porting vars is deleted: 12/10/12 R.Yoshida
       call af_driver(           &
            ADM_gall_in,         & !--- IN
            rho(:,:,l),          &  !--- IN
            pre(:,:,l),          &  !--- IN
            tem(:,:,l),          &  !--- IN
            vx(:,:,l),           &  !--- IN
            vy(:,:,l),           &  !--- IN
            vz(:,:,l),           &  !--- IN
            w(:,:,l),            &  !--- IN
            lat(:,l),            &  !--- IN
            z(:,:,l),            &  !--- IN   ! add 11/08/14 A.Noda
            frhog(:,:,l),        &  !--- INOUT
            frhogvx(:,:,l),      &  !--- INOUT
            frhogvy(:,:,l),      &  !--- INOUT
            frhogvz(:,:,l),      &  !--- INOUT
            frhogw(:,:,l),       &  !--- INOUT
            frhoge(:,:,l),       &  !--- INOUT
            frhogetot(:,:,l),    &  !--- INOUT
            frhogq(:,:,l,:),     &  !--- INOUT
            GSGAM2(:,:,l),       &  !--- in
            GSGAM2H(:,:,l),      &  !--- in
            add_type='APPEND'    &  !--- IN
            )
    end do ! l-loop
    ADM_l_me = 0 ! reset: 2010.5.11. M.Satoh
    
    rhog = rhog + TIME_DTL*frhog
    rhogvx = rhogvx + TIME_DTL*frhogvx
    rhogvy = rhogvy + TIME_DTL*frhogvy
    rhogvz = rhogvz + TIME_DTL*frhogvz
    rhogw = rhogw + TIME_DTL*frhogw
    rhoge = rhoge + TIME_DTL*frhoge
    rhogq = rhogq + TIME_DTL*frhogq
    !
    !
    call prgvar_set_in( &  
         rhog,          &
         rhogvx,        &
         rhogvy,        &
         rhogvz,        &
         rhogw,         &
         rhoge,         &
         rhogq)
    !
    ! sfcvar_set_in re-install for PS: R.Yoshida 20120801
    call sfcvar_set_in(      &
         pre_sfc,            &  !--- IN : surface variable
         vid = I_PRE_SFC     &  !--- IN : variable ID
         )
    !
    !
    return
  end subroutine forcing
  !
  !------------------------------------------------  ! 12/10/12 R.Yoshida Imported from mod_sv_driver
  !-----------------------------------------------------------------------------  ! 2010.5.17 M.Satoh
  subroutine sv_pre_sfc(     &
       ijdim,               & !--- IN : number of horizontal grid
       rho,                 & !--- IN : density
       pre,                 & !--- IN : pressure
       z,                   & !--- IN : height
       z_srf,               & !--- IN : surface height
       rho_srf,             & !--- OUT : density at the surface 11/08/16b M.Satoh
       pre_srf,             & !--- OUT : pressure at the surface
       tem,                 & !--- IN : pressure
       qd,                  & !--- IN : pressure
       qv,                  & !--- IN : pressure
       !rho_srf2,            & !--- OUT : density at the surface 11/08/16b M.Satoh
       !pre_srf2,            & !--- OUT : pressure at the surface
       tem_sfc              & !--- OUT : density at the surface 11/08/16b M.Satoh
       !rho_srf3,            & !--- OUT : density at the surface 11/08/16b M.Satoh
       !pre_srf3             & !--- OUT : pressure at the surface
       )
    !
    use mod_adm, only :  &
         kdim => ADM_kall,    &
         kmin => ADM_kmin
    use mod_cnst, only :  &
       CNST_RAIR, &
         CNST_EGRAV
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    implicit none
    !
    integer, intent(in) :: ijdim
    real(8), intent(in) :: rho(1:ijdim,1:kdim)
    real(8), intent(in) :: pre(1:ijdim,1:kdim)
    real(8), intent(in) :: z(1:ijdim,1:kdim)
    real(8), intent(in) :: z_srf(1:ijdim)
    !
    real(8), intent(out) :: rho_srf(1:ijdim)
    real(8), intent(out) :: pre_srf(1:ijdim)

    real(8), intent(in) :: tem(1:ijdim,1:kdim)
    real(8), intent(in) :: qd (1:ijdim,1:kdim)
    real(8), intent(in) :: qv (1:ijdim,1:kdim)

    !real(8), intent(out) :: rho_srf2(1:ijdim)
    !real(8), intent(out) :: pre_srf2(1:ijdim)

    real(8), intent(out) :: tem_sfc(1:ijdim)
    !real(8), intent(out) :: rho_srf3(1:ijdim)
    !real(8), intent(out) :: pre_srf3(1:ijdim)

    real(8) :: LASPdry = 9.771958146487295D-3
    real(8) :: preh, rhoh, temh, zh
    real(8) :: logrho, logrho_sfc

    real(8) :: zz,z1,z2,z3,p1,p2,p3
    real(8) :: lag_intpl
    lag_intpl(zz,z1,p1,z2,p2,z3,p3)              &
         = ((zz-z2)*(zz-z3))/((z1-z2)*(z1-z3))*p1&
         + ((zz-z1)*(zz-z3))/((z2-z1)*(z2-z3))*p2&
         + ((zz-z1)*(zz-z2))/((z3-z1)*(z3-z2))*p3


    integer ::  ij
    real(8) :: pre_s, pre_sfc, f, df
    integer :: ite, itelim = 100
    real(8) :: criteria = 0.1D0

    !--- surface density ( extrapolation )
    do ij = 1, ijdim
       rho_srf(ij) = lag_intpl(z_srf(ij),       &
            z(ij,kmin  ),rho(ij,kmin  ),&
            z(ij,kmin+1),rho(ij,kmin+1),&
            z(ij,kmin+2),rho(ij,kmin+2) &
            )
       rho_srf(ij) = rho(ij,kmin) - ( rho(ij,kmin+1)-rho(ij,kmin) ) / ( z(ij,kmin+1)-z(ij,kmin) ) * ( z(ij,kmin)-z_srf(ij) )
    end do

    !--- surface pressure ( hydrostatic balance )
    do ij = 1, ijdim
       pre_srf(ij) = pre(ij,kmin) + 0.5D0 * ( rho_srf(ij)+rho(ij,kmin) ) * CNST_EGRAV * ( z(ij,kmin)-z_srf(ij) )
    enddo

    do ij = 1, ijdim
       tem_sfc(ij)  = tem(ij,kmin) + LASPdry * ( z(ij,kmin)-z_srf(ij) )

       !preh = 0.5D0 * ( GRD_afac(kmin+1) * pre(ij,kmin+1) + GRD_bfac(kmin+1) * pre(ij,kmin) )
       !temh = 0.5D0 * ( GRD_afac(kmin+1) * tem(ij,kmin+1) + GRD_bfac(kmin+1) * tem(ij,kmin) )
       !rhoh = 0.5D0 * ( GRD_afac(kmin+1) * rho(ij,kmin+1) + GRD_bfac(kmin+1) * rho(ij,kmin) )
       !zh   = 0.5D0 * ( GRD_afac(kmin+1) * z  (ij,kmin+1) + GRD_bfac(kmin+1) * z  (ij,kmin) )

       !pre_srf2(ij) = preh + CNST_EGRAV * ( zh - z_srf(ij) ) * rho(ij,kmin)
       !rho_srf2(ij) = pre_srf2(ij) / ( tem_sfc(ij) * CNST_RAIR )

       !logrho_sfc   = log( rhoh ) + ( CNST_EGRAV / CNST_RAIR * ( zh-z_srf(ij) ) + ( temh-tem_sfc(ij) ) ) / tem(ij,kmin)
       !rho_srf3(ij) = exp(logrho_sfc)
       !pre_srf3(ij) = rho_srf3(ij) * ( tem_sfc(ij) * CNST_RAIR )
    enddo

    !do ij = 1, ijdim
    !   pre_s   = 0.D0
    !   pre_sfc = pre(ij,kmin) ! first guess

    !   ! Newton-Lapson
    !   do ite = 1, itelim
    !      if( abs(pre_sfc-pre_s) <= criteria ) exit

    !      pre_s = pre_sfc

    !      f  = ( log( pre(ij,kmin) ) - log( pre_sfc ) ) / ( z(ij,kmin)-z_srf(ij) ) &
    !         + CNST_EGRAV / CNST_RAIR * 2.D0 / ( tem(ij,kmin) + tem_sfc(ij) )

    !      df = 1.D0 / ( z(ij,kmin)-z_srf(ij) ) / pre_sfc

    !      pre_sfc = pre_s + f/df
    !   enddo

    !   if ( ite > itelim ) then
    !      write(*,*) 'xxx iteration not converged!', pre_sfc-pre_s, pre_sfc, tem(ij,kmin), tem_sfc(ij)
    !   endif

    !   pre_srf3(ij) = pre_sfc
    !enddo

    return
  end subroutine sv_pre_sfc

end module mod_forcing
