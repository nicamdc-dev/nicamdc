!-------------------------------------------------------------------------------
!
!+  history variables
!
!-------------------------------------------------------------------------------
module mod_history_vars
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module is for the calculation of variables 
  !       based on mod_history_vars_cfmip which is derived from 
  !       the old mod_history_vars(version 06-09-07) ( W. Yanase )
  !       
  ! 
  !++ Current Corresponding Author : W.Yanase, S.Iga, H.Tomita, Y.Niwa
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !                07-06-27   Y.Niwa: Imported from mod_history_vars_cfmip
  !                                   add new variables
  !                                    u, v, vx, vy, vz, land vars, ocean vars
  !                                   add history_vars_setup for 
  !                                    calculating options.
  !                                   albedo => albedo_sfc, add REAL albedo
  !                07-07-31   K.Suzuki: add SPRINTARS variables
  !                07-08-06   Y.Niwa: bug fix  'th' => 'ml_th' 'rh' => 'ml_rh'
  !                07-08-06   K.Suzuki: move SPRINTARS variables from mod_postphystep
  !                07-08-20   W.Yanase: debug in "call getr1_in_mat"
  !                07-11-22   T.Mitsui: Effective Radius is calculated in rd_driver => mp_driver
  !                                     and some variables are moved from aerovar 
  !                07-11-30   T.Mitsui: del double use(trivial)
  !                07-12-05   T.Mitsui: add radiative flux categolized as ISCCP-D series
  !                08-02-19   T.Mitsui: Add output for slice-cloud normalized by cot 
  !                                     and trivial fix for output option
  !                08-03-06   H.Tomita: add mp_wsm3_wcdiag.
  !                08-05-24   T.Mitsui: add nc, nr, ni, ns, ng for 2moment option
  !                08-06-13   T.Mitsui: change adiabat_diag => adiabat_diag2
  !                                     and add arguments (positive only cape, and all cin)
  !                08-07-30   W.Yanase: add sl_u10m, sl_v10m, sl_tauu, sl_tauv
  !                09-01-23   H.Tomita: add ml_omg, ml_hgt
  !                09-04-14   H.Tomita: bugfix : zero clear of lwpqc,lwpqr.
  !                09-04-14   T.Mitsui: sl_nc,r,i,s,g and sl_tauc,r,i,s,g
  !                                     sl_qi, sl_qg
  !                09-07-10   H.Tomita: 1. Some quantities related with energy conservation
  !                                     were added.
  !                                     2. Some tag names were changed to avoid the confusion.
  !                                       'sl_swi' -> 'sl_swd_sfc'
  !                                       'sl_sswr' -> 'sl_swu_sfc'
  !                                       'sl_slwd' -> 'sl_lwd_sfc'
  !                                       'sl_slwu' -> 'sl_lwu_sfc'
  !                                       'sl_slh' ->  'sl_lh_sfc'
  !                                       'sl_ssh' ->  'sl_sh_sfc'
  !                                       'sl_sw_toai' ->  'sl_swd_toa'
  !                                       'sl_sw_toar' ->  'sl_swu_toa'
  !                                       'sl_lw_toa' ->  'sl_lwu_toa'
  !                09-07-13   S.Iga: Bug fix:  'sl_lw_toa' ->  'sl_lwu_toa'
  !                10-06-19   A.T.Noda:
  !                                Allow to use a convection parameterization
  !                                with an advanced microphysics schemes, 
  !                                such as G98, NSW?,
  !                10-08-20   C.Kodama: land model output is filled with undef values over the ocean.
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_MAXFNAME, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: history_vars_setup
  public :: history_vars

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  logical, private, save :: out_cld_frac = .false.
  logical, private, save :: out_cldw = .false.
  logical, private, save :: out_qr = .false.
  logical, private, save :: out_qs = .false.
  logical, private, save :: out_qi = .true.
  logical, private, save :: out_qg = .true.
  logical, private, save :: out_cldi = .false.
  logical, private, save :: out_slh = .false.
  logical, private, save :: out_tauucos = .false.
  logical, private, save :: out_tauvcos = .false.
  logical, private, save :: out_ucos10m = .false.
  logical, private, save :: out_vcos10m = .false.
  logical, private, save :: out_tauu = .false.
  logical, private, save :: out_tauv = .false.
  logical, private, save :: out_u10m = .false.
  logical, private, save :: out_v10m = .false.
  logical, private, save :: out_vap_atm = .false.
  logical, private, save :: out_tem_atm = .false.
  logical, private, save :: out_albedo = .false.
  logical, private, save :: out_slp = .false.
  logical, private, save :: out_cape_cin = .false.
  logical, private, save :: out_ucos_vcos = .false.
  logical, private, save :: out_u_v = .false.
  logical, private, save :: out_vx_vy_vz = .false.
  logical, private, save :: out_rh = .false.
  logical, private, save :: out_th = .false.
  logical, private, save :: out_omg = .false.
  logical, private, save :: out_hgt = .false.
  logical, private, save :: out_isccp_rad = .false.
  logical, private, save :: out_landvars = .false.
  logical, private, save :: out_oceanvars = .false.
  logical, private, save :: out_tauclk = .false.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine history_vars_setup
    use mod_history, only: &
       HIST_req_nmax, &
       item_save
!    use mod_runconf, only : &
!         LAND_TYPE,  &
!         OCEAN_TYPE
!    use mod_landvar_bc, only : &
!         GLVNAME_bc => GLVNAME, &
!         I_ALL_bc => I_ALL
!    use mod_landvar_matsiro, only : &
!         GLVNAME_mat => GLVNAME, &
!         I_ALL_mat => I_ALL
!    use mod_oceanvar_mixedlayer, only : &
!         GOVNAME,  &
!         I_ALL_ocn => I_ALL
    implicit none

    integer :: n
    !---------------------------------------------------------------------------

    do n = 1, HIST_req_nmax

       if(      trim(item_save(n)) == 'sl_cld_frac'        ) out_cld_frac  = .true.

       if(      trim(item_save(n)) == 'sl_cldw'            &
           .OR. trim(item_save(n)) == 'sl_lwpqc'           &
           .OR. trim(item_save(n)) == 'sl_lwpqr'           &
           .OR. trim(item_save(n)) == 'sl_refmw'           ) out_cldw      = .true. 

       if(      trim(item_save(n)) == 'sl_qr'              &
           .OR. trim(item_save(n)) == 'sl_lwpqc'           &
           .OR. trim(item_save(n)) == 'sl_lwpqr'           &
           .OR. trim(item_save(n)) == 'sl_refmw'           ) out_qr        = .true.

       if(      trim(item_save(n)) == 'sl_qs'              ) out_qs        = .true.
       if(      trim(item_save(n)) == 'sl_qi'              ) out_qi        = .true.
       if(      trim(item_save(n)) == 'sl_qg'              ) out_qg        = .true.
       if(      trim(item_save(n)) == 'sl_cldi'            ) out_cldi      = .true.
       if(      trim(item_save(n)) == 'sl_slh'             ) out_slh       = .true.
       if(      trim(item_save(n)) == 'sl_tauucos'         ) out_tauucos   = .true.
       if(      trim(item_save(n)) == 'sl_tauvcos'         ) out_tauvcos   = .true.
       if(      trim(item_save(n)) == 'sl_ucos10m'         ) out_ucos10m   = .true.
       if(      trim(item_save(n)) == 'sl_vcos10m'         ) out_vcos10m   = .true.

       if(      trim(item_save(n)) == 'sl_tauu'            ) out_tauu      = .true.
       if(      trim(item_save(n)) == 'sl_tauv'            ) out_tauv      = .true.
       if(      trim(item_save(n)) == 'sl_u10m'            ) out_u10m      = .true.
       if(      trim(item_save(n)) == 'sl_v10m'            ) out_v10m      = .true.
       if(      trim(item_save(n)) == 'sl_vap_atm'         ) out_vap_atm   = .true.
       if(      trim(item_save(n)) == 'sl_tem_atm'         ) out_tem_atm   = .true.
       if(      trim(item_save(n)) == 'sl_slp'             ) out_slp       = .true.
       if(      trim(item_save(n)) == 'sl_cape'            &
           .OR. trim(item_save(n)) == 'sl_cin'             ) out_cape_cin  = .true.
       if(      trim(item_save(n)) == 'ml_ucos'            &
           .OR. trim(item_save(n)) == 'ml_vcos'            ) out_ucos_vcos = .true.
       if(      trim(item_save(n)) == 'ml_u'               &
           .OR. trim(item_save(n)) == 'ml_v'               ) out_u_v       = .true.
       if(      trim(item_save(n)) == 'vx'                 &
           .OR. trim(item_save(n)) == 'vy'                 &
           .OR. trim(item_save(n)) == 'vz'                 ) out_vx_vy_vz  = .true.
       if(      trim(item_save(n)) == 'ml_rh'              ) out_rh        = .true.
       if(      trim(item_save(n)) == 'ml_th'              ) out_th        = .true.
       if(      trim(item_save(n)) == 'ml_omg'             ) out_omg       = .true.
       if(      trim(item_save(n)) == 'ml_hgt'             ) out_hgt       = .true.
       if(      trim(item_save(n)) == 'sl_albedo'          ) out_albedo    = .true.

       if(      trim(item_save(n)) == 'dfq_isccp2_sw_toar' &
           .OR. trim(item_save(n)) == 'dfq_isccp2_lw_toa'  ) out_isccp_rad = .true.
       if(      trim(item_save(n)) == 're_cl'              &
           .OR. trim(item_save(n)) == 'nc_cl'              &
           .OR. trim(item_save(n)) == 'lwc_cl'             &
           .OR. trim(item_save(n)) == 'z_cl'               ) out_tauclk = .true.       

!       if ( trim(LAND_TYPE) == 'BUCKET' ) then
!          do idx=1, I_ALL_bc
!             if( trim(item_save(n)) == trim(GLVNAME_bc(idx)) ) out_landvars = .true.
!          end do
!       elseif(trim(LAND_TYPE) == 'MATSIRO') then
!          do idx=1, I_ALL_mat
!             if( trim(item_save(n)) == trim(GLVNAME_mat(idx)) ) out_landvars = .true.
!          enddo
!       endif
!
!       if ( trim(OCEAN_TYPE) == 'MIXEDLAYER' ) then
!          do idx=1, I_ALL_ocn
!             if( trim(item_save(n)) == trim(GOVNAME(idx)) ) out_oceanvars = .true.
!          enddo
!       endif

    enddo

    return
  end subroutine history_vars_setup

  !----------------------------------------------------------------------
  subroutine history_vars
    use mod_adm, only :         &
         ADM_LOG_FID,           &
         ADM_GALL_PL,           &
         ADM_LALL_PL,           &
         ADM_prc_me,            &
         ADM_prc_pl,            &
         ADM_gall,              &
         ADM_kall,              &
         ADM_lall,              &
         ADM_kmin,ADM_kmax,     &
         ADM_KNONE,             &
         ADM_gall_in             ! Y.Niwa add 070627
    use mod_gmtr, only :        &
         GMTR_area,             &
         GMTR_area_pl,          &
         GMTR_P_var,            &
         GMTR_P_var_pl,         &
         GMTR_P_IX,             &
         GMTR_P_IY,             &
         GMTR_P_IZ,             &
         GMTR_P_JX,             &
         GMTR_P_JY,             &
         GMTR_P_JZ,             &
         GMTR_P_LAT
    use mod_gtl, only :         &
         GTL_output_var2,       &
         GTL_output_var2_da,    &
         GTL_generate_uv,       &
         GTL_clip_region_1layer   ! [add] 2010.08.20 C.Kodama
    use mod_vmtr, only :        &
         VMTR_GSGAM2,           &
         VMTR_GSGAM2_pl,        &
         VMTR_GSGAM2H,          &
         VMTR_GSGAM2H_pl,       &
         VMTR_VOLUME,           &
         VMTR_VOLUME_pl
    use mod_prgvar, only :     &
         prgvar_get
    use mod_sfcvar, only :    &
         sfcvar_get,          &
         sfcvar_get1,         &
         sfcvar_get2,         &
         I_ALBEDO_SFC,        &
         I_PRECIP_TOT,        &
         I_PRECIP_CP,         &
         I_EVAP_SFC,          &
         I_PRE_SFC,           &
         I_TEM_SFC,           &
         I_SH_FLUX_SFC,       &
         I_LH_FLUX_SFC,       &
         I_SFCRAD_ENERGY,     &
         I_TOARAD_ENERGY,     &
         I_EVAP_ENERGY,       &
         I_PRECIP_ENERGY,     &
         I_TAUX_SFC,          &
         I_TAUY_SFC,          &
         I_TAUZ_SFC,          &
         I_RFLUXS_SD,        &
         I_RFLUXS_SU,        &
         I_RFLUXS_LD,        &
         I_RFLUXS_LU,        &
         I_RFLUX_TOA_SD,     &
         I_RFLUX_TOA_SU,     &
         I_RFLUX_TOA_LD,     &
         I_RFLUX_TOA_LU,     &
         I_RFLUX_TOA_SD_C,   &
         I_RFLUX_TOA_SU_C,   &
         I_RFLUX_TOA_LD_C,   &
         I_RFLUX_TOA_LU_C,   &
         I_CUMFRC,           &
         I_VX10,             &
         I_VY10,             &
         I_VZ10,             &
         I_T2,               &
         I_Q2,               &
         INDEX_SEA
    use mod_runconf, only : &
         LHV,LHS,           &
         RAIN_TYPE,         &
         CP_TYPE,           & ! 10/06/10 A.T.Noda
         NQW_MAX,           &
         TRC_VMAX,          &
         TRC_name,          &
         MP_TYPE,           & ! 07/05/08 H.tomita
         AE_TYPE,           & ! 07/07/30 K.Suzuki
         opt_2moment_water, & ! 08/05/24 T.Mitsui
         opt_incloud_aerosol,& ! 12/02/01 T.Seiki
         I_NC, I_NR, I_NI,  & ! 08/05/24 T.Mitsui
         I_NS, I_NG,        & ! 08/05/24 T.Mitsui
         NNW_STR, NNW_END,  & ! 08/05/24 T.Mitsui
         I_QV,              &
         I_QC,              &
         I_QR,              &
         I_QI,              &
         I_QS,              &
         I_QG,              &
         I_QH,              &
         NQW_STR,NQW_END,   &
         opt_carb_on, opt_dust_on, &! 12/02/01 T.Seiki
         opt_salt_on, opt_sulf_on, &! 12/02/01 T.Seiki
         NQDU_STR, NQDU_END,&       ! 12/02/01 T.Seiki
         NQDUIN_STR, NQDUIN_END, &  ! 12/02/01 T.Seiki
         NQCB_STR, NQCB_END, &      ! 12/02/01 T.Seiki
         NQCBIN_STR, NQCBIN_END, &  ! 12/02/01 T.Seiki
         NQSU_SO4,          &       ! 12/02/01 T.Seiki
         NQSU_SO2,          &       ! 12/02/01 T.Seiki
         NQSU_DMS,          &       ! 12/02/01 T.Seiki
         NQSO4IN,           &       ! 12/02/01 T.Seiki
         NQSA_STR, NQSA_END, &      ! 12/02/01 T.Seiki
         NQSAIN_STR, NQSAIN_END, &  ! 12/02/01 T.Seiki
         NRDIR,             &
         NRBND,             &
         NRBND_VIS,         &
         NRBND_NIR,         &
         NRBND_IR,          &
         NRDIR_DIRECT,      &
         NRDIR_DIFFUSE,     &
         NCRF,              &
         RAD_TYPE,          &  ! 07/12/05 T.Mitsui
         NTAU_ISCCP,        &  ! 07/12/05 T.Mitsui
         NPRES_ISCCP,       &  ! 07/12/05 T.Mitsui
         LAND_TYPE,         &  ! Y.Niwa add 070627
         OCEAN_TYPE            ! Y.Niwa add 070627
    use mod_bsstate, only : &
         phi, phi_pl
    use mod_cnvvar, only :  &
         cnvvar_p2d
    use mod_thrmdyn, only : &
         thrmdyn_th
    use mod_cnst, only : &
         CNST_EGRAV,     &
         CNST_RAIR,      &
         CNST_RVAP,      &
         CNST_TMELT,     &
         CNST_EPS_ZERO,  &
         CNST_UNDEF,     &
         CNST_DWATR
    use mod_bndcnd, only : &
         bndcnd_thermo
    use mod_diagvar, only : &
         diagvar_get,       &
         I_GDCFRC,          &
         I_CUMCLW,          &
         I_UNCCN                 ! 08/02/19 T.Mitsui
    use mod_grd, only:          & 
         grd_zsfc,              &
         grd_zs,                &
         grd_zs_pl,             &
         GRD_vz,                &
         GRD_Z,                 &
         GRD_VEGINDX              ! [add] 2010.08.20 C.Kodama
    use mod_history, only : &
       history_in


!    use mod_landvar_bc, only : &
!         landvar_getr1_in_bc => landvar_getr1_in, &
!         KMAX_bc => KMAX, &
!         I_ALL_bc => I_ALL, &
!         GLVNAME_bc => GLVNAME
!    use mod_landvar_matsiro, only : &
!         landvar_getr1_in_mat => landvar_getr1_in, &
!         KMAX_mat => KMAX,  &
!         I_ALL_mat => I_ALL, &
!         GLVNAME_mat => GLVNAME
!    use mod_oceanvar_mixedlayer, only : &
!         oceanvar_getr1_in_mixedlayer => oceanvar_getr1_in, &
!         KMAX_ocn => KMAX, &
!         I_ALL_ocn => I_ALL, &
!         GOVNAME
!    use mod_radvar, only :  &
!         radvar_get1,       &
!         radvar_get2,       &  ! 07/12/05 [Add] T.Mitsui
!         radvar_getr1,      &  ! 09/04/14 [Add] T.Mitsui
!         I_RCEFF,           &
!         I_RCEFF_CLD,       & ! 08/02/19 [Add] T.Mitsui
!         I_RWTOP,           &
!         I_RCTOP,           &
!         I_TAUCL,           &
!         I_TAUCI,           &
!         I_TAUCLK,          & ! 08/02/19 [Add] T.Mitsui
!         I_TAUCIK,          & ! 08/02/19 [Add] T.Mitsui
!         I_TAUC_ALL,        & ! 09/04/14 [Add] T.Mitsui
!         I_TCTOP,           &
!         I_DFQ_ISCCP           ! 07/12/05 [Add] T.Mitsui
    
    implicit none

    integer, parameter :: num = 0 

    real(8) :: rhog(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhoge(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogq(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    !
    real(8) :: pre(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rho(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: tem(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: th(ADM_gall,ADM_kall,ADM_lall)
    !
    real(8) :: q(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: vx(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vy(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vz(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: w(ADM_gall,ADM_kall,ADM_lall)
    !
    real(8) :: ucos(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vcos(ADM_gall,ADM_kall,ADM_lall)
    !
    real(8) :: rfluxs_sd(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rfluxs_su(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rfluxs_ld(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rfluxs_lu(ADM_gall,ADM_KNONE,ADM_lall)
    !
    real(8) :: rflux_toa_sd(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_su(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_ld(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_lu(ADM_gall,ADM_KNONE,ADM_lall)

    real(8) :: rflux_toa_sd_c(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_su_c(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_ld_c(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_lu_c(ADM_gall,ADM_KNONE,ADM_lall)
    !
    !
    real(8) :: precip(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: precip_cp(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: ptmp(ADM_gall,ADM_KNONE,ADM_lall,2)
    !
    real(8) :: sh_flux_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: lh_flux_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: sfcrad(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: toarad(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: precip_energy(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: evap_energy(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: taux(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: tauy(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: tauz(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: evap(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: pre_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: tem_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: vx10m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: vy10m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: vz10m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: t2m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: q2m(ADM_gall,ADM_KNONE,ADM_lall)
    !
    real(8) :: rh(ADM_gall,ADM_kall,ADM_lall)
    !
    real(8) :: cloud_frac(ADM_gall,ADM_kall,ADM_lall)
    !
    real(8) :: q_cumclw(ADM_gall,ADM_kall,ADM_lall)
    !
    real(8) :: CUMFRC(ADM_gall,ADM_KNONE,ADM_lall)
    !
    real(8) :: q_clw(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: q_cli(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: qtot(ADM_gall,ADM_kall,ADM_lall) ! qc+qr+qi+qs+qg [add] H.Yashiro 20120901
    !
    ! 09/04/14 T.Mitsui
    !--- for rd_mstrnx_ar5
    real(8) :: sl_tauc_all(ADM_gall,ADM_KNONE,ADM_lall,I_QC:I_QH)
    !
    real(8) :: rw(ADM_gall,ADM_lall)
    !
    real(8) :: albedo_sfc(ADM_gall,ADM_KNONE,ADM_lall,1:NRDIR,1:NRBND)
    !
    real(8) :: slp(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: zs(ADM_gall,ADM_lall)
    real(8) :: cape(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: cin(ADM_gall,ADM_KNONE,ADM_lall)
    ! 08/06/13 [Add] T.Mitsui, cape(positive only) and cin(negative only from kp to LNB)
    real(8) :: cape_posi(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: cin_all(ADM_gall,ADM_KNONE,ADM_lall)
    !
    ! 07/12/05 [Add] T.Mitsui
    ! notice, radvar_set2_in, get2_in and radvar_set2, get2 deal different dimension
    real(8) :: isccp1(ADM_gall,NTAU_ISCCP*NPRES_ISCCP,ADM_lall) ! isccp type data
    real(8) :: isccp2(ADM_gall,NTAU_ISCCP*NPRES_ISCCP,ADM_lall) ! isccp type data
    real(8) :: isccp3(ADM_gall,NTAU_ISCCP*NPRES_ISCCP,ADM_lall) ! isccp type data
    real(8) :: isccp1_pl(ADM_gall_pl,NTAU_ISCCP*NPRES_ISCCP,ADM_lall_pl) 
    !
    real(8) :: v2d(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: v3d(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: v2d_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL)
    real(8) :: v3d_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d2_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d3_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d4_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d5_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: v3d6_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
!
    real(8) :: ptmp_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL,2)
    real(8) :: albedo_sfc_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,1:NRDIR,1:NRBND) ! 05/11/13 S.Iga
    real(8) :: rhogq_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,TRC_VMAX)
    integer :: l,k,nq,n,ij
    integer :: ntau, npres ! 07/12/05 T.Mitsui
    !Y.Niwa add 070627 =>
    real(8) :: u(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: v(ADM_gall,ADM_kall,ADM_lall)
    !H.Tomita add 081114 =>
    real(8) :: omg(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: hgt(ADM_gall,ADM_kall,ADM_lall)
    !
    !--- for land variables
    integer, parameter :: LLKMAX = 6
    integer :: idx
    real(8) :: vll(ADM_gall_in,LLKMAX)
    !
    !--- for ocean variables
    integer, parameter :: OLKMAX = 2
    real(8) :: vol(ADM_gall_in,OLKMAX)
    !<= Y.Niwa add 070627
    ! 08/02/19, T.Mitsui
    real(8) :: rceff_cld(ADM_gall,ADM_kall,ADM_lall)
    ! for cloud and aerosol: 07/07/30 K.Suzuki add
    real(8) :: rceff(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: taucl(ADM_gall,ADM_KNONE,ADM_lall), tauci(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rwtop(ADM_gall,ADM_KNONE,ADM_lall), rctop(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: tctop(ADM_gall,ADM_KNONE,ADM_lall), lwpqc(ADM_gall), lwpqr(ADM_gall)
    ! 08/02/19 T.Mitsui, application for comparison among 
    ! Re derived by satellites using two-channel method.
    ! See Nakajima and King(1990), JAS, vol.47, no.15, p.1878-1893
    integer, parameter :: cot_binmax=11
    real(8) :: cot_bin(cot_binmax)=(/0.05d0,0.10d0,0.20d0,0.30d0,0.40d0, 0.50d0, 0.60d0, 0.70d0, 0.80d0, 0.90d0, 0.95d0/)
    real(8) :: tauclk(ADM_gall, ADM_kall, ADM_lall)     !  cloud optical thickness 
    real(8) :: nc(ADM_gall, ADM_kall, ADM_lall)         !  Number of Cloud Condensation Nuclei
    real(8) :: z_cl(ADM_gall, cot_binmax, ADM_lall)     !  height for each cloud layer
    real(8) :: re_cl(ADM_gall, cot_binmax, ADM_lall)    !  Effective Radius(qc+qr)
    real(8) :: rec_cl(ADM_gall, cot_binmax, ADM_lall)   !  Effective Radius(qc)
    real(8) :: nc_cl(ADM_gall, cot_binmax, ADM_lall)    !  Nc for each cloud layer
    real(8) :: lwc_cl(ADM_gall, cot_binmax, ADM_lall)   !  LWC for each cloud layer
    real(8) :: cwc_cl(ADM_gall, cot_binmax, ADM_lall)   !  CWC for each cloud layer
    !
    real(8) :: z_cb(ADM_gall, ADM_knone, ADM_lall)      !  height @ cloud base
    real(8) :: re_cb(ADM_gall, ADM_knone,  ADM_lall)    !  re     @ cloud base
    real(8) :: rec_cb(ADM_gall, ADM_knone,  ADM_lall)   !  rec    @ cloud base
    real(8) :: nc_cb(ADM_gall, ADM_knone,  ADM_lall)    !  nc     @ cloud base
    real(8) :: lwc_cb(ADM_gall, ADM_knone,  ADM_lall)   !  lwc    @ cloud base
    real(8) :: cwc_cb(ADM_gall, ADM_knone,  ADM_lall)   !  cwc    @ cloud base
    real(8) :: tau_work(ADM_gall)
    integer :: kctop(ADM_gall, ADM_KNONE)               ! cloud top layer
    integer :: kcbase(ADM_gall, ADM_KNONE)              ! cloud base layer
    logical :: flag_thick(ADM_gall, cot_binmax)
    integer :: ibin
    !
    real(8) :: gindex_sfc(ADM_gall_in,ADM_KNONE,ADM_lall) ! [add] 2010.08.20 C.Kodama
    integer :: index_sfc(ADM_gall_in,ADM_KNONE,ADM_lall)  !
    character(len=2)  :: wc ![Add] 2012/02/01 T.Seiki

    !--- H.Tomita 090414
    lwpqc(:)=0.0D0
    lwpqr(:)=0.0D0
    !<--- H.Tomita 090414
    !
    ! 2D variables
    call sfcvar_get(          &
         precip, v2d_pl,   &  !--- out
         vid = I_PRECIP_TOT   &  !--- in
         )
    call sfcvar_get1(         &
         ptmp, ptmp_pl,       &  !--- out
         vid = I_PRECIP_CP,   &  !--- in
         mdim = 2             &  !--- in
         )
    precip_cp(:,:,:) = ptmp(:,:,:,1) + ptmp(:,:,:,2)


    call sfcvar_get(         &
         sfcrad, v2d_pl,      & !--- IN 
         vid = I_SFCRAD_ENERGY   & !--- IN
         )
    call sfcvar_get(         &
         toarad, v2d_pl,      & !--- IN 
         vid = I_TOARAD_ENERGY   & !--- IN
         )
    call sfcvar_get(         &
         evap_energy, v2d_pl,      & !--- IN 
         vid = I_EVAP_ENERGY  & !--- IN
         )
    call sfcvar_get(         &
         precip_energy, v2d_pl,      & !--- IN 
         vid = I_PRECIP_ENERGY  & !--- IN
         )
    call sfcvar_get(         &
         sh_flux_sfc, v2d_pl,      & !--- IN 
         vid = I_SH_FLUX_SFC  & !--- IN
         )
    call sfcvar_get(         &
         lh_flux_sfc, v2d_pl,      & !--- IN 
         vid = I_LH_FLUX_SFC  & !--- IN
         )
    call sfcvar_get(          &
         taux,v2d_pl,        &  !--- out
         vid = I_TAUX_SFC     &  !--- in
         )
    call sfcvar_get(           &
         tauy,v2d_pl,         & !--- out
         vid = I_TAUY_SFC      & !--- in
         )
    call sfcvar_get(           &
         tauz,v2d_pl,         & !--- out
         vid = I_TAUZ_SFC      & !--- in
         )
    call sfcvar_get(           &
         evap, v2d_pl,        &  !--- out
         vid = I_EVAP_SFC      &  !--- in
         )
    call sfcvar_get(           &
         pre_sfc, v2d_pl,  &  !--- out
         vid = I_PRE_SFC       &  !--- in
         )
    call sfcvar_get(           &
         tem_sfc, v2d_pl,  &  !--- out
         vid = I_TEM_SFC       &  !--- in
         )
    call sfcvar_get(           &
         vx10m, v2d_pl,      &  !--- out
         vid = I_VX10          &  !--- in
         )
    call sfcvar_get(           &
         vy10m, v2d_pl,      &  !--- out
         vid = I_VY10          &  !--- in
         )
    call sfcvar_get(           &
         vz10m, v2d_pl,      &  !--- out
         vid = I_VZ10          &  !--- in
         )
    call sfcvar_get(           &
         q2m, v2d_pl,          &  !--- out
         vid = I_Q2            &  !--- in
         )
    call sfcvar_get(           &
         t2m, v2d_pl,          &  !--- out
         vid = I_T2            &  !--- in
         )
    call sfcvar_get(               &
         rfluxs_sd, v2d_pl,  &  !--- out
         vid = I_RFLUXS_SD         &  !--- in
         )
    call sfcvar_get(               &
         rfluxs_su, v2d_pl,  &  !--- out
         vid = I_RFLUXS_SU         &  !--- in
         )
    call sfcvar_get(               &
         rfluxs_ld, v2d_pl,  &  !--- out
         vid = I_RFLUXS_LD         &  !--- in
         )
    call sfcvar_get(               &
         rfluxs_lu, v2d_pl,  &  !--- out
         vid = I_RFLUXS_LU         &  !--- in
         )
    call sfcvar_get(                     &
         rflux_toa_sd, v2d_pl,  &  !--- out
         vid = I_RFLUX_TOA_SD            &  !--- in
         )
    call sfcvar_get(                     &
         rflux_toa_su, v2d_pl,  &  !--- out
         vid = I_RFLUX_TOA_SU            &  !--- in
         )
    call sfcvar_get(                     &
         rflux_toa_ld, v2d_pl,  &  !--- out
         vid = I_RFLUX_TOA_LD            &  !--- in
         )
    call sfcvar_get(                     &
         rflux_toa_lu, v2d_pl,  &  !--- out
         vid = I_RFLUX_TOA_LU            &  !--- in
         )
    call sfcvar_get(                     &
         rflux_toa_sd_c, v2d_pl,  &  !--- out
         vid = I_RFLUX_TOA_SD_C            &  !--- in
         )
    call sfcvar_get(                     &
         rflux_toa_su_c, v2d_pl,  &  !--- out
         vid = I_RFLUX_TOA_SU_C            &  !--- in
         )
    call sfcvar_get(                     &
         rflux_toa_ld_c, v2d_pl,  &  !--- out
         vid = I_RFLUX_TOA_LD_C            &  !--- in
         )
    call sfcvar_get(                     &
         rflux_toa_lu_c, v2d_pl,  &  !--- out
         vid = I_RFLUX_TOA_LU_C            &  !--- in
         )
    call sfcvar_get2   (  &
       albedo_sfc, albedo_sfc_pl, &  !--- OUT : surface variable
       I_ALBEDO_SFC,      &  !--- IN : variable ID
       NRDIR, NRBND       &  !--- IN
       )

    !--- get variables
    call prgvar_get(         &
         rhog, v3d_pl,      &  !--- out
         rhogvx, v3d2_pl,  &  !--- out
         rhogvy, v3d3_pl,  &  !--- out
         rhogvz, v3d4_pl,  &  !--- out
         rhogw,  v3d5_pl,   &  !--- out
         rhoge,  v3d6_pl,   &  !--- out
         rhogq,  rhogq_pl,   &  !--- out
         num )                  !--- in
!
!   if(trim(RAIN_TYPE)=='CLOUD_PARAM') then
    if(trim(CP_TYPE)/='NONE') then ! 10/06/10 A.T.Noda
       call diagvar_get(              &
            cloud_frac, v3d_pl,&
            vid = I_GDCFRC            &
            )
       call diagvar_get(              &
            q_cumclw, v3d_pl,    &
            vid = I_CUMCLW            &
            )
       call sfcvar_get(       &
            CUMFRC, v2d_pl,&
            vid = I_CUMFRC )
    else ! 05/11/08 M.Satoh
       q_cumclw(:,:,:) = 0.0d0
    end if
    !
    !--- convert variables
    do l=1,ADM_lall
       call cnvvar_p2d(         &
            ADM_gall,           &  !--- in
            rho(:,:,l),         &  !--- out
            pre(:,:,l),         &  !--- out
            tem(:,:,l),         &  !--- out
            vx(:,:,l),          &  !--- out
            vy(:,:,l),          &  !--- out
            vz(:,:,l),          &  !--- out
            w(:,:,l),           &  !--- out
            q(:,:,l,:),         &  !--- out
            rhog(:,:,l),        &  !--- in
            rhogvx(:,:,l),      &  !--- in
            rhogvy(:,:,l),      &  !--- in
            rhogvz(:,:,l),      &  !--- in
            rhogw(:,:,l),       &  !--- in
            rhoge(:,:,l),       &  !--- in
            rhogq(:,:,l,:),     &  !--- in
            VMTR_GSGAM2(:,:,l), &  !--- in
            VMTR_GSGAM2H(:,:,l) &  !--- in
            )
       call bndcnd_thermo( &
            ADM_gall,      &  !--- IN
            tem(:,:,l),    &  !--- INOUT
            rho(:,:,l),    &  !--- INOUT
            pre(:,:,l),    &  !--- INOUT
            phi(:,:,l)     &  !--- IN
            )     
       ! 09/04/14 T.Mitsui
       if(trim(RAIN_TYPE)=='CLOUD_PARAM') then
          q_clw(:,:,l) = q(:,:,l,I_QC)
          ! 09/04/14 T.Mitsui
       else !--- CRM
          ![Mod] T.Mitsui 08.02.19
          if( trim(MP_TYPE) == 'G98') then
          else
             ! 05/12/12 M.Satoh (suggested by T.Mitsui) 
             q_clw(:,:,l) = 0.0D0
             q_cli(:,:,l) = 0.0D0
             do nq = NQW_STR, NQW_END
                if ( nq == I_QC ) then
                   q_clw(:,:,l) = q_clw(:,:,l) + q(:,:,l,nq)
                else if ( nq == I_QR ) then
                   q_clw(:,:,l) = q_clw(:,:,l) + q(:,:,l,nq)
                else if ( nq == I_QI ) then
                   q_cli(:,:,l) = q_cli(:,:,l) + q(:,:,l,nq)
                   ! 09/04/14 T.Mitsui
                else if ( nq == I_QS ) then
                   q_cli(:,:,l) = q_cli(:,:,l) + q(:,:,l,nq)
                   ! 09/04/14 T.Mitsui
                else if ( nq == I_QG ) then
                   q_cli(:,:,l) = q_cli(:,:,l) + q(:,:,l,nq)
                   ! 09/04/14 T.Mitsui
                else if ( nq == I_QH ) then
                   q_cli(:,:,l) = q_cli(:,:,l) + q(:,:,l,nq)
                end if
             end do
          end if

          qtot(:,:,l) = q_clw(:,:,l) + q_cli(:,:,l)

          cloud_frac(:,:,l) = 0.0D0
          where((q_clw(:,:,l)+q_cli(:,:,l))>=0.005D-3) !--- tompkins & craig
             cloud_frac(:,:,l)=1.0D0
          end where
       end if
       !
       if(out_th) then  ! Y.Niwa add 070627
          call thrmdyn_th( &
               ADM_gall,   &  !--- in
               th(:,:,l),  &  !--- out
               tem(:,:,l), &  !--- in
               pre(:,:,l)  &  !--- in
               )
       end if ! Y.Niwa add 070627
       !
    end do
    !
    if(out_ucos_vcos) then ! Y.Niwa add 070627
       !--- zonal and meridonal wind generation with cos(phi)
       call GTL_generate_uv(&
            ucos,v3d_pl,   &
            vcos,v3d2_pl,   &
            vx, v3d3_pl,      &
            vy, v3d4_pl,      &
            vz, v3d5_pl,      &
            icos = 1)
    end if ! Y.Niwa add 070627
    !
    ! Y.Niwa add 070627 =>
    if(out_u_v) then 
       !--- zonal and meridonal wind generation
       call GTL_generate_uv( &
            u,  v3d_pl,   &
            v,  v3d2_pl,  &
            vx, v3d3_pl,  &
            vy, v3d4_pl,  &
            vz, v3d5_pl,  &
            icos = 0)
    end if
    !
    !--- 2D-1column ( Single Level )
    !
    if(out_slp) then ! Y.Niwa add 070627
       ! slp
       zs(:,:)= GRD_zs(:,ADM_KNONE,:,GRD_ZSFC)
       slp(:,ADM_KNONE,:) = pre_sfc(:,ADM_KNONE,:)  &
            *exp(cnst_egrav*zs(:,:)/                 &
            (cnst_rair*(t2m(:,ADM_KNONE,:)+0.5*0.005*(zs(:,:)+2))))
    end if ! Y.Niwa add 070627
    !
    do l=1,ADM_lall
!       call history_in ( 'sl_sw_toai', rflux_toa_sd(:,:,l) )  ! kmax = 1
!       call history_in ( 'sl_sw_toar', rflux_toa_su(:,:,l) )  ! kmax = 1
!       call history_in ( 'sl_lw_toa',  rflux_toa_lu(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_swd_toa', rflux_toa_sd(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_swu_toa', rflux_toa_su(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_lwu_toa',  rflux_toa_lu(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_lwd_toa',  rflux_toa_ld(:,:,l) )  ! kmax = 1
       
       !
       if(out_cldw) then ! Y.Niwa add 070627
          v2d(:,ADM_KNONE,l) = 0.0D0
          do k= ADM_kmin,ADM_kmax
             v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l) &
                  + rho(:,k,l)*VMTR_VOLUME(:,k,l)*q_clw(:,k,l)
          end do
          v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l)/GMTR_area(:,l)
          call history_in ( 'sl_cldw',  v2d(:,:,l) )  ! kmax = 1
          !
          ! [Del] 07/11/22 T.Mitsui change
!!$          if ( trim(AE_TYPE) == 'SPRINTARS' ) then ! K.Suzuki add 070730
          lwpqc( : ) = v2d( :,ADM_KNONE,l )
          do ij = 1, ADM_gall
             if ( tctop( ij,ADM_KNONE,l ) > CNST_TMELT ) then
                v2d( ij,ADM_KNONE,l ) = lwpqc( ij )
             else
                v2d( ij,ADM_KNONE,l ) = CNST_UNDEF
             end if
          end do
          call history_in( 'sl_lwpqc', v2d(:,:,l) )
          ! [Del] 07/11/22 T.Mitsui change
!!$       end if
          !
       end if ! Y.Niwa add 070627
       !
       if(out_qr) then ! Y.Niwa add 070627
          v2d(:,ADM_KNONE,l) = 0.0D0
          do k= ADM_kmin,ADM_kmax
             v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l) &
                  + rho(:,k,l)*VMTR_VOLUME(:,k,l)*q(:,k,l,I_QR)
          end do
          v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l)/GMTR_area(:,l)
          call history_in ( 'sl_qr',  v2d(:,:,l) )  ! kmax = 1
          !
          ! [Del] 07/11/22 T.Mitsui change
!!$          if ( trim(AE_TYPE) == 'SPRINTARS' ) then ! K.Suzuki add 070730
          !
          lwpqr( : ) = v2d( :,ADM_KNONE,l )
          do ij = 1, ADM_gall
             if ( tctop( ij,ADM_KNONE,l ) > CNST_TMELT ) then
                v2d( ij,ADM_KNONE,l ) = lwpqr( ij )
             else
                v2d( ij,ADM_KNONE,l ) = CNST_UNDEF
             end if
          end do
          call history_in ( 'sl_lwpqr', v2d(:,:,l) )
          !
          do ij = 1, ADM_gall
             if ( tctop( ij,ADM_KNONE,l ) > CNST_TMELT ) then
                v2d( ij,ADM_KNONE,l ) = lwpqc( ij ) + lwpqr( ij )
             else
                v2d( ij,ADM_KNONE,l ) = CNST_UNDEF
             end if
          end do
          call history_in ( 'sl_lwpmw', v2d(:,:,l) )
          !
          ! [Del] 07/11/22 T.Mitsui change
!!$          end if
          !
       end if ! Y.Niwa add 070627
       !
       !
       if(out_cldi) then ! Y.Niwa add 070627
          v2d(:,ADM_KNONE,l) = 0.0D0
          do k= ADM_kmin,ADM_kmax
             v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l) &
                  + rho(:,k,l)*VMTR_VOLUME(:,k,l)*q_cli(:,k,l)
          end do
          v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l)/GMTR_area(:,l)
          call history_in ( 'sl_cldi',  v2d(:,:,l) )  ! kmax = 1
       end if ! Y.Niwa add 070627
       call history_in ( 'sl_tppn',  precip(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_evap',  evap(:,:,l) )  ! kmax = 1
!       call history_in ( 'sl_swi',  rfluxs_sd(:,:,l) )  ! kmax = 1
!       call history_in ( 'sl_sswr',  rfluxs_su(:,:,l) )  ! kmax = 1
!       call history_in ( 'sl_slwd',  rfluxs_ld(:,:,l) )  ! kmax = 1
!       call history_in ( 'sl_slwu',  rfluxs_lu(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_swd_sfc',  rfluxs_sd(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_swu_sfc',  rfluxs_su(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_lwd_sfc',  rfluxs_ld(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_lwu_sfc',  rfluxs_lu(:,:,l) )  ! kmax = 1
       !
!        if(out_slh) then ! Y.Niwa add 070627 
! !           where(t2m(:,ADM_KNONE,l)>CNST_TMELT)
! !              v2d(:,ADM_KNONE,l) = evap(:,ADM_KNONE,l) * LHV
! !           elsewhere
! !              v2d(:,ADM_KNONE,l) = evap(:,ADM_KNONE,l) * LHS
! !           end where
! !           !
!           v2d(:,ADM_KNONE,l) = lh_flux_sfc(:,ADM_KNONE,l)
! !          call history_in ( 'sl_slh',  v2d(:,:,l) )  ! kmax = 1
!           call history_in ( 'sl_lh_sfc',  v2d(:,:,l) )! kmax = 1
!        end if ! Y.Niwa add 070627
       !
       call history_in ( 'sl_lh_sfc',  lh_flux_sfc(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_sh_sfc',  sh_flux_sfc(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_evap_energy',  evap_energy(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_tppn_energy',  precip_energy(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_rad_toa',  toarad(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_rad_sfc',  sfcrad(:,:,l) )  ! kmax = 1

       call history_in ( 'sl_t2m', t2m(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_q2m', q2m(:,:,l) )  ! kmax = 1
       !
       if(out_tauucos) then ! Y.Niwa add 070627
          !---sl_tauucos       with cos(lat)
          v2d(:,ADM_KNONE,l) &
               = (taux(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IX) &
               + tauy(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IY) &
               + tauz(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IZ))&
               * cos(GMTR_P_var(:,ADM_KNONE,l,GMTR_P_LAT))
          call history_in ( 'sl_tauucos', v2d(:,:,l) )  ! kmax = 1
          !
       end if  ! Y.Niwa add 070627
       if(out_tauu) then ! W.Yanase  08/07/25
          v2d(:,ADM_KNONE,l) &
               = (taux(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IX) &
               + tauy(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IY) &
               + tauz(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IZ))
          call history_in ( 'sl_tauu', v2d(:,:,l) )  ! kmax = 1
       end if  
       if(out_tauvcos) then  ! Y.Niwa add 070627
          !---sl_tauvcos       with cos(lat)
          v2d(:,ADM_KNONE,l) &
               = (taux(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JX) &
               + tauy(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JY) &
               + tauz(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JZ))&
               * cos(GMTR_P_var(:,ADM_KNONE,l,GMTR_P_LAT))
          call history_in ( 'sl_tauvcos', v2d(:,:,l) )  ! kmax = 1
          !
       end if   ! Y.Niwa add 070627
       if(out_tauv) then  ! W.Yanase   
          v2d(:,ADM_KNONE,l) &
               = (taux(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JX) &
               + tauy(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JY) &
               + tauz(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JZ))
          call history_in ( 'sl_tauv', v2d(:,:,l) )  ! kmax = 1
       end if   
       if(out_ucos10m) then  ! Y.Niwa add 070627
          !---sl_ucos10m       with cos(lat)
             v2d(:,ADM_KNONE,l) &
                  = (vx10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IX) &
                  + vy10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IY) &
                  + vz10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IZ))&
                  * cos(GMTR_P_var(:,ADM_KNONE,l,GMTR_P_LAT))
             call history_in ( 'sl_ucos10m', v2d(:,:,l) )  ! kmax = 1
          !
       end if  ! Y.Niwa add 070627
       if(out_u10m) then  ! W. Yanase   08/07/25
             v2d(:,ADM_KNONE,l) &
                  = (vx10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IX) &
                  + vy10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IY) &
                  + vz10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IZ))
             call history_in ( 'sl_u10m', v2d(:,:,l) )  
       end if  
       if(out_vcos10m) then  ! Y.Niwa add 070627
          !---sl_vcos10m       with cos(lat)
          v2d(:,ADM_KNONE,l) &
               = (vx10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JX) &
               + vy10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JY) &
               + vz10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JZ))&
               * cos(GMTR_P_var(:,ADM_KNONE,l,GMTR_P_LAT))
          call history_in ( 'sl_vcos10m', v2d(:,:,l) )  ! kmax = 1
          !
       end if  ! Y.Niwa add 070627
       if(out_v10m) then  ! W.Yanase
          v2d(:,ADM_KNONE,l) &
               = (vx10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JX) &
               + vy10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JY) &
               + vz10m(:,ADM_KNONE,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JZ))
          call history_in ( 'sl_v10m', v2d(:,:,l) )  ! kmax = 1
          !
       end if 
       !
       !
       call history_in ( 'sl_ps',  pre_sfc(:,:,l) )  ! kmax = 1
       !
       if(out_vap_atm) then ! Y.Niwa add 070627
          !--- sl_vap_atm      
          v2d(:,ADM_KNONE,l) = 0.0D0
          do k= ADM_kmin,ADM_kmax
             v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l) &
                  + rho(:,k,l)*VMTR_VOLUME(:,k,l)*q(:,k,l,I_QV)
          end do
          v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l)/GMTR_area(:,l)
          call history_in ( 'sl_vap_atm', v2d(:,:,l) )  ! kmax = 1
          !
       end if  ! Y.Niwa add 070627
       !
       if(out_tem_atm) then ! Y.Niwa add 070627
          !--- sl_tem_atm       
          v2d(:,ADM_KNONE,l) = 0.0D0
          rw(:,:) = 0.0D0 
          do k= ADM_kmin,ADM_kmax
             v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l) &
                  + rho(:,k,l)*VMTR_VOLUME(:,k,l)*tem(:,k,l)
             rw(:,l) = rw(:,l) + rho(:,k,l)*VMTR_VOLUME(:,k,l)
          end do
          v2d(:,ADM_KNONE,l) = v2d(:,ADM_KNONE,l)/rw(:,l)
          call history_in ( 'sl_tem_atm', v2d(:,:,l) )  ! kmax = 1
          !
       end if  ! Y.Niwa add 070627
       !
       call history_in ( 'sl_tem_sfc', tem_sfc(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_cppn', precip_cp(:,:,l) )  ! kmax = 1
!       call history_in ( 'sl_lw_toa_c', rflux_toa_lu_c(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_lwu_toa_c', rflux_toa_lu_c(:,:,l) )  ! kmax = 1
!       call history_in ( 'sl_sw_toar_c', rflux_toa_su_c(:,:,l) )  ! kmax = 1
       call history_in ( 'sl_swu_toa_c', rflux_toa_su_c(:,:,l) )  ! kmax = 1
       
       ! Y.Niwa add 070627 suggested by M.Satoh
       if(out_albedo) then
          v2d(:,ADM_KNONE,l) &
               = min( rflux_toa_su(:,ADM_KNONE,l) &
               / max( rflux_toa_sd(:,ADM_KNONE,l), CNST_EPS_ZERO ), 1.0d0 )
          call history_in ( 'sl_albedo', v2d(:,:,l) ) ! kmax = 1
       end if
       !
       call history_in ( 'sl_albedo_sfc', albedo_sfc(:,:,l,NRDIR_DIRECT,NRBND_VIS) )  ! kmax = 1
       !
       if(out_slp) then ! Y.Niwa add 070627
          call history_in ( 'sl_slp', slp(:,:,l) )  ! kmax = 1
       end if ! Y.Niwa add 070627
       !
       if(out_cape_cin) then ! Y.Niwa add 070627
          call history_in ( 'sl_cape',cape(:,:,l) )  ! kmax = 1
          call history_in ( 'sl_cin', cin(:,:,l) )  ! kmax = 1
          !
          call history_in ( 'sl_cape_posi',cape_posi(:,:,l) )  ! kmax = 1
          call history_in ( 'sl_cin_all', cin_all(:,:,l) )  ! kmax = 1
       end if
       !
       ! aerosol parameters
       ! [Del] 07/11/22 T.Mitsui change
!!$       if ( trim(AE_TYPE) == 'SPRINTARS' ) then ! 07/07/30 K.Suzuki add
       ! rceff
       call history_in ( 'ml_rceff', rceff(:,:,l) )
       ! [Add] 09/04/14 T.Mitsui
       ! cloud optical thickness for each species
       call history_in ( 'sl_tauc', sl_tauc_all(:,:,l,I_QC) )
       call history_in ( 'sl_taur', sl_tauc_all(:,:,l,I_QR) )
       call history_in ( 'sl_taui', sl_tauc_all(:,:,l,I_QI) )
       call history_in ( 'sl_taus', sl_tauc_all(:,:,l,I_QS) )
       call history_in ( 'sl_taug', sl_tauc_all(:,:,l,I_QG) )
       call history_in ( 'sl_tauh', sl_tauc_all(:,:,l,I_QH) )       
       ! taucl: optical thickness of liquid cloud
       call history_in ( 'sl_taucl', taucl(:,:,l) )
       
       ! tauci: optical thickness of ice cloud
       call history_in ( 'sl_tauci', tauci(:,:,l) )
       
       ! taucw: COT only for Tc > 273.15K
       do ij = 1, ADM_gall
          if ( tctop( ij,ADM_KNONE,l ) > CNST_TMELT ) then
             v2d(ij,ADM_KNONE,l) = taucl( ij,ADM_KNONE,l )
          else
             v2d(ij,ADM_KNONE,l) = CNST_UNDEF
          end if
       end do
       call history_in ( 'sl_taucw', v2d(:,:,l) )
       
       ! rwtop: effective particle radius at cloud top > 273.15K
       call history_in ( 'sl_rwtop', rwtop(:,:,l) )
       
       ! rctop: effective particle radius at cloud top
       call history_in ( 'sl_rctop', rctop(:,:,l) )
       
       ! tctop: cloud top temperature
       call history_in ( 'sl_tctop', tctop(:,:,l) )
       
       ! column effective radius
       do ij = 1, ADM_gall
          ![Mod] 07/11/22 T.Mitsui threshold based on ISCCP_simulator
          ! tau < 0.3 : sub-visible cloud ?
!!$          if ( taucl( ij,ADM_KNONE,l ) > 1.D0 .and. &
          if ( taucl( ij,ADM_KNONE,l ) > 0.3D0 .and. &
               tctop( ij,ADM_KNONE,l ) > CNST_TMELT ) then
             v2d( ij,ADM_KNONE,l ) = 1.5D0*( lwpqc( ij )+lwpqr( ij ) ) &
                  /taucl( ij,ADM_KNONE,l )/CNST_DWATR
          else
             v2d( ij,ADM_KNONE,l ) = CNST_UNDEF
          end if
       end do
       call history_in ( 'sl_refmw', v2d(:,:,l) )
       ! [Del] 07/11/22 T.Mitsui change    
!!$      end if
    
    end do
    !
    ! 3D variables
    do l=1,ADM_lall
       if(out_ucos_vcos) then ! Y.Niwa add 070627
          call history_in ( 'ml_ucos', ucos(:,:,l) )  
          call history_in ( 'ml_vcos', vcos(:,:,l) ) 
       end if ! Y.Niwa add 070627
       ! Y.Niwa add 070627 =>
       if(out_u_v) then 
          call history_in ( 'ml_u', u(:,:,l) )
          call history_in ( 'ml_v', v(:,:,l) )
       end if
       if(out_vx_vy_vz) then
          call history_in ( 'vx', vx(:,:,l) )
          call history_in ( 'vy', vy(:,:,l) )
          call history_in ( 'vz', vz(:,:,l) )
       end if
       !<= Y.Niwa add 070627
       !
       do k=ADM_kmin,ADM_kmax
          v3d(:,k,l) = 0.5D0*(w(:,k,l)+w(:,k+1,l))
       end do
       !
       !
       !==== bug fix
       v3d(:,ADM_kmin-1,l) = 0.0D0
       v3d(:,ADM_kmax+1,l) = 0.0D0
       !
       if(out_omg) then
          omg(:,:,l) = -rho(:,:,l)*CNST_EGRAV*v3d(:,:,l)
       end if
       if(out_hgt) then
          hgt(:,:,l) = GRD_vz(:,:,l,GRD_Z)
       end if
       ! <= H.Tomita add 081114 : NOTE Hydrostatic assumption
       !
       call history_in ( 'ml_w', v3d(:,:,l) ) 
       call history_in ( 'ml_omg', omg(:,:,l) ) 
       call history_in ( 'ml_hgt', hgt(:,:,l) ) 
       call history_in ( 'ml_rho', rho(:,:,l) ) 
       call history_in ( 'ml_tem', tem(:,:,l) ) 
       call history_in ( 'ml_pres', pre(:,:,l) ) 
       !
       v3d(:,:,l)= 0.0D0
       if(NQW_MAX==6) then !--- cold rain
          call history_in ( 'ml_qv', q(:,:,l,I_QV) ) 
          call history_in ( 'ml_qc', q(:,:,l,I_QC) ) 
          call history_in ( 'ml_qr', q(:,:,l,I_QR) ) 
          call history_in ( 'ml_qi', q(:,:,l,I_QI) ) 
          call history_in ( 'ml_qs', q(:,:,l,I_QS) ) 
          call history_in ( 'ml_qg', q(:,:,l,I_QG) ) 
          ! 08/05/24 [Add] T.Mitsui
          if( opt_2moment_water ) then
             call history_in ( 'ml_nc', rho(:,:,l)*q(:,:,l,I_NC) ) 
             call history_in ( 'ml_nr', rho(:,:,l)*q(:,:,l,I_NR) ) 
             call history_in ( 'ml_ni', rho(:,:,l)*q(:,:,l,I_NI) ) 
             call history_in ( 'ml_ns', rho(:,:,l)*q(:,:,l,I_NS) ) 
             call history_in ( 'ml_ng', rho(:,:,l)*q(:,:,l,I_NG) )
          end if
          !
       else if(NQW_MAX==3) then
          call history_in ( 'ml_qv', q(:,:,l,I_QV) ) 
          call history_in ( 'ml_qc', q(:,:,l,I_QC) ) 
          call history_in ( 'ml_qr', q(:,:,l,I_QR) ) 
          !
          ![Mod] T.Mitsui 08.02.19
          if(trim(MP_TYPE) == 'G98') then
!             call history_in ( 'ml_qg', v3d(:,:,l) ) 
             !
          else
!             call history_in ( 'ml_qi', v3d(:,:,l) ) 
!             call history_in ( 'ml_qs', v3d(:,:,l) ) 
!             call history_in ( 'ml_qg', v3d(:,:,l) ) 
             !
          end if
       else if(NQW_MAX==2) then
          call history_in ( 'ml_qv', q(:,:,l,I_QV) ) 
          call history_in ( 'ml_qc', q(:,:,l,I_QC) ) 
!          call history_in ( 'ml_qr', v3d(:,:,l) ) 
!          call history_in ( 'ml_qi', v3d(:,:,l) ) 
!          call history_in ( 'ml_qs', v3d(:,:,l) ) 
!          call history_in ( 'ml_qg', v3d(:,:,l) ) 
          !
       endif
       call history_in ( 'ml_cumclw', q_cumclw(:,:,l) ) 
       call history_in ( 'ml_qtot', qtot(:,:,l) ) ! [add] H.Yashiro 20120901
       if(out_rh) then ! Y.Niwa add 070627
          call history_in ( 'ml_rh', rh(:,:,l) )
       end if ! Y.Niwa add 070627
       if(out_th) then ! Y.Niwa add 070627
          call history_in ( 'ml_th', th(:,:,l) ) 
       end if ! Y.Niwa add 070627
       !
       !
       ! 2012/02/01 T.Seiki
       !
       if( opt_salt_on )then
          n=0
          do nq=NQSA_STR, NQSA_END
             n=n+1
             write(wc,'(i2.2)') n
             call history_in ('ml_qa_sa'//trim(wc), q(:,:,l,nq))             
          end do
          if( opt_incloud_aerosol )then
             n=0
             do nq=NQSAIN_STR, NQSAIN_END
                n=n+1
                write(wc,'(i2.2)') n
                call history_in ('ml_qain_sa'//trim(wc), q(:,:,l,nq))             
             end do
          end if
       end if
       !
       if( opt_dust_on )then
          n=0
          do nq=NQDU_STR, NQDU_END
             n=n+1
             write(wc,'(i2.2)') n
             call history_in ('ml_qa_du'//trim(wc), q(:,:,l,nq))             
          end do
          if( opt_incloud_aerosol )then
             n=0
             do nq=NQDUIN_STR, NQDUIN_END
                n=n+1
                write(wc,'(i2.2)') n
                call history_in ('ml_qain_du'//trim(wc), q(:,:,l,nq))             
             end do
          end if
       end if
       !
       if( opt_carb_on )then
          n=0
          do nq=NQCB_STR, NQCB_END
             n=n+1
             write(wc,'(i2.2)') n
             call history_in ('ml_qa_cb'//trim(wc), q(:,:,l,nq))             
          end do
          if( opt_incloud_aerosol )then
             n=0
             do nq=NQCBIN_STR, NQCBIN_END
                n=n+1
                write(wc,'(i2.2)') n
                call history_in ('ml_qain_cb'//trim(wc), q(:,:,l,nq))             
             end do
          end if
       end if
       !
       if( opt_sulf_on )then
          call history_in ('ml_qa_so4', q(:,:,l,NQSU_SO4))
          call history_in ('ml_qa_so2', q(:,:,l,NQSU_SO2))
          call history_in ('ml_qa_dms', q(:,:,l,NQSU_DMS))
          if( opt_incloud_aerosol )then
             call history_in ('ml_qain_so4', q(:,:,l,NQSO4IN))          
          end if
       end if

       !=> [ADD] 20121107 H.Yashiro
       do nq = 1, TRC_vmax
          call history_in ( TRC_name(nq), q(:,:,l,nq) )
       enddo

    enddo

!    if(out_landvars) then
!       call GTL_clip_region_1layer( GRD_zs(:,ADM_KNONE,:,GRD_VEGINDX), gindex_sfc)
!       index_sfc(:,:,:) = nint ( gindex_sfc(:,:,:) )
!
!       if(trim(LAND_TYPE) == 'BUCKET') then
!          do l=1, ADM_lall
!             do idx=1, I_ALL_bc
!                call landvar_getr1_in_bc(    &
!                     vll(:,1:KMAX_bc(idx)),  &  !--- out
!                     idx,                    &  !--- in
!                     KMAX_bc(idx),           &  !--- in
!                     l                       &  !--- in
!                     )
!                do k=1, KMAX_bc(idx)
!                   where ( index_sfc(:,ADM_KNONE,l) <= INDEX_SEA )
!                      vll(:,k) = CNST_UNDEF
!                   end where
!                end do
!
!                call history_in( trim(GLVNAME_bc(idx)), vll(:,1:KMAX_bc(idx)) )
!                vll = 0.0
!
!             end do
!          end do
!
!       else if(trim(LAND_TYPE) == 'MATSIRO') then
!          do l=1, ADM_lall
!             do idx=1, I_ALL_mat
!                call landvar_getr1_in_mat(     &
!                     vll(:,1:KMAX_mat(idx)),  &  !--- out
!                     idx,                     &  !--- in
!                     KMAX_mat(idx),           &  !--- in
!                     l                        &  !--- in
!                  )
!                do k=1, KMAX_mat(idx)
!                   where ( index_sfc(:,ADM_KNONE,l) <= INDEX_SEA )
!                      vll(:,k) = CNST_UNDEF
!                   end where
!                end do
!                call history_in( trim(GLVNAME_mat(idx)), vll(:,1:KMAX_mat(idx)) )
!                vll = 0.0
!                !
!             end do
!          end do
!       endif
!    endif
!
!    if ( out_oceanvars ) then
!       if ( trim(OCEAN_TYPE) == 'MIXEDLAYER' ) then
!          do l = 1, ADM_lall
!             do idx = 1, I_ALL_ocn
!                vol(:,:) = 0.D0
!                call oceanvar_getr1_in_mixedlayer( vol(:,1:KMAX_ocn(idx)), idx, KMAX_ocn(idx), l )
!                call history_in( trim(GOVNAME(idx)), vol(:,1:KMAX_ocn(idx)), l )
!             enddo
!          enddo
!       endif
!    endif

    return
  end subroutine history_vars

end module mod_history_vars
!-------------------------------------------------------------------------------
