!-------------------------------------------------------------------------------
!
!+  Module for diagnosys
!
!-------------------------------------------------------------------------------
module mod_diagvar
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module provides diagnosys
  !       in non-hydrostatic model.
  !       
  ! 
  !++ Current Corresponding Author : M.Satoh
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-11-02   NICAM nicam_start041102.tgz
  !      0.01      04-12-07   M.Satoh: add diagvar_checkvalues
  !                05-11-15   M.Satoh: diagvar_comm
  !                07-07-17   A.T.Noda: 
  !                           1. Add I_QKEd, I_TSQd, I_QSQd, I_COVd, I_CFRACP 
  !                              I_QCLWB, I_QCLIWB for turbulent module.
  !                           2. Add diagvar_timeinfo
  !                07-07-23   K.Suzuki: Add I_DFE for use in SPRINTARS
  !                07-11-06   T.Mitsui: Save memory and add flag_diagset(debug for SPRINTARS)
  !                07-12-05   T.Mitsui: add UNCCN for restart and diag_cname(order cname)
  !                08-03-10   T.Mitsui: add output of intermediate restart file
  !                08-05-30   T.Mitsui: call diagvar_comm in diagvar_restart_output
  !                09-01-28   A.T.Noda: Implement mynn
  !                09-08-18   T.Mitsui: clear up intent(in) for -debug option @ SR11000
  !                10-04-26   M.Satoh: add diagvar1, 
  !                                    diagvar_get_in_1layer,
  !                                    diagvar_set_in_1layer
  !                                    diagvar_get_in_1layer_region,
  !                                    diagvar_set_in_1layer_region
  !                10-05-05   M.Satoh: add 
  !                                    diagvar_get_in_region,
  !                                    diagvar_set_in_region
  !                10-06-19   A.T.Noda:
  !                             1. Allow to use a convection parameterization
  !                                with advanced microphysics schemes (G98, NSW?, ...)
  !                             2. bug fix in call of diagvar_timeinfo
  !                10-08-03   T.Mitsui: fix for gfortran
  !                11-05-07   Y.Yamada: Implementation of ES tuning cord by NEC.
  !                             Modified line: (2011/03/02 NEC)
  !                11-08-16   M.Satoh: bug fix for TDK: conv => tendency
  !                11-08-16   M.Satoh: implementation for Chikira scheme: 
  !                                    introduce diagvarn for n-layers diagnostic variables
  !                11-08-16b  M.Satoh : update TIEDKE, cp_driver
  !                11-08-17   M.Satoh : bug fix
  !                11-08-30   A.Noda  : bug fix
  !                11-09-03   H.Yashiro : New I/O
  !                11-11-02   H.Yashiro : bug fix
  !                11-11-25   Y.Yamada : bug fix
  !                12-03-09   S.Iga: tuned (phase4-1)
  !                12-03-28   T.Seiki  : bug fix (doesn't affect at this moment)
  !                12-06-07   T.Seiki  : trivial modification for multi-job run
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only : &
       ADM_MAXFNAME,  &
       ADM_LOG_FID,       &
       ADM_CTL_FID,       &
       ADM_GALL_PL,       &
       ADM_LALL_PL,       &
       ADM_prc_me,        &
       ADM_prc_pl,        &
       ADM_kall,          &
       ADM_KNONE,         &
       ADM_gall,          &
       ADM_gmin,          &
       ADM_gmax,          &
       ADM_lall,          &
       ADM_gall_1d,       &
       ADM_IopJop_nmax,   &
       ADM_IopJop,        &
       ADM_GIoJo,         &
       ADM_rgnid_npl_mng, &
       ADM_rgnid_spl_mng, &
       ADM_prc_tab,       &
       ADM_prc_npl,       &
       ADM_prc_spl,       &
       ADM_GSLF_PL,       &
       ADM_NPL,           &
       ADM_SPL,           &
       ADM_COMM_RUN_WORLD,&
       ADM_rgn2prc,       &
       ADM_proc_stop,     &
       ADM_NSYS
  use mod_comm, only : &
       COMM_data_transfer

  implicit none
  
  logical, private,save ::  flag_diagset           ! [Add] 07.11.06 T.Mitsui
  !
  integer, public, save ::  DIAG_VMAX
  integer, private      ::  I_STA          
  integer, private      ::  I_END           

  ! [Add] 10.04.26 M.Satoh
  integer, public, save ::  DIAG_VMAX_1LAYER
  integer, private      ::  I_STA_1LAYER ! 2010.5.5 M.Satoh [add]
  integer, private      ::  I_END_1LAYER

  ! [Add] 11.07.25 M.Satoh
  integer, public, save ::  DIAG_VMAX_NLAYER
  integer, private      ::  I_STA_NLAYER
  integer, private      ::  I_END_NLAYER
  integer, public, save ::  DIAG_KTOT_NLAYER

  ! only for CP_TYPE == CLOUD_PARAM
  ! 07/12/05 [Mod] T.Mitsui
!!$  integer, private,save ::  I_CLOUD_PARAM_VMAX ! 11/08/16b M.Satoh [del]
  integer, public, save ::  I_CUMCLW       = -999
  integer, public, save ::  I_GDCLW        = -999
  integer, public, save ::  I_GDCFRC       = -999
  integer, public, save ::  I_CBMFX        = -999
  integer, public, save ::  I_RHOGQV_CONV  = -999
  integer, public, save ::  I_QV_DYN_TEND  = -999 ! 11/08/16 M.Satoh
  integer, public, save ::  I_QV_TB_TEND  = -999  ! 11/08/16 M.Satoh
!!$  integer, public, save ::  I_RHOGQV_CONV_TB = -999 ! 10/05/06 M.Satoh
  integer, public, save ::  I_EVAP_SFC     = -999 ! 10/05/22 M.Satoh
  integer, public, save ::  I_SH_FLUX_SFC  = -999 ! 2011/08/16b M.Satoh
  ! 11/08/16 M.Satoh
  !   I_CBMFX_CHIKIRA is specific to Chikira scheme, 
  !   while I_CBMFX is generic for any cumulus parameterization.
  !   The variable for I_CBMFX_CHIKIRA has NCTP layers, 
  !   while the variable for I_CBMFX has ADM_kall layers.
  integer, public, save ::  I_CBMFX_CHIKIRA = -999 ! 11/08/16 M.Satoh
  ! only for TB_TYPE == MY2MOIST 
  ! 07/12/05 [Mod] T.Mitsui
!!$  integer, private, parameter :: I_MY2MOIST_VMAX = 7
! integer, private,save ::  I_MY2MOIST_VMAX
  integer, private,save ::  I_TURB_VMAX
  integer, public, save ::  I_QKEd         = -999
  integer, public, save ::  I_TSQd         = -999
  integer, public, save ::  I_QSQd         = -999
  integer, public, save ::  I_COVd         = -999
  integer, public, save ::  I_CFRACP       = -999
  integer, public, save ::  I_QCLWB        = -999
  integer, public, save ::  I_QCLIB        = -999
  ! only for AE_TYPE == SPRINTARS
  ! 07/12/05 [Mod] T.Mitsui
!!$  integer, private, parameter :: I_SPRINTARS_VMAX = 1
  integer, private,save ::  I_SPRINTARS_VMAX
  integer, public, save ::  I_DFE          = -999
  ! 07/12/05 [Add] T.Mitsui
  ! for Cloudmicrophysics
  integer, private,save ::  I_MP_VMAX
  integer, public, save ::  I_UNCCN        = -999

  ! 10/04/26 [Add] M.Satoh
  ! for roughness_seq_yqw
!!$  integer, private,save ::  I_ROUGHNESS_SEA_MAX    ! 10/05/05 [add] M.Satoh; 11/08/16 M.Satoh [del]
  integer, public, save ::  I_ROUGHNESS_SEA = -999 ! 10/04/26 [add] M.Satoh
  !
  real(8), public, allocatable, save :: diagvar(:,:,:,:)
  real(8), public, allocatable, save :: diagvar_pl(:,:,:,:)

  real(8), public, allocatable, save :: diagvar1(:,:,:,:) ! 10/04/26 [add] M.Satoh
  real(8), public, allocatable, save :: diagvar1_pl(:,:,:,:) ! 10/04/26 [add] M.Satoh
  real(8), public, allocatable, save :: diagvarn(:,:,:,:) ! 11/08/16 [add] M.Satoh
  real(8), public, allocatable, save :: diagvarn_pl(:,:,:,:) ! 11/08/16 [add] M.Satoh
  !
  integer, private, parameter :: DIAG_VMAX_DEF = 256 ! 2010/05/05 M.Satoh [add]
  character(len=ADM_NSYS), public, save :: & 
       diag_cname(DIAG_VMAX_DEF)             ! 07/12/05 [Add] T.Mitsui
  character(len=ADM_NSYS), public, save :: &
       diag_cname_1layer(DIAG_VMAX_DEF)      ! 10/04/29 [Add] M.Satoh
  character(len=ADM_NSYS), public, save :: &
       diag_cname_nlayer(DIAG_VMAX_DEF)      ! 11/08/16 [Add] M.Satoh
  integer, public, save :: &
       diag_knum_nlayer(DIAG_VMAX_DEF)       ! 11/08/16 [Add] M.Satoh
  integer, public, save :: &
       diag_ksta_nlayer(DIAG_VMAX_DEF)       ! 11/08/16 [Add] M.Satoh
  integer, public, save :: &
       diag_kend_nlayer(DIAG_VMAX_DEF)       ! 11/08/16 [Add] M.Satoh
  !
  public :: diagvar_setup
  public :: diagvar_get
  public :: diagvar_set
  public :: diagvar_get_in
  public :: diagvar_set_in
  public :: diagvar_get_in_region ! 2010.5.5 M.Satoh
  public :: diagvar_set_in_region ! 2010.5.5 M.Satoh
  public :: diagvar_get_in_1layer ! 2010.4.26 M.Satoh
  public :: diagvar_set_in_1layer ! 2010.4.26 M.Satoh
  public :: diagvar_get_in_1layer_region ! 2010.5.5 M.Satoh
  public :: diagvar_set_in_1layer_region ! 2010.5.5 M.Satoh
  public :: diagvar_get_in_1layer_k ! 2011/03/02 NEC
  public :: diagvar_get_in_nlayer_region ! 2011.8.16 M.Satoh
  public :: diagvar_set_in_nlayer_region ! 2011.8.16 M.Satoh
  public :: diagvar_comm
  public :: diagvar_checkvalues
  public :: diagvar_restart_output

  character(LEN=ADM_MAXFNAME), private, save :: output_basename = '' ! [add] H.Yashiro 20120512
  ! [Add] 2012/06/07 T.Seiki, for multi-job system with LEGACY formatted data
  character(ADM_MAXFNAME), private, save  :: output_basename_CBMFX = ''
  character(ADM_MAXFNAME), private, save  :: output_basename_MP    = ''
  character(ADM_MAXFNAME), private, save  :: output_basename_QV_TB_TEND = ''
  character(ADM_MAXFNAME), private, save  :: output_basename_EVAP_SFC = '' 
  character(ADM_MAXFNAME), private, save  :: output_basename_SH_FLUX_SFC = '' 
  character(ADM_MAXFNAME), private, save  :: output_basename_ROUGHNESS_SEA = ''
  character(ADM_MAXFNAME), private, save  :: output_basename_CBMFX_CHIKIRA = ''
  
  character(ADM_MAXFNAME), private, save  :: CBMFX_fname = 'NONE'
  character(ADM_MAXFNAME), private, save  :: TB_fname    = 'NONE' ! 07/07/05 A.T.Noda
  character(ADM_MAXFNAME), private, save  :: MP_fname    = 'NONE' ! 07/12/05  T.Mitsui
  character(ADM_MAXFNAME), private, save  :: QV_TB_TEND_fname = 'NONE' ! 11/08/16 M.Satoh
!!$  character(ADM_MAXFNAME), private, save  :: RHOGQV_CONV_TB_fname = 'NONE' ! 10/05/06  M.Satoh
  character(ADM_MAXFNAME), private, save  :: EVAP_SFC_fname = 'NONE' ! 10/05/22  M.Satoh
  character(ADM_MAXFNAME), private, save  :: SH_FLUX_SFC_fname = 'NONE' ! 11/08/16b  M.Satoh
  character(ADM_MAXFNAME), private, save  :: ROUGHNESS_SEA_fname = 'NONE' ! 10/04/28 M.Satoh
  character(ADM_MAXFNAME), private, save  :: CBMFX_CHIKIRA_fname = 'NONE' ! 11/08/16 M.Satoh
  
  logical, private, save :: input_direct_access = .false.
  logical, private, save :: output_direct_access = .false.
  
  character(LEN=ADM_MAXFNAME), private, save :: input_io_mode  = 'LEGACY' ! [add] H.Yashiro 20110819
  character(LEN=ADM_MAXFNAME), private, save :: output_io_mode = 'LEGACY' ! [add] H.Yashiro 20110819
  character(LEN=ADM_MAXFNAME), private, save :: restart_layername = ''    ! [add] H.Yashiro 20110826

  integer :: NCTP = 14 !! No. of cloud types for CHIKIRA scheme: 11/08/16 M.Satoh
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine diagvar_setup
    !
    use mod_misc, only : &
         MISC_msg_nmerror
    use mod_gtl, only : &
         GTL_input_var2,  &
         GTL_input_var2_da
    use mod_runconf, only : &
         RAIN_TYPE,         &
         TB_TYPE,           &    ! 07/07/05 A.T.Noda
         RAD_CLOUD_TYPE,    &    ! 07/07/05 A.T.Noda
         AE_TYPE,           &    ! 07/11/06 T.Mitsui
         MP_TYPE,           &    ! 07/12/05 T.Mitsui
         CP_TYPE,           &    ! 10/05/05 M.Satoh
         ROUGHNESS_SEA_TYPE      ! 10/04/28 M.Satoh
    use mod_comm, only : &
         comm_var
    use mod_fio, only : & ! [add] H.Yashiro 20110819
         FIO_input
    implicit none
    !
    namelist / DIAGVARPARAM / &
         CBMFX_fname,         &
         TB_fname,            &  ! 07/07/05 A.T.Noda
         MP_fname,            &  ! 07/12/05 T.Mitsui
         QV_TB_TEND_fname,    &  ! 11/08/16 M.Satoh
!!$         RHOGQV_CONV_TB_fname, & ! 10/05/22 M.Satoh
         EVAP_SFC_fname,      &  ! 10/05/22 M.Satoh
         SH_FLUX_SFC_fname,   &  ! 11/08/16 M.Satoh
         ROUGHNESS_SEA_fname, &  ! 10/04/28 M.Satoh
         CBMFX_CHIKIRA_fname, &  ! 11/08/16 M.Satoh
         input_direct_access, &
         output_direct_access, &
         input_io_mode,       & !--- [add] H.Yashiro 20110819
         output_io_mode,      & !--- [add] H.Yashiro 20110819
         restart_layername,   & !--- [add] H.Yashiro 20110826
         output_basename, &     !--- [add] H.Yashiro 20120512
         output_basename_CBMFX,         & ! [Add] 2012/06/07 T.Seiki
         output_basename_MP,            & ! [Add] 2012/06/07 T.Seiki
         output_basename_QV_TB_TEND,    & ! [Add] 2012/06/07 T.Seiki
         output_basename_EVAP_SFC,      & ! [Add] 2012/06/07 T.Seiki
         output_basename_SH_FLUX_SFC,   & ! [Add] 2012/06/07 T.Seiki
         output_basename_ROUGHNESS_SEA, & ! [Add] 2012/06/07 T.Seiki
         output_basename_CBMFX_CHIKIRA    ! [Add] 2012/06/07 T.Seiki

    NAMELIST / NM_CP_CHIKIRA_SETUP / &
         NCTP
    !
    integer :: ierr
    integer :: iv         ! 07/12/05 [Add] T.Mitsui
    integer :: ksta, kend ! 11/07/26 [Add] M.Satoh
    !---------------------------------------------------------------------------

    !
    ! [Add] T.Mitsui 07/11/06
    I_STA        = 0
    I_END        = 0
    flag_diagset = .false.
    diag_cname(:) = ''
    diag_cname_1layer(:) = '' ! 10/04/29 M.Satoh
    diag_cname_nlayer(:) = '' ! 11/08/16 M.Satoh
    diag_knum_nlayer(:) = 0   ! 11/08/16 M.Satoh
    diag_ksta_nlayer(:) = 0   ! 11/08/16 M.Satoh
    diag_kend_nlayer(:) = 0   ! 11/08/16 M.Satoh
    ! [Add] 07/12/05 Mitsui
    if(  trim(MP_TYPE) /= 'NONE' ) then
       flag_diagset  = .true.
       I_STA         = I_END + 1
       I_UNCCN       = I_STA
       I_END         = I_STA
       I_MP_VMAX     = I_END - I_STA + 1
       diag_cname(I_UNCCN) = 'unccn'
    end if
    if(  trim(CP_TYPE)/='NONE' ) then ! 10/06/10 A.T.Noda
       flag_diagset  = .true.
       I_STA         = I_END + 1
       I_CUMCLW      = I_STA
       I_GDCLW       = I_STA + 1
       I_GDCFRC      = I_STA + 2
       I_CBMFX       = I_STA + 3
       I_END         = I_STA + 3
       ![Add] 07/12/05 Mitsui
       diag_cname(I_CUMCLW)      = 'cumclw'
       diag_cname(I_GDCLW)       = 'gdclw'
       diag_cname(I_GDCFRC)      = 'gdcfrc'
       diag_cname(I_CBMFX)       = 'cbmfx'

       ! 2010/05/11 M.Satoh
       if ( trim(CP_TYPE) == 'KUO' ) then
          I_RHOGQV_CONV    = I_STA + 4
          I_END            = I_STA + 4
          diag_cname(I_RHOGQV_CONV) = 'rhogqv_conv'
       end if
       if ( trim(CP_TYPE) == 'TDK' ) then
          ! 2011.08.16 M.Satoh,bug fix: conv => tendency
          I_QV_DYN_TEND   = I_STA + 4
          I_QV_TB_TEND    = I_STA + 5
!!$          I_RHOGQV_CONV    = I_STA + 4
!!$          I_RHOGQV_CONV_TB = I_STA + 5
          I_END            = I_STA + 5
          diag_cname(I_QV_DYN_TEND) = 'qv_dyn_tend'
          diag_cname(I_QV_TB_TEND) = 'qv_tb_tend'
!!$          diag_cname(I_RHOGQV_CONV) = 'rhogqv_conv'
!!$          diag_cname(I_RHOGQV_CONV_TB) = 'rhogqv_conv_tb'
       end if

!      I_CLOUD_PARAM_VMAX = I_END - I_STA + 1 ! not used ? ! 11/08/16b M.Satoh [del]
    end if

    if(     trim(TB_TYPE)=='MY2MOIST' .or. trim(TB_TYPE)=='MYNN2' &     !==> 09/01/28 A.T.Noda
       .or. trim(TB_TYPE)=='MYNN2.5'  .or. trim(TB_TYPE)=='MYNN3' ) then
       flag_diagset  = .true.
       I_STA = I_END + 1
       if( trim(TB_TYPE)=='MY2MOIST' .or. trim(TB_TYPE)=='MYNN2' .or. trim(TB_TYPE)=='MYNN2.5')then
         if( trim(TB_TYPE)=='MY2MOIST' .or. trim(TB_TYPE)=='MYNN2' )then
           I_END         = I_END + 1
           I_QKEd        = I_END
           diag_cname(I_QKEd)   = 'qke'
         endif
         I_END         = I_END + 1
         I_TSQd        = I_END
         I_END         = I_END + 1
         I_QSQd        = I_END
         I_END         = I_END + 1
         I_COVd        = I_END
         ![Add] 07/12/05 Mitsui
!        diag_cname(I_QKEd)   = 'qke'
         diag_cname(I_TSQd)   = 'tsq'
         diag_cname(I_QSQd)   = 'qsq'
         diag_cname(I_COVd)   = 'cov'
       endif
       ![Mod] 07/12/05 Mitsui, for memory saving
       if ( trim(RAD_CLOUD_TYPE)=='PARTIAL_CLOUD') then
!         I_CFRACP      = I_STA + 4  
!         I_QCLWB       = I_STA + 5  
!         I_QCLIB       = I_STA + 6  
!         I_END         = I_STA + 6
          I_END         = I_END + 1
          I_CFRACP      = I_END
          I_END         = I_END + 1
          I_QCLWB       = I_END
          I_END         = I_END + 1
          I_QCLIB       = I_END
          diag_cname(I_CFRACP) = 'cfrac_pdf'
          diag_cname(I_QCLWB)  = 'q_clw_before'
          diag_cname(I_QCLIB)  = 'q_cli_before'
       end if
       I_TURB_VMAX      = I_END - I_STA + 1   ! renamed 090128 A.T.Noda
    end if                                    !<== 09/01/28 A.T.Noda

    if(  trim(AE_TYPE) == 'SPRINTARS' ) then
       flag_diagset  = .true.
       I_STA         = I_END + 1
       I_DFE         = I_STA
       I_END         = I_STA
       I_SPRINTARS_VMAX  = I_END - I_STA + 1       
       ![Add] 07/12/05 Mitsui
       diag_cname(I_DFE) = 'dfe'
    end if
    DIAG_VMAX = I_END

    ! 2010.04.26 M.Satoh
    I_STA_1LAYER = 0
    I_END_1LAYER = 0
    if ( ROUGHNESS_SEA_TYPE == 'YQW' ) then
       flag_diagset  = .true.
       I_STA_1LAYER = I_END_1LAYER + 1 ! 2010.5.5 M.Satoh [add]
       I_END_1LAYER = I_STA_1LAYER     ! 2010.5.5 M.Satoh [mod]
       I_ROUGHNESS_SEA = I_END_1LAYER
! 2011/08/16b M.Satoh [del]
!       I_ROUGHNESS_SEA_MAX &
!            = I_END_1LAYER - I_STA_1LAYER + 1 ! 2010.5.5 M.Satoh [add]
       diag_cname_1layer(I_ROUGHNESS_SEA) = 'z0_roughness_sea'
    end if
    if ( CP_TYPE == 'TDK' ) then ! 2010.5.22 M.Satoh [add]
       flag_diagset  = .true.
       I_STA_1LAYER  = I_END_1LAYER + 1
       I_EVAP_SFC    = I_STA_1LAYER
       I_SH_FLUX_SFC = I_STA_1LAYER + 1
       I_END_1LAYER  = I_STA_1LAYER + 1
!       I_ROUGHNESS_SEA_MAX = I_END_1LAYER - I_STA_1LAYER + 1 ! <=??
!       diag_cname_1layer(I_ROUGHNESS_SEA) = 'evap_sfc'
       diag_cname_1layer(I_EVAP_SFC) = 'evap_sfc'
       diag_cname_1layer(I_SH_FLUX_SFC) = 'sh_flux_sfc'
    end if
    DIAG_VMAX_1LAYER = I_END_1LAYER

    ! 2011.07.25 M.Satoh
    I_STA_NLAYER = 0
    I_END_NLAYER = 0

    if ( CP_TYPE == 'CHIKIRA' ) then ! 2011.8.16 M.Satoh [add]

       rewind(ADM_CTL_FID)
       read(ADM_CTL_FID, NM_CP_CHIKIRA_SETUP,iostat=ierr)
       call MISC_msg_nmerror( &
            ierr,             & !--- in
            ADM_LOG_FID,      & !--- in
            'NM_CP_CHIKIRA_SETUP',   & !--- in
            'diagvar_setup',  & !--- in
            'diagvar'         & !--- in
            )
       write(ADM_LOG_FID,NM_CP_CHIKIRA_SETUP)
       write(ADM_LOG_FID,*) '# CHIKIRA scheme: NCTP=', NCTP

       flag_diagset  = .true.
       I_STA_NLAYER = I_END_NLAYER + 1
       I_END_NLAYER = I_STA_NLAYER
       I_CBMFX_CHIKIRA = I_END_NLAYER
       diag_cname_nlayer(I_CBMFX_CHIKIRA) = 'cbmfx_chikira'
       diag_knum_nlayer(I_CBMFX_CHIKIRA) = NCTP
    end if
    DIAG_VMAX_NLAYER = I_END_NLAYER
    if ( DIAG_VMAX_NLAYER > 0 ) then
       diag_ksta_nlayer(1) = 1
       diag_kend_nlayer(1) = diag_knum_nlayer(1)
       if ( DIAG_VMAX_NLAYER > 1 ) then
          do iv=2, DIAG_VMAX_NLAYER
             diag_ksta_nlayer(iv) = diag_kend_nlayer(iv-1) + 1
             diag_kend_nlayer(iv) &
                  = diag_ksta_nlayer(iv) + diag_knum_nlayer(iv) - 1
          end do
       end if
       DIAG_KTOT_NLAYER = diag_kend_nlayer(DIAG_VMAX_NLAYER)
    end if


    ! 2010.5.5. M.Satoh
    if ( DIAG_VMAX > DIAG_VMAX_DEF ) then
       write(ADM_LOG_FID,*) &
            "DIAG_VMAX exceeds the default value", DIAG_VMAX, DIAG_VMAX_DEF 
    end if
    if ( DIAG_VMAX_1LAYER > DIAG_VMAX_DEF ) then
       write(ADM_LOG_FID,*) &
            "DIAG_VMAX_1LAYER exceeds the default value", &
            DIAG_VMAX_1LAYER, DIAG_VMAX_DEF 
    end if
    if ( DIAG_VMAX_NLAYER > DIAG_VMAX_DEF ) then ! 11/08/16 M.Satoh
       write(ADM_LOG_FID,*) &
            "DIAG_VMAX_NLAYER exceeds the default value", &
            DIAG_VMAX_NLAYER, DIAG_VMAX_DEF 
    end if

    !
    write(ADM_LOG_FID,'(a  )') 'Msg : Sub[diagvar_setup]/Mod[diagvar]'
    write(ADM_LOG_FID,*      ) "DIAG_VMAX   = ", DIAG_VMAX
    write(ADM_LOG_FID,*      ) "DIAG_VMAX_1LAYER = ", DIAG_VMAX_1LAYER ! 10/04/30 M.Satoh
    write(ADM_LOG_FID,*      ) "DIAG_VMAX_NLAYER = ", DIAG_VMAX_NLAYER ! 11/08/16 M.Satoh
    write(ADM_LOG_FID,*      ) "flag_diagset= ", flag_diagset
    ! [Add] 07/12/05 Mitsui
    do iv=1, DIAG_VMAX
       write(ADM_LOG_FID,'(a20, i5 )') trim(diag_cname(iv)),iv
    end do
    ! [Add] 10/04/30 M.Satoh
    if ( DIAG_VMAX_1LAYER > 0 ) then
       do iv=1, DIAG_VMAX_1LAYER
          write(ADM_LOG_FID,'(a20, i5 )') trim(diag_cname_1layer(iv)),iv
       end do
    end if
    ! [Add] 11/08/16 M.Satoh
    if ( DIAG_VMAX_NLAYER > 0 ) then
       do iv=1, DIAG_VMAX_NLAYER
          write(ADM_LOG_FID,'(a20, i5 )') trim(diag_cname_nlayer(iv)),iv
       end do
    end if

    ! -> [mod] H.Yashiro 20111102 change position of config file reading  
    if( flag_diagset )then
       rewind(ADM_CTL_FID)
       read(ADM_CTL_FID, DIAGVARPARAM,iostat=ierr)
       call MISC_msg_nmerror( &
            ierr,             & !--- in
            ADM_LOG_FID,      & !--- in
            'DIAGVARPARAM',   & !--- in
            'diagvar_setup',  & !--- in
            'diagvar'         & !--- in
            )
       write(ADM_LOG_FID,DIAGVARPARAM)
    endif

    ! -> [add] H.Yashiro 20110819
    if ( input_io_mode == 'ADVANCED' ) then
       write(ADM_LOG_FID,*) '*** io_mode for restart,input : ',trim(input_io_mode)
       input_direct_access = .true.
    elseif( input_io_mode == 'LEGACY' ) then
       write(ADM_LOG_FID,*) '*** io_mode for restart,input : ',trim(input_io_mode)
    else
       write(ADM_LOG_FID,*) 'xxx Invalid input_io_mode!',trim(input_io_mode)
       call ADM_proc_stop
    endif
    if ( output_io_mode == 'ADVANCED' ) then
       write(ADM_LOG_FID,*) '*** io_mode for restart,output: ',trim(output_io_mode)
       output_direct_access = .true.
    elseif( output_io_mode == 'LEGACY' ) then
       write(ADM_LOG_FID,*) '*** io_mode for restart,output: ',trim(output_io_mode)
    else
       write(ADM_LOG_FID,*) 'xxx Invalid output_io_mode!',trim(output_io_mode)
       call ADM_proc_stop
    endif
    ! <- [add] H.Yashiro 20110819

    ! [Mod] T.Mitsui 07/11/06, bug fix for SPRINTARS
    if( flag_diagset )then
       ! -> [mod] H.Yashiro 20111102 change position of config file reading  

       allocate(diagvar(     &
            ADM_gall,        &
            ADM_kall,        &
            ADM_lall,        &
            DIAG_VMAX       &
            ))
       allocate(diagvar_pl(   &
            ADM_GALL_PL,     &
            ADM_kall,        &
            ADM_LALL_PL,     &
            DIAG_VMAX       &
            ))
       diagvar = 0.0D0
       diagvar_pl = 0.0D0
       
       ! [Add] 10/04/26 M.Satoh
       if ( DIAG_VMAX_1LAYER > 0 ) then 
          allocate(diagvar1(    &
               ADM_gall,        &
               ADM_KNONE,       &
               ADM_lall,        &
               DIAG_VMAX_1LAYER &
               ))
          allocate(diagvar1_pl( &
               ADM_GALL_PL,     &
               ADM_KNONE,       &
               ADM_LALL_PL,     &
               DIAG_VMAX_1LAYER &
               ))
          diagvar1 = -999.0D0
          diagvar1_pl = -999.0D0
       end if
          
       ! [Add] 11/08/16 M.Satoh
       if ( DIAG_VMAX_NLAYER > 0 ) then 
          allocate(diagvarn(    &
               ADM_gall,        &
               DIAG_KTOT_NLAYER, &
               ADM_lall,        &
               1                &
               ))
          allocate(diagvarn_pl( &
               ADM_GALL_PL,     &
               DIAG_KTOT_NLAYER, &
               ADM_LALL_PL,     &
               1                &
               ))
          diagvarn(:,:,:,:) = -999.0D0
          diagvarn_pl(:,:,:,:) = -999.0D0
       end if
       
       ! [Add] 07/12/05 Mitsui
       if(trim(MP_TYPE) /= 'NONE') then
          if(trim(MP_fname)=='NONE') then
             diagvar(:,:,:,I_UNCCN)     = 100.d6 ! 100 [/cm3]
             diagvar_pl (:,:,:,I_UNCCN) = 100.d6 ! 100 [/cm3]
          else
             ! -> [add] H.Yashiro 20110826
             if ( input_io_mode == 'ADVANCED' ) then

                call FIO_input( diagvar(:,:,:,I_UNCCN),MP_fname,'UNCCN', &
                                restart_layername,1,ADM_kall,1           )

             elseif( input_io_mode == 'LEGACY' ) then
             ! <- [add] H.Yashiro 20110826

                if(input_direct_access) then
                   call GTL_input_var2_da( trim(MP_fname),    &
                        diagvar(:,:,:,I_UNCCN), 1, ADM_kall,  &
                        recnum=1, input_size=8 )
                else
                   call GTL_input_var2( trim(MP_fname),                    &
                        diagvar(:,:,:,I_UNCCN), diagvar_pl(:,:,:,I_UNCCN), &
                       1, ADM_kall )
                endif !--- direct/sequencial

             endif !--- io_mode
          endif
       endif
       
       if( trim(CP_TYPE)=='AS'.or.trim(CP_TYPE)=='PAS' ) then ! 10/06/10 A.T.Noda
          if(CBMFX_fname=='NONE') then
             !--- nothing
          else
             ! -> [add] H.Yashiro 20110826
             if ( input_io_mode == 'ADVANCED' ) then

                call FIO_input( diagvar(:,:,:,I_CBMFX),CBMFX_fname,'CBMFX', &
                                restart_layername,1,ADM_kall,1           )

             elseif( input_io_mode == 'LEGACY' ) then
             ! <- [add] H.Yashiro 20110826

                if(input_direct_access) then
                   call GTL_input_var2_da( trim(CBMFX_fname), &
                        diagvar(:,:,:,I_CBMFX), 1, ADM_kall,  &
                        recnum=1, input_size=8 )
                else
                   call GTL_input_var2( trim(CBMFX_fname),                 &
                        diagvar(:,:,:,I_CBMFX), diagvar_pl(:,:,:,I_CBMFX), &
                        1, ADM_kall )
                endif !--- direct/sequencial

             endif !--- io_mode
          end if
       end if
       
       ! 2011/08/16 M.Satoh
       if ( trim(CP_TYPE) == 'CHIKIRA' ) then
          if(CBMFX_CHIKIRA_fname=='NONE') then
             diagvarn(:,ksta:kend,:,1) = 0.d0
             diagvarn_pl(:,ksta:kend,:,1) = 0.d0
          else
             ksta = diag_ksta_nlayer(I_CBMFX_CHIKIRA)
             kend = diag_kend_nlayer(I_CBMFX_CHIKIRA)
             ! 11/08/17 M.Satoh bug fix: ksta-kend+1 => kend-ksta+1

             ! -> [add] H.Yashiro 20110826
             if ( input_io_mode == 'ADVANCED' ) then

                call FIO_input( diagvarn(:,ksta:kend,:,1),CBMFX_CHIKIRA_fname, &
                                'CBMFX_CHIKIRA',                               &
                                'GCCKR',1,ADM_kall,1                           )

             elseif( input_io_mode == 'LEGACY' ) then
             ! <- [add] H.Yashiro 20110826

                if(input_direct_access) then
                   call GTL_input_var2_da( trim(CBMFX_CHIKIRA_fname), &
                        diagvarn(:,ksta:kend,:,1), 1, kend-ksta+1,  &
                        recnum=1, input_size=8 )
                else
                   call GTL_input_var2( trim(CBMFX_CHIKIRA_fname), &
                        diagvarn(:,ksta:kend,:,1), &
                        diagvarn_pl(:,ksta:kend,:,1), &
                        1, kend-ksta+1 )
                endif !--- direct/sequencial

             endif !--- io_mode
          end if
       end if
       
       ! 2010/05/11 M.Satoh
       ! 2011/08/16 M.Satoh, change names
       if ( trim(CP_TYPE) == 'TDK' ) then
          if(QV_TB_TEND_fname=='NONE') then
             diagvar(:,:,:,I_QV_TB_TEND) = 0.d0
             diagvar_pl(:,:,:,I_QV_TB_TEND) = 0.d0 ! 11/08/16 M.Satoh
          else
             ! -> [add] H.Yashiro 20110826
             if ( input_io_mode == 'ADVANCED' ) then

                call FIO_input( diagvar(:,:,:,I_QV_TB_TEND),QV_TB_TEND_fname, &
                                'QV_TB_TEND',                                 &
                                restart_layername,1,ADM_kall,1                )

             elseif( input_io_mode == 'LEGACY' ) then
             ! <- [add] H.Yashiro 20110826

                if(input_direct_access) then
                   call GTL_input_var2_da( trim(QV_TB_TEND_fname), &
                        diagvar(:,:,:,I_QV_TB_TEND), 1, ADM_kall,  &
                        recnum=1, input_size=8 )
                else
                   call GTL_input_var2( trim(QV_TB_TEND_fname),    &
                        diagvar(:,:,:,I_QV_TB_TEND), &
                        diagvar_pl(:,:,:,I_QV_TB_TEND), &
                        1, ADM_kall )
                endif !--- direct/sequencial

             endif !--- io_mode
          end if
          ! 2010.5.22 M.Satoh
          if( EVAP_SFC_fname == 'NONE' ) then
             ! set initial value of evap_sfc
             diagvar1(:,:,:,I_EVAP_SFC) = 0.0d0
             diagvar1_pl(:,:,:,I_EVAP_SFC) = 0.0d0 ! 11/08/16 M.Satoh
          else
             ! -> [add] H.Yashiro 20110826
             if ( input_io_mode == 'ADVANCED' ) then

                call FIO_input( diagvar(:,:,:,I_EVAP_SFC),EVAP_SFC_fname, &
                                'EVAP_SFC', 'ZSSFC1',1,1,1                )

             elseif( input_io_mode == 'LEGACY' ) then
             ! <- [add] H.Yashiro 20110826

                if(input_direct_access) then
                   call GTL_input_var2_da( trim(EVAP_SFC_fname), &
                        diagvar1(:,:,:,I_EVAP_SFC), 1, ADM_KNONE,  &
                        recnum=1, input_size=8 )
                else
                   call GTL_input_var2( trim(EVAP_SFC_fname), &
                        diagvar1(:,:,:,I_EVAP_SFC), &
                        diagvar1_pl(:,:,:,I_EVAP_SFC), &
                        1, ADM_KNONE )
                endif !--- direct/sequencial

             endif !--- io_mode

             write(ADM_LOG_FID,*) 'restart in: evap_sfc', &
                  maxval(diagvar1(:,:,:,I_EVAP_SFC)), &
                  minval(diagvar1(:,:,:,I_EVAP_SFC))
          end if

          ! 2011/08/16b M.Satoh
          if( SH_FLUX_SFC_fname == 'NONE' ) then
             ! set initial value of sh_flux_sfc
             diagvar1(:,:,:,I_SH_FLUX_SFC) = 0.0d0
             diagvar1_pl(:,:,:,I_SH_FLUX_SFC) = 0.0d0
          else
             ! -> [add] H.Yashiro 20110826
             if ( input_io_mode == 'ADVANCED' ) then

                call FIO_input( diagvar(:,:,:,I_SH_FLUX_SFC),SH_FLUX_SFC_fname, &
                                'SH_FLUX_SFC','ZSSFC1',1,1,1                    )

             elseif( input_io_mode == 'LEGACY' ) then
             ! <- [add] H.Yashiro 20110826

                if(input_direct_access) then
                   call GTL_input_var2_da( trim(SH_FLUX_SFC_fname), &
                        diagvar1(:,:,:,I_SH_FLUX_SFC), 1, ADM_KNONE,  &
                        recnum=1, input_size=8 )
                else
                   call GTL_input_var2( trim(SH_FLUX_SFC_fname), &
                        diagvar1(:,:,:,I_SH_FLUX_SFC), &
                        diagvar1_pl(:,:,:,I_SH_FLUX_SFC), &
                        1, ADM_KNONE )
                endif !--- direct/sequencial

             endif !--- io_mode

             write(ADM_LOG_FID,*) 'restart in: sh_sfc', &
                  maxval(diagvar1(:,:,:,I_SH_FLUX_SFC)), &
                  minval(diagvar1(:,:,:,I_SH_FLUX_SFC))
          end if

       end if
       
       ! -> [add&mod] H.Yashiro 20110826
       !==> 09/01/28 A.T.Noda
       if(TB_fname=='NONE') then
          !--- nothing
       else
          if (      trim(TB_TYPE)=='MY2MOIST' &
               .or. trim(TB_TYPE)=='MYNN2'    &
               .or. trim(TB_TYPE)=='MYNN2.5'  &
               .or. trim(TB_TYPE)=='MYNN3'    ) then
             !
             if (      trim(TB_TYPE)=='MY2MOIST' &
                  .or. trim(TB_TYPE)=='MYNN2'    &
                  .or. trim(TB_TYPE)=='MYNN2.5'  ) then
                !
                if (      trim(TB_TYPE)=='MY2MOIST' &
                     .or. trim(TB_TYPE)=='MYNN2'    ) then
                   !
                   if ( input_io_mode == 'ADVANCED' ) then
                      !
                      call FIO_input( diagvar(:,:,:,I_QKEd),TB_fname,       &
                                      'qked',restart_layername,1,ADM_kall,1 )

                   elseif( input_io_mode == 'LEGACY' ) then
                      !
                      if (input_direct_access) then
                         call GTL_input_var2_da( trim(TB_fname), &
                              diagvar(:,:,:,I_QKEd), 1, ADM_kall,  &
                              recnum=1, input_size=8 )
                      else
                         call GTL_input_var2( trim(TB_fname),                  &
                              diagvar(:,:,:,I_QKEd), diagvar_pl(:,:,:,I_QKEd), &
                              1, ADM_kall )
                      endif !--- direct/sequencial

                   endif  !--- io_mode

                   write(ADM_LOG_FID,*) 'restart in: qked', &
                        maxval(diagvar(:,:,:,I_QKEd)), minval(diagvar(:,:,:,I_QKEd))
                endif !--- MY2M, MYNN2

                !
                if ( input_io_mode == 'ADVANCED' ) then
                   !
                   call FIO_input( diagvar(:,:,:,I_TSQd),TB_fname,       &
                                   'tsqd',restart_layername,1,ADM_kall,1 )
                   call FIO_input( diagvar(:,:,:,I_QSQd),TB_fname,       &
                                   'qsqd',restart_layername,1,ADM_kall,1 )
                   call FIO_input( diagvar(:,:,:,I_COVd),TB_fname,       &
                                   'covd',restart_layername,1,ADM_kall,1 )

                elseif( input_io_mode == 'LEGACY' ) then
                   !
                   if (input_direct_access) then
                      call GTL_input_var2_da( trim(TB_fname), &
                           diagvar(:,:,:,I_TSQd), 1, ADM_kall,  &
                           recnum=2, input_size=8 )
                      call GTL_input_var2_da( trim(TB_fname), &
                           diagvar(:,:,:,I_QSQd), 1, ADM_kall,  &
                           recnum=3, input_size=8 )
                      call GTL_input_var2_da( trim(TB_fname), &
                           diagvar(:,:,:,I_COVd), 1, ADM_kall,  &
                           recnum=4, input_size=8 )
                   else
                      call GTL_input_var2( trim(TB_fname),                  &
                           diagvar(:,:,:,I_TSQd), diagvar_pl(:,:,:,I_TSQd), &
                           1, ADM_kall )
                      call GTL_input_var2( trim(TB_fname),                  &
                           diagvar(:,:,:,I_QSQd), diagvar_pl(:,:,:,I_QSQd), &
                           1, ADM_kall )
                      call GTL_input_var2( trim(TB_fname),                  &
                           diagvar(:,:,:,I_COVd), diagvar_pl(:,:,:,I_COVd), &
                           1, ADM_kall )
                   endif !--- direct/sequencial

                endif  !--- io_mode

                write(ADM_LOG_FID,*) 'restart in: tsqd', &
                     maxval(diagvar(:,:,:,I_TSQd)), minval(diagvar(:,:,:,I_TSQd))
                write(ADM_LOG_FID,*) 'restart in: qsqd', &
                     maxval(diagvar(:,:,:,I_QSQd)), minval(diagvar(:,:,:,I_QSQd))
                write(ADM_LOG_FID,*) 'restart in: covd', &
                     maxval(diagvar(:,:,:,I_COVd)), minval(diagvar(:,:,:,I_COVd))
             endif !--- MY2M, MYNN2, MYNN2.5

             if ( RAD_CLOUD_TYPE=='PARTIAL_CLOUD' )then
                !
                if ( input_io_mode == 'ADVANCED' ) then
                   !
                   call FIO_input( diagvar(:,:,:,I_QCLWB),TB_fname,              &
                                   'q_clw_before',restart_layername,1,ADM_kall,1 )
                   call FIO_input( diagvar(:,:,:,I_QCLIB),TB_fname,              &
                                   'q_cli_before',restart_layername,1,ADM_kall,1 )
                   call FIO_input( diagvar(:,:,:,I_CFRACP),TB_fname,          &
                                   'cfrac_pdf',restart_layername,1,ADM_kall,1 )

                elseif( input_io_mode == 'LEGACY' ) then
                   !
                   if (input_direct_access) then
                   call GTL_input_var2_da( trim(TB_fname), &
                        diagvar(:,:,:,I_QCLWB ), 1, ADM_kall,  &
                        recnum=5, input_size=8 )
                   call GTL_input_var2_da( trim(TB_fname), &
                        diagvar(:,:,:,I_QCLIB ), 1, ADM_kall,  &
                        recnum=6, input_size=8 )
                   call GTL_input_var2_da( trim(TB_fname), &
                        diagvar(:,:,:,I_CFRACP ), 1, ADM_kall,  &
                        recnum=7, input_size=8 )
                   else
                   call GTL_input_var2( trim(TB_fname),                  &
                        diagvar(:,:,:,I_QCLWB), diagvar_pl(:,:,:,I_QCLWB), &
                        1, ADM_kall )
                   call GTL_input_var2( trim(TB_fname),                  &
                        diagvar(:,:,:,I_QCLIB), diagvar_pl(:,:,:,I_QCLIB), &
                        1, ADM_kall )
                   call GTL_input_var2( trim(TB_fname),                  &
                        diagvar(:,:,:,I_CFRACP), diagvar_pl(:,:,:,I_CFRACP), &
                        1, ADM_kall )
                   endif !--- direct/sequencial

                endif  !--- io_mode

                write(ADM_LOG_FID,*) 'restart in: q_clw_before', &
                     maxval(diagvar(:,:,:,I_QCLWB)), minval(diagvar(:,:,:,I_QCLWB))
                write(ADM_LOG_FID,*) 'restart in: q_cli_before', &
                     maxval(diagvar(:,:,:,I_QCLIB)), minval(diagvar(:,:,:,I_QCLIB))
                write(ADM_LOG_FID,*) 'restart in: cfrac_pdf', &
                     maxval(diagvar(:,:,:,I_CFRACP)), minval(diagvar(:,:,:,I_CFRACP))
             endif !--- PARTIAL_CLOUD

          endif !--- MY2M, MYNN2, MYNN2.5
       endif
       !<=== 09/01/28 A.T.Noda
       ! <- [add&mod] H.Yashiro 20110826
       
       ! 10/04/28 M.Satoh
       if( ROUGHNESS_SEA_TYPE == 'YQW' ) then 
          if( ROUGHNESS_SEA_fname == 'NONE' ) then 
             ! set initial value in roughness_sea_init
             diagvar1(:,:,:,I_ROUGHNESS_SEA) = -999.0d0
             diagvar1_pl(:,:,:,I_ROUGHNESS_SEA) = -999.0d0 ! 11/08/16 M.Satoh
          else
             ! -> [add] H.Yashiro 20110826
             if ( input_io_mode == 'ADVANCED' ) then
                ! -> [mod] Y.Yamada 20111125
                !                call FIO_input( diagvar(:,:,:,I_ROUGHNESS_SEA),ROUGHNESS_SEA_fname, &
                !                                'ROUGHNESS_SEA','ZSSFC1',1,1,1                      )
                call FIO_input( diagvar1(:,:,:,I_ROUGHNESS_SEA),ROUGHNESS_SEA_fname, &
                     'ROUGHNESS_SEA','ZSSFC1',1,1,1                      )
                ! <- [mod] Y.Yamada 20111125
             elseif( input_io_mode == 'LEGACY' ) then
             ! <- [add] H.Yashiro 20110826

                if(input_direct_access) then
                   call GTL_input_var2_da( trim(ROUGHNESS_SEA_fname), &
                        diagvar1(:,:,:,I_ROUGHNESS_SEA), 1, ADM_KNONE,  &
                        recnum=1, input_size=8 )
                else
                   call GTL_input_var2( trim(ROUGHNESS_SEA_fname), &
                        diagvar1(:,:,:,I_ROUGHNESS_SEA), &
                        diagvar1_pl(:,:,:,I_ROUGHNESS_SEA), &
                        1, ADM_KNONE )
                endif !--- direct/sequencial

             endif !--- io_mode

             write(ADM_LOG_FID,*) 'restart in: z0_roughness_sea', &
                  maxval(diagvar1(:,:,:,I_ROUGHNESS_SEA)), &
                  minval(diagvar1(:,:,:,I_ROUGHNESS_SEA))
          end if
       end if
       ! <= 10/04/28 M.Satoh
       
       ! [Add] T.Mitsui 07/11/06, call at once 
       if(input_direct_access) then
          call comm_var(            &
               diagvar,diagvar_pl,  &
               ADM_kall,            &
               DIAG_VMAX,           &
               comm_type = 2,       &
               NSval_fix=.true.     &
               )
          if ( DIAG_VMAX_1LAYER > 0 ) then ! [add] 2010.4.28 M.Satoh 
             call comm_var(            &
                  diagvar1,diagvar1_pl,&
                  ADM_KNONE,           &
                  DIAG_VMAX_1LAYER,    &
                  comm_type = 2,       &
                  NSval_fix=.true.     &
                  )
          end if
          if ( DIAG_VMAX_NLAYER > 0 ) then ! [add] 2011.8.16 M.Satoh 
             call comm_var(            &
                  diagvarn,diagvarn_pl,&
                  DIAG_KTOT_NLAYER,    &
                  1,                   &
                  comm_type = 2,       &
                  NSval_fix=.true.     &
                  )
          end if
       end if
       
    endif ! flag_diagset

   return
  end subroutine diagvar_setup

  !-----------------------------------------------------------------------------
  subroutine diagvar_get(   &
     dv, dv_pl,          &  !--- OUT : surface variable
     vid                 &  !--- IN  : variable ID
     )
  !------
  !------ get prognostic variables from diag[num].
  !------
  implicit none
  real(8), intent(out) :: dv(ADM_gall,ADM_kall,ADM_lall)
  real(8), intent(out) :: dv_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
  integer:: ij,k,l
  integer, intent(in)  :: vid 
  ! [Add] 07/11/06 T.Mitsui for check
  if(.not. flag_diagset .or. vid == -999 )then 
     write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
     write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
     call ADM_proc_stop
  end if
  do l=1,ADM_lall
     do k=1,ADM_kall
        do ij=1,ADM_gall
           dv(ij,k,l) = diagvar(ij,k,l,vid)
        enddo
     enddo
  enddo


  if(ADM_prc_me==ADM_prc_pl) then
     dv_pl(1:ADM_gall_pl,:,:)   = diagvar_pl(1:ADM_gall_pl,:,:,vid)
  end if
  !
  return
  !
  end subroutine diagvar_get

  !-----------------------------------------------------------------------------
  subroutine diagvar_set(   &
     dv, dv_pl,          &  !--- IN : surface variable
     vid                 &  !--- IN : variable ID
     )
    !------
    !------ set diagnostic variables
    !------ and COMMUNICATION.
    !------ 
    implicit none
    real(8), intent(in) :: dv(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: dv_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    integer, intent(in) :: vid
    !
    integer :: i,j,suf,ij,k,l
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    ! [Add] 07/11/06 T.Mitsui for check
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if

  do l=1,ADM_lall
     do k=1,ADM_kall
        do ij=1,ADM_gall
           diagvar(ij,k,l,vid)   = dv(ij,k,l)
        enddo
     enddo
  enddo

    if(ADM_prc_me==ADM_prc_pl) then
       diagvar_pl(1:ADM_gall_pl,:,:,vid)   = dv_pl(1:ADM_gall_pl,:,:)
    end if
    !
    diagvar(suf(ADM_gall_1d,1),:,:,vid) = diagvar(suf(ADM_gmax+1,ADM_gmin),:,:,vid)
    diagvar(suf(1,ADM_gall_1d),:,:,vid) = diagvar(suf(ADM_gmin,ADM_gmax+1),:,:,vid)
    !
    return
    !
  end subroutine diagvar_set

  !-----------------------------------------------------------------------------
  subroutine diagvar_get_in(&
       sv,                 &  !--- OUT : surface variable
       vid                 &  !--- IN : variable ID
       )
    !------
    !------ get diagnostic variables
    !------
    implicit none
    real(8), intent(out) :: sv(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn,k
    !
    ! [Add] 07/11/06 T.Mitsui for check
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do l=1, ADM_lall
       do k = 1, ADM_kall
          do n=1, ADM_IopJop_nmax
             nn = ADM_IopJop(n,ADM_GIoJo)
             sv(n,k,l)  = diagvar(nn,k,l,vid)
          end do
       end do
    end do
    !
    return
    !
  end subroutine diagvar_get_in

  !-----------------------------------------------------------------------------
  subroutine diagvar_set_in( &
       sv,                  &  !--- IN : surface variable
       vid                  &  !--- IN : variable ID
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------ 
    implicit none
    real(8), intent(in) :: sv(ADM_IopJop_nmax,ADM_kall,ADM_lall)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn,k
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    ! [Add] 07/11/06 T.Mitsui for check
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do l=1, ADM_lall
       do k = 1, ADM_kall
          do n=1, ADM_IopJop_nmax
             nn = ADM_IopJop(n,ADM_GIoJo)
             diagvar(nn,k,l,vid) = sv(n,k,l)
          end do
       end do
    end do
    !
    diagvar(suf(ADM_gall_1d,1),:,:,vid) = diagvar(suf(ADM_gmax+1,ADM_gmin),:,:,vid)
    diagvar(suf(1,ADM_gall_1d),:,:,vid) = diagvar(suf(ADM_gmin,ADM_gmax+1),:,:,vid)
    !
  end subroutine diagvar_set_in

  !-----------------------------------------------------------------------------
  ! 2010.5.5 M.Satoh [add]
  subroutine diagvar_get_in_region (&
       sv,                 &  !--- OUT : surface variable
       vid,                &  !--- IN : variable ID
       l_region            &  !--- IN
       )
    !------
    !------ get diagnostic variables
    !------
    implicit none
    real(8), intent(out) :: sv(ADM_IopJop_nmax,ADM_kall)
    integer, intent(in) :: vid
    integer, intent(in) :: l_region
    !
    integer :: l,n,nn,k
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do k = 1, ADM_kall
       do n=1, ADM_IopJop_nmax
          nn = ADM_IopJop(n,ADM_GIoJo)
          sv(n,k)  = diagvar(nn,k,l_region,vid)
       end do
    end do
    !
    return
    !
  end subroutine diagvar_get_in_region

  !-----------------------------------------------------------------------------
  ! 2010.5.5 M.Satoh [add]
  subroutine diagvar_set_in_region ( &
       sv,                  &  !--- IN : surface variable
       vid,                 &  !--- IN : variable ID
       l_region             &  !--- IN
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------ 
    implicit none
    real(8), intent(in) :: sv(ADM_IopJop_nmax,ADM_kall)
    integer, intent(in) :: vid
    integer, intent(in) :: l_region
    !
    integer :: l,n,nn,k
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do k = 1, ADM_kall
       do n=1, ADM_IopJop_nmax
          nn = ADM_IopJop(n,ADM_GIoJo)
          ![fix] 12/03/28 T.Seiki
!!$       diagvar(nn,k,l,vid) = sv(n,k)
          diagvar(nn,k,l_region,vid) = sv(n,k)
       end do
    end do
    !
    diagvar(suf(ADM_gall_1d,1),:,l_region,vid) &
         = diagvar(suf(ADM_gmax+1,ADM_gmin),:,l_region,vid)
    diagvar(suf(1,ADM_gall_1d),:,l_region,vid) &
         = diagvar(suf(ADM_gmin,ADM_gmax+1),:,l_region,vid)
    !
  end subroutine diagvar_set_in_region

  !-----------------------------------------------------------------------------
  ! 2010.4.26 M.Satoh
  subroutine diagvar_get_in_1layer ( &
       sv,                 &  !--- OUT : surface variable
       vid                 &  !--- IN : variable ID
       )
    !------
    !------ get diagnostic variables
    !------
    implicit none
    real(8), intent(out) :: sv(ADM_IopJop_nmax,ADM_lall)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do l=1, ADM_lall
       do n=1, ADM_IopJop_nmax
          nn = ADM_IopJop(n,ADM_GIoJo)
          sv(n,l)  = diagvar1(nn,ADM_KNONE,l,vid)
       end do
    end do
    !
    return
    !
  end subroutine diagvar_get_in_1layer

  !-----------------------------------------------------------------------------
  ! 2011/03/02 NEC   
  subroutine diagvar_get_in_1layer_k ( &
       sv,                 &  !--- OUT : surface variable
       ksize ,             &  !--- IN : 
       k ,                 &  !--- IN :
       vid                 &  !--- IN : variable ID
       )
    !------
    !------ get diagnostic variables
    !------
    implicit none
    integer, intent(in) :: ksize,k
    real(8), intent(out) :: sv(ADM_IopJop_nmax,ksize,ADM_lall)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do l=1, ADM_lall
       do n=1, ADM_IopJop_nmax
          nn = ADM_IopJop(n,ADM_GIoJo)
          sv(n,k,l)  = diagvar1(nn,ADM_KNONE,l,vid)
       end do
    end do
    !
    return
    !
  end subroutine diagvar_get_in_1layer_k
 
  !-----------------------------------------------------------------------------
  subroutine diagvar_set_in_1layer ( &
       sv,                  &  !--- IN : surface variable
       vid                  &  !--- IN : variable ID
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------ 
    implicit none
    real(8), intent(in) :: sv(ADM_IopJop_nmax,ADM_lall)
    integer, intent(in) :: vid
    !
    integer :: l,n,nn
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do l=1, ADM_lall
       do n=1, ADM_IopJop_nmax
          nn = ADM_IopJop(n,ADM_GIoJo)
          diagvar1(nn,ADM_KNONE,l,vid) = sv(n,l)
       end do
    end do
    !
    diagvar1(suf(ADM_gall_1d,1),:,:,vid) &
         = diagvar1(suf(ADM_gmax+1,ADM_gmin),:,:,vid)
    diagvar1(suf(1,ADM_gall_1d),:,:,vid) &
         = diagvar1(suf(ADM_gmin,ADM_gmax+1),:,:,vid)
    !
  end subroutine diagvar_set_in_1layer

  !-----------------------------------------------------------------------------
  subroutine diagvar_get_in_1layer_region ( &
       sv,                 &  !--- OUT : surface variable
       vid,                &  !--- IN : variable ID
       l_region            &  !--- IN
       )
    !------
    !------ get diagnostic variables
    !------
    implicit none
    real(8), intent(out) :: sv(ADM_IopJop_nmax)
    integer, intent(in) :: vid
    integer, intent(in) :: l_region
    !
    integer :: n,nn
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do n=1, ADM_IopJop_nmax
       nn = ADM_IopJop(n,ADM_GIoJo)
       sv(n)  = diagvar1(nn,ADM_KNONE,l_region,vid)
    end do
    !
    return
    !
  end subroutine diagvar_get_in_1layer_region

  !-----------------------------------------------------------------------------
  subroutine diagvar_set_in_1layer_region ( &
       sv,                  &  !--- IN : surface variable
       vid,                 &  !--- IN : variable ID
       l_region             &  !--- IN
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------ 
    implicit none
    real(8), intent(in) :: sv(ADM_IopJop_nmax)
    integer, intent(in) :: vid
    integer, intent(in) :: l_region
    !
    integer :: n,nn
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    do n=1, ADM_IopJop_nmax
       nn = ADM_IopJop(n,ADM_GIoJo)
       diagvar1(nn,ADM_KNONE,l_region,vid) = sv(n)
    end do
    !
!!$    diagvar1(suf(ADM_gall_1d,1),:,:,vid) &
!!$         = diagvar1(suf(ADM_gmax+1,ADM_gmin),:,:,vid)
!!$    diagvar1(suf(1,ADM_gall_1d),:,:,vid) &
!!$         = diagvar1(suf(ADM_gmin,ADM_gmax+1),:,:,vid)
    ! 2010.5.5. M. Satoh [mod]
    diagvar1(suf(ADM_gall_1d,1),:,l_region,vid) &
         = diagvar1(suf(ADM_gmax+1,ADM_gmin),:,l_region,vid)
    diagvar1(suf(1,ADM_gall_1d),:,l_region,vid) &
         = diagvar1(suf(ADM_gmin,ADM_gmax+1),:,l_region,vid)
    !
  end subroutine diagvar_set_in_1layer_region

  !-----------------------------------------------------------------------------
  ! 2011.8.16 M.Satoh [add]
  subroutine diagvar_get_in_nlayer_region (&
       sv,                 &  !--- OUT : variable
       vid,                &  !--- IN : variable ID
       knum,               &  !--- IN : number of layers
       l_region            )  !--- IN
    implicit none

    integer, intent(in) :: knum
    real(8), intent(out) :: sv(ADM_IopJop_nmax,knum)
    integer, intent(in) :: vid
    integer, intent(in) :: l_region
    !
    integer :: l,n,nn,k
    integer :: ksta
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    if ( DIAG_VMAX_NLAYER >= vid &
         .and. diag_knum_nlayer(vid) == knum ) then
       ksta = diag_ksta_nlayer(vid)
       do k = 1, knum
          do n=1, ADM_IopJop_nmax
             nn = ADM_IopJop(n,ADM_GIoJo)
             sv(n,k)  = diagvarn(nn,ksta+k-1,l_region,1)
          end do
       end do
    else
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR illegal diag_vmax_nlayer", &
            DIAG_VMAX_NLAYER
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if

    return
    !
  end subroutine diagvar_get_in_nlayer_region

  !-----------------------------------------------------------------------------
  ! 2011.8.16 M.Satoh [add]
  subroutine diagvar_set_in_nlayer_region ( &
       sv,                  &  !--- IN : variable
       vid,                 &  !--- IN : variable ID
       knum,                &  !--- IN : number of layers
       l_region             &  !--- IN
       )
    !------
    !------ set prognostic variables to diag[num]
    !------ and COMMUNICATION.
    !------ 
    implicit none
    real(8), intent(in) :: sv(ADM_IopJop_nmax,ADM_kall)
    integer, intent(in) :: vid
    integer, intent(in) :: knum
    integer, intent(in) :: l_region
    !
    integer :: l,n,nn,k
    integer :: ksta, kend
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    if(.not. flag_diagset .or. vid == -999 )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if
    if ( DIAG_VMAX_NLAYER >= vid &
         .and. diag_knum_nlayer(vid) == knum ) then
       ksta = diag_ksta_nlayer(vid)
       kend = diag_kend_nlayer(vid)
!      do k = 1, ADM_kall  ! 110830 A.Noda
       do k = 1, knum
          do n=1, ADM_IopJop_nmax
             nn = ADM_IopJop(n,ADM_GIoJo)
             diagvarn(nn,ksta+k-1,l_region,1) = sv(n,k)
          end do
       end do
       diagvarn(suf(ADM_gall_1d,1),ksta:kend,l_region,1) &
            = diagvar(suf(ADM_gmax+1,ADM_gmin),ksta:kend,l_region,1)
       diagvarn(suf(1,ADM_gall_1d),ksta:kend,l_region,1) &
            = diagvar(suf(ADM_gmin,ADM_gmax+1),ksta:kend,l_region,1)
    else
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR illegal diag_vmax_nlayer", &
            DIAG_VMAX_NLAYER
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if

    !
  end subroutine diagvar_set_in_nlayer_region

  !-----------------------------------------------------------------------------
  subroutine diagvar_comm( &
       comm_type          & !--- IN : communication type
       )
    !
    !--- comm_type : 1 ( region -> pole )
    !---           : 2 ( region -> pole -> regular communication )
    !---           : 3 ( regular communication only )
    use mod_adm, only :     &
         ADM_LOG_FID,       &
         ADM_KNONE,         &
         ADM_lall,          &
         ADM_proc_stop,     &
         ADM_comm_run_world
    use mod_comm, only :    &
         comm_var
    !
    implicit none
    integer, intent(in) :: comm_type

    integer :: ierr

!!$    !
!!$    integer :: ireq,istat(MPI_STATUS_SIZE),ierr
!!$    integer :: l
!!$    !
!!$    real(8) :: v_npl(ADM_kall,DIAG_VMAX)
!!$    real(8) :: v_spl(ADM_kall,DIAG_VMAX)
!!$    !
!!$    integer :: rgnid
!!$    !
!!$    integer :: i,j,suf
!!$    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
!!$    !
    ! [Add] 07/11/06 T.Mitsui for check
    if(.not. flag_diagset )then 
       write(ADM_LOG_FID,*) "*** Error, DIAGVAR has not been setup."
       write(ADM_LOG_FID,*) "    Check [sub: diagvar_setup] "
       call ADM_proc_stop
    end if















    call comm_var(  &
         diagvar,   & !--- INOUT : variables
         diagvar_pl,& !--- INOUT : variables at poles
         adm_kall,  & !--- IN : number of layers
         DIAG_VMAX, & !--- IN : number of variables
         comm_type, & !--- IN : communication type
         NSval_fix=.true.&
         )















    ! 2010.4.26 M.Satoh
    if ( DIAG_VMAX_1LAYER > 0 ) then 















       call comm_var(  &
            diagvar1,         & !--- INOUT : variables
            diagvar1_pl,      & !--- INOUT : variables at poles
            adm_KNONE,        & !--- IN : number of layers
            DIAG_VMAX_1LAYER, & !--- IN : number of variables
            comm_type,        & !--- IN : communication type
            NSval_fix=.true.  &
            )















    end if

    ! 2011.8.16 M.Satoh
    if ( DIAG_VMAX_NLAYER > 0 ) then 
       call comm_var(  &
            diagvarn,         & !--- INOUT : variables
            diagvarn_pl,      & !--- INOUT : variables at poles
            DIAG_KTOT_NLAYER, & !--- IN : number of layers
            1,                & !--- IN : number of variables
            comm_type,        & !--- IN : communication type
            NSval_fix=.true.  &
            )
    end if

    return

  end subroutine diagvar_comm

  !-----------------------------------------------------------------------------
  subroutine diagvar_checkvalues ( item, var, vmin, vmax )

    real(8), intent(in) :: var(:,:)
    real(8), intent(in) :: vmax
    real(8), intent(in) :: vmin
    character(LEN=*), intent(in) :: item

    real(8) :: varmax
    real(8) :: varmin

    varmax = maxval( var )
    varmin = minval( var )

    ! 04/09/07 M.Satoh: bug fix: (varmin > vmin) => (varmin < vmin)
    if ( varmax > vmax .or. varmin < vmin ) then
       write(ADM_LOG_FID,*) &
            '@DIAGVAR_CHECKVALUES: invalid values', &
            item, &
            ' :max=', varmax, &
            ' :min=', varmin, &
            vmax, vmin
    end if

  end subroutine diagvar_checkvalues

  !-----------------------------------------------------------------------------
  subroutine diagvar_restart_output( cdate ) !08/03/10 T.Mitsui [Add] cdate
    use mod_gtl, only : &
         GTL_output_var2,&
         GTL_output_var2_da
    use mod_adm, only :     &               ! 07/07/05 A.T.Noda
         ADM_NSYS
    use mod_runconf, only : &
         RAIN_TYPE,         &
         MP_TYPE,           &               ! 07/12/05  T.Mitsui
         CP_TYPE,           &               ! 10/05/24  M.Satoh
         AE_TYPE,           &               ! 07/12/05  T.Mitsui
         TB_TYPE,           &               ! 07/07/05 A.T.Noda
         RAD_CLOUD_TYPE,    &               ! 07/07/05 A.T.Noda
         ROUGHNESS_SEA_TYPE                 ! 10/04/28 M.Satoh
    use mod_fio, only : & ! [add] H.Yashiro 20110819
         FIO_output, &
         FIO_HMID,   &
         FIO_REAL8
    use mod_time, only : &
         TIME_CTIME
    implicit none
    integer :: ksta, kend
    !
    character(len=14), intent(in) :: cdate !08/03/10 T.Mitsui [Add] 
    ! 07/12/05 [Del] T.Mitsui
!!$    character(LEN=ADM_NSYS) :: cname(10)    ! 07/07/05 A.T.Noda
    integer :: iv     ! 07/12/05 [Add] T.Mitsui
!!$    ! [Add] work T.Mitsui 10/08/03
!!$    character(len=256) :: cwork(256)
!!$    character(len=256) :: cwork_1layer(256)

    ! -> [add] H.Yashiro 20110819
    character(LEN=ADM_MAXFNAME) :: basename
    character(LEN=FIO_HMID)     :: desc = 'INITIAL/RESTART DATA of DIAGNOSTIC VARIABLES'
    ! <- [add] H.Yashiro 20110819
    !---------------------------------------------------------------------------

    !
    ! [Add] 08/05/30 T.Mitsui 
    call diagvar_comm( comm_type=2 )

!   if(RAIN_TYPE/='CLOUD_PARAM') return 
!   if(RAIN_TYPE=='CLOUD_PARAM')then        ! 07/07/05 A.T.Noda
    if(trim(CP_TYPE)=='AS'.or.trim(CP_TYPE)=='PAS')then ! 10/06/10 A.T.Noda

       ! -> [add] H.Yashiro 20110819
       if ( output_io_mode == 'ADVANCED' ) then
          if ( trim(output_basename) == "" ) then ! [add] H.Yashiro 20120512
             basename='restart_CBMFX'//trim(cdate)
          else
             basename=trim(output_basename)//trim(cdate)
          endif

          call FIO_output( diagvar(:,:,:,I_CBMFX), basename, desc, '',     &
                           'CBMFX', 'Cloud-base Mass Flux', '', 'kg/m2/s', &
                           FIO_REAL8, restart_layername, 1, ADM_kall,      &
                           1, TIME_CTIME, TIME_CTIME                       )

       elseif( output_io_mode == 'LEGACY' ) then
       ! <- [add] H.Yashiro 20110819
          ! [Add] 2012/06/07 T.Seiki
          if ( trim(output_basename_CBMFX) == "" ) then 
             basename='restart_CBMFX'//trim(cdate)
          else
             basename=trim(output_basename_CBMFX)//trim(cdate)
          endif
          
          if (output_direct_access) then    
             call GTL_output_var2_da( basename, &
                                      diagvar(:,:,:,I_CBMFX), 1, ADM_kall, &
                                      recnum=1, output_size=8 )
          else
             call GTL_output_var2( basename,        &
                                   diagvar(:,:,:,I_CBMFX), diagvar_pl(:,:,:,I_CBMFX), &
                                   1, ADM_kall )
          endif !--- direct/sequencial

       endif  !--- io_mode
    endif

    ! 07/12/05 [Add] T.Mitsui
    if(  (MP_TYPE/='NONE') .and. (AE_TYPE/='NONE') )then 

       ! -> [add] H.Yashiro 20110819
       if ( output_io_mode == 'ADVANCED' ) then
          if ( trim(output_basename) == "" ) then ! [add] H.Yashiro 20120512
             basename='restart_mp'//trim(cdate)
          else
             basename=trim(output_basename)//trim(cdate)
          endif

          call FIO_output( diagvar(:,:,:,I_UNCCN), basename, desc, '', &
                           'UNCCN', 'CCN Concentration', '', '1/m3',   &
                           FIO_REAL8, restart_layername, 1, ADM_kall,  &
                           1, TIME_CTIME, TIME_CTIME                   )

       elseif( output_io_mode == 'LEGACY' ) then
       ! <- [add] H.Yashiro 20110819
          ! [Add] 2012/06/07 T.Seiki
          if ( trim(output_basename_MP) == "" ) then 
             basename='restart_mp'//trim(cdate)
          else
             basename=trim(output_basename_MP)//trim(cdate)
          endif
          if(output_direct_access) then    
             call GTL_output_var2_da( basename, &
                  diagvar(:,:,:,I_UNCCN), 1, ADM_kall, &
                  recnum=1, output_size=8 )
          else
             call GTL_output_var2( basename,           &
                  diagvar(:,:,:,I_UNCCN), diagvar_pl(:,:,:,I_UNCCN),    &
                  1, ADM_kall )
          endif !--- direct/sequencial

       endif  !--- io_mode
    endif

    ! 2010.5.24 M.Satoh
    if ( trim(CP_TYPE) == 'TDK' ) then

       ! -> [add] H.Yashiro 20110819
       if ( output_io_mode == 'ADVANCED' ) then
          if ( trim(output_basename) == "" ) then ! [add] H.Yashiro 20120512
             basename='restart_QV_TB_TEND'//trim(cdate)
          else
             basename=trim(output_basename)//trim(cdate)
          endif

          call FIO_output( diagvar(:,:,:,I_QV_TB_TEND), basename, desc, '', &
                           'QV_TB_TEND', 'QV_TB_TEND', '', '',              &
                           FIO_REAL8, restart_layername, 1, ADM_kall,       &
                           1, TIME_CTIME, TIME_CTIME                        )

          if ( trim(output_basename) == "" ) then ! [add] H.Yashiro 20120512
             basename='restart_EVAP_SFC'//trim(cdate)
          else
             basename=trim(output_basename)//trim(cdate)
          endif

          call FIO_output( diagvar(:,:,:,I_EVAP_SFC), basename, desc, '', &
                           'EVAP_SFC', 'EVAP_SFC', '', '',                &
                           FIO_REAL8, 'ZSSFC1', 1, 1,                     &
                           1, TIME_CTIME, TIME_CTIME                      )

          if ( trim(output_basename) == "" ) then ! [add] H.Yashiro 20120512
             basename='restart_SH_FLUX_SFC'//trim(cdate)
          else
             basename=trim(output_basename)//trim(cdate)
          endif

          call FIO_output( diagvar(:,:,:,I_SH_FLUX_SFC), basename, desc, '', &
                           'SH_FLUX_SFC', 'SH_FLUX_SFC', '', '',             &
                           FIO_REAL8, 'ZSSFC1', 1, 1,                        &
                           1, TIME_CTIME, TIME_CTIME                         )

       elseif( output_io_mode == 'LEGACY' ) then
       ! <- [add] H.Yashiro 20110819
          ! [Add] 2012/06/07 T.Seiki
          if ( trim(output_basename_QV_TB_TEND) == "" ) then 
             basename='restart_QV_TB_TEND'//trim(cdate)
          else
             basename=trim(output_basename_QV_TB_TEND)//trim(cdate)
          endif
          ! 2011/08/16 M.Satoh, change names
          if(output_direct_access) then
             call GTL_output_var2_da( basename, &
                  diagvar(:,:,:,I_QV_TB_TEND), 1, ADM_kall,  &
                  recnum=1, output_size=8 )
          else
             call GTL_output_var2( basename, &
                  diagvar(:,:,:,I_QV_TB_TEND), &
                  diagvar_pl(:,:,:,I_QV_TB_TEND), &
                  1, ADM_kall )
          endif !--- direct/sequencial
          ! [Add] 2012/06/07 T.Seiki
          if ( trim(output_basename_EVAP_SFC) == "" ) then 
             basename='restart_EVAP_SFC'//trim(cdate)
          else
             basename=trim(output_basename_EVAP_SFC)//trim(cdate)
          endif
          if(output_direct_access) then
             call GTL_output_var2_da( basename, &
                  diagvar1(:,:,:,I_EVAP_SFC), 1, ADM_KNONE,  &
                  recnum=1, output_size=8 )
          else
             call GTL_output_var2( basename, &
                  diagvar1(:,:,:,I_EVAP_SFC), &
                  diagvar1_pl(:,:,:,I_EVAP_SFC), &
                  1, ADM_KNONE )
          endif !--- direct/sequencial
          ! [Add] 2012/06/07 T.Seiki
          if ( trim(output_basename_SH_FLUX_SFC) == "" ) then 
             basename='restart_EVAP_SFC'//trim(cdate)
          else
             basename=trim(output_basename_EVAP_SFC)//trim(cdate)
          endif          
          ! 2011/08/16b M.Satoh
          if(output_direct_access) then
             call GTL_output_var2_da( basename, &
                  diagvar1(:,:,:,I_SH_FLUX_SFC), 1, ADM_KNONE,  &
                  recnum=1, output_size=8 )
          else
             call GTL_output_var2( basename, &
                  diagvar1(:,:,:,I_SH_FLUX_SFC), &
                  diagvar1_pl(:,:,:,I_SH_FLUX_SFC), &
                  1, ADM_KNONE )
          endif !--- direct/sequencial
          
       endif  !--- io_mode
    endif


    ! 2011.8.16 M.Satoh
    if ( trim(CP_TYPE) == 'CHIKIRA' ) then
       ksta = diag_ksta_nlayer(I_CBMFX_CHIKIRA)
       kend = diag_kend_nlayer(I_CBMFX_CHIKIRA)
       ! 11/08/17 M.Satoh bug fix: ksta-kend+1 => kend-ksta+1

       ! -> [add] H.Yashiro 20110819
       if ( output_io_mode == 'ADVANCED' ) then
          if ( trim(output_basename) == "" ) then ! [add] H.Yashiro 20120512
             basename='restart_CBMFX_CHIKIRA'//trim(cdate)
          else
             basename=trim(output_basename)//trim(cdate)
          endif

          call FIO_output( diagvarn(:,ksta:kend,:,1), basename, desc, '',             &
                           'CBMFX_CHIKIRA', 'Cloud-base Mass Flux (Chikira)', '', '', &
                           FIO_REAL8, 'GCCKR', 1, kend-ksta+1,                        &
                           1, TIME_CTIME, TIME_CTIME                                  )

       elseif( output_io_mode == 'LEGACY' ) then
       ! <- [add] H.Yashiro 20110819
          ! [Add] 2012/06/07 T.Seiki
          if ( trim(output_basename_CBMFX_CHIKIRA) == "" ) then 
             basename='restart_CBMFX_CHIKIRA'//trim(cdate)
          else
             basename=trim(output_basename_CBMFX_CHIKIRA)//trim(cdate)
          endif          

          if(output_direct_access) then
             call GTL_output_var2_da( basename, &
                                      diagvarn(:,ksta:kend,:,1), 1, kend-ksta+1,  &
                                      recnum=1, output_size=8 )
          else
             call GTL_output_var2( basename, &
                                   diagvarn(:,ksta:kend,:,1), &
                                   diagvarn_pl(:,ksta:kend,:,1), &
                                   1, kend-ksta+1 )
          endif !--- direct/sequencial

       endif  !--- io_mode
    endif

! -> [add&mod] H.Yashiro 20110826
!==> 09/01/28 A.T.Noda
!   if(TB_TYPE=='MY2MOIST')then             ! 07/07/05 A.T.Noda
    if (      trim(TB_TYPE)=='MY2MOIST' &
         .or. trim(TB_TYPE)=='MYNN2'    &
         .or. trim(TB_TYPE)=='MYNN2.5'  &
         .or. trim(TB_TYPE)=='MYNN3'    ) then

       if ( trim(output_basename) == "" ) then ! [add] H.Yashiro 20120512
          basename='restart_diagvar'//trim(cdate)
       else
          basename=trim(output_basename)//trim(cdate)
       endif
       !
       if (      trim(TB_TYPE)=='MY2MOIST' &
            .or. trim(TB_TYPE)=='MYNN2'    &
            .or. trim(TB_TYPE)=='MYNN2.5'  ) then
          !
          if (      trim(TB_TYPE)=='MY2MOIST' &
               .or. trim(TB_TYPE)=='MYNN2'    ) then
             !
             if ( output_io_mode == 'ADVANCED' ) then
                !
                call FIO_output( diagvar(:,:,:,I_QKEd), basename, desc, '',         &
                                 'qked', 'Turbulent Kinetic Energy*2', '', 'J/kg',  &
                                 FIO_REAL8, restart_layername, 1, ADM_kall,         &
                                 1, TIME_CTIME, TIME_CTIME                          )

             elseif( output_io_mode == 'LEGACY' ) then
                !
                if (output_direct_access) then
                   call GTL_output_var2_da( basename,                           &
                                            diagvar(:,:,:,I_QKEd), 1, ADM_kall, &
                                            recnum=1, output_size=8             )
                else
                   call GTL_output_var2( basename,                 &
                                         diagvar(:,:,:,I_QKEd),    &
                                         diagvar_pl(:,:,:,I_QKEd), &
                                         1, ADM_kall               )
                endif !--- direct/sequencial

             endif  !--- io_mode

             write(ADM_LOG_FID,*) 'restart out: qked',           &
                                  maxval(diagvar(:,:,:,I_QKEd)), &
                                  minval(diagvar(:,:,:,I_QKEd))
          endif !--- MY2M, MYNN2

          !
          if ( output_io_mode == 'ADVANCED' ) then
             !
             call FIO_output( diagvar(:,:,:,I_TSQd), basename, desc, '',         &
                              'tsqd', 'Variance of theta_l', '', 'K2',           &
                              FIO_REAL8, restart_layername, 1, ADM_kall,         &
                              1, TIME_CTIME, TIME_CTIME                          )
             call FIO_output( diagvar(:,:,:,I_QSQd), basename, desc, '',         &
                              'qsqd', 'Variance of Total Water(qw)', '', '',     &
                              FIO_REAL8, restart_layername, 1, ADM_kall,         &
                              1, TIME_CTIME, TIME_CTIME                          )
             call FIO_output( diagvar(:,:,:,I_COVd), basename, desc, '',         &
                              'covd', 'Covarriance of qw and theta_l', '', 'K',  &
                              FIO_REAL8, restart_layername, 1, ADM_kall,         &
                              1, TIME_CTIME, TIME_CTIME                          )

          elseif( output_io_mode == 'LEGACY' ) then
             !
             if (output_direct_access) then
                call GTL_output_var2_da( basename,                           &
                                         diagvar(:,:,:,I_TSQd), 1, ADM_kall, &
                                         recnum=2, output_size=8             )
                call GTL_output_var2_da( basename,                           &
                                         diagvar(:,:,:,I_QSQd), 1, ADM_kall, &
                                         recnum=3, output_size=8             )
                call GTL_output_var2_da( basename,                           &
                                         diagvar(:,:,:,I_COVd), 1, ADM_kall, &
                                         recnum=4, output_size=8             )
             else
                call GTL_output_var2( basename,                 &
                                      diagvar(:,:,:,I_TSQd),    &
                                      diagvar_pl(:,:,:,I_TSQd), &
                                      1, ADM_kall               )
                call GTL_output_var2( basename,                 &
                                      diagvar(:,:,:,I_QSQd),    &
                                      diagvar_pl(:,:,:,I_QSQd), &
                                      1, ADM_kall               )
                call GTL_output_var2( basename,                 &
                                      diagvar(:,:,:,I_COVd),    &
                                      diagvar_pl(:,:,:,I_COVd), &
                                      1, ADM_kall               )
             endif !--- direct/sequencial

          endif  !--- io_mode

          write(ADM_LOG_FID,*) 'restart out: tsqd',           &
                               maxval(diagvar(:,:,:,I_TSQd)), &
                               minval(diagvar(:,:,:,I_TSQd))
          write(ADM_LOG_FID,*) 'restart out: qsqd',           &
                               maxval(diagvar(:,:,:,I_QSQd)), &
                               minval(diagvar(:,:,:,I_QSQd))
          write(ADM_LOG_FID,*) 'restart out: covd',           &
                               maxval(diagvar(:,:,:,I_COVd)), &
                               minval(diagvar(:,:,:,I_COVd))
       endif !--- MY2M, MYNN2, MYNN2.5

       if ( RAD_CLOUD_TYPE=='PARTIAL_CLOUD' )then
          !
          if ( output_io_mode == 'ADVANCED' ) then
             !
             call FIO_output( diagvar(:,:,:,I_QCLWB), basename, desc, '',        &
                              'q_clw_before', 'q_clw_before', '', 'kg/kg',       &
                              FIO_REAL8, restart_layername, 1, ADM_kall,         &
                              1, TIME_CTIME, TIME_CTIME                          )
             call FIO_output( diagvar(:,:,:,I_QCLIB), basename, desc, '',        &
                              'q_cli_before', 'q_cli_before', '', 'kg/kg',       &
                              FIO_REAL8, restart_layername, 1, ADM_kall,         &
                              1, TIME_CTIME, TIME_CTIME                          )
             call FIO_output( diagvar(:,:,:,I_CFRACP), basename, desc, '0-1',    &
                              'cfrac_pdf', 'cfrac_pdf', '', '',                  &
                              FIO_REAL8, restart_layername, 1, ADM_kall,         &
                              1, TIME_CTIME, TIME_CTIME                          )

          elseif( output_io_mode == 'LEGACY' ) then
             !
             if (output_direct_access) then
                call GTL_output_var2_da( basename,                             &
                                         diagvar(:,:,:,I_QCLWB), 1, ADM_kall,  &
                                         recnum=5, output_size=8               )
                call GTL_output_var2_da( basename,                             &
                                         diagvar(:,:,:,I_QCLIB), 1, ADM_kall,  &
                                         recnum=6, output_size=8               )
                call GTL_output_var2_da( basename,                             &
                                         diagvar(:,:,:,I_CFRACP), 1, ADM_kall, &
                                         recnum=7, output_size=8               )
             else
                call GTL_output_var2( basename,                   &
                                      diagvar(:,:,:,I_QCLWB),     &
                                      diagvar_pl(:,:,:,I_QCLWB),  &
                                      1, ADM_kall                 )
                call GTL_output_var2( basename,                   &
                                      diagvar(:,:,:,I_QCLIB),     &
                                      diagvar_pl(:,:,:,I_QCLIB),  &
                                      1, ADM_kall                 )
                call GTL_output_var2( basename,                   &
                                      diagvar(:,:,:,I_CFRACP),    &
                                      diagvar_pl(:,:,:,I_CFRACP), &
                                      1, ADM_kall                 )
             endif !--- direct/sequencial

          endif  !--- io_mode

          write(ADM_LOG_FID,*) 'restart out: q_clw_before',    &
                               maxval(diagvar(:,:,:,I_QCLWB)), &
                               minval(diagvar(:,:,:,I_QCLWB))
          write(ADM_LOG_FID,*) 'restart out: q_cli_before',    &
                               maxval(diagvar(:,:,:,I_QCLIB)), &
                               minval(diagvar(:,:,:,I_QCLIB))
          write(ADM_LOG_FID,*) 'restart out: cfrac_pdf',        &
                               maxval(diagvar(:,:,:,I_CFRACP)), &
                               minval(diagvar(:,:,:,I_CFRACP))
       endif !--- PARTIAL_CLOUD

    endif !--- MY2M, MYNN2, MYNN2.5, MYNN3
!<== 09/01/28 A.T.Noda
! <- [add&mod] H.Yashiro 20110826

    ! 10/04/28 M.Satoh
    if( ROUGHNESS_SEA_TYPE == 'YQW' ) then 

       ! -> [add] H.Yashiro 20110819
       if ( output_io_mode == 'ADVANCED' ) then
          if ( trim(output_basename) == "" ) then ! [add] H.Yashiro 20120512
             basename='restart_ROUGHNESS_SEA'//trim(cdate)
          else
             basename=trim(output_basename)//trim(cdate)
          endif
          ! -> [mod] Y.Yamada 20111125
          !          call FIO_output( diagvar(:,:,:,I_ROUGHNESS_SEA), basename, desc, '', &
          !                           'ROUGHNESS_SEA', 'Sea Roughness Length', '', 'm',   &
          !                           FIO_REAL8, 'ZSSFC1', 1, 1,                          &
          !                           1, TIME_CTIME, TIME_CTIME                           )
          call FIO_output( diagvar1(:,:,:,I_ROUGHNESS_SEA), basename, desc, '', &
               'ROUGHNESS_SEA', 'Sea Roughness Length', '', 'm',   &
                           FIO_REAL8, 'ZSSFC1', 1, 1,                          &
                           1, TIME_CTIME, TIME_CTIME                           )
          ! <- [mod] Y.Yamada 20111125
       elseif( output_io_mode == 'LEGACY' ) then
       ! <- [add] H.Yashiro 20110819
          ! [Add] 2012/06/07 T.Seiki
          if ( trim(output_basename_ROUGHNESS_SEA) == "" ) then 
             basename='restart_ROUGHNESS_SEA'//trim(cdate)
          else
             basename=trim(output_basename_ROUGHNESS_SEA)//trim(cdate)
          endif
          
          if (output_direct_access) then
             call GTL_output_var2_da( basename, &
                                      diagvar1(:,:,:,I_ROUGHNESS_SEA), 1, ADM_KNONE,    &
                                      recnum=1, output_size=8 )
          else
             call GTL_output_var2( basename,   &
                                   diagvar1(:,:,:,I_ROUGHNESS_SEA), &
                                   diagvar1_pl(:,:,:,I_ROUGHNESS_SEA),    &
                                   1, ADM_KNONE )
          endif !--- direct/sequencial

       endif  !--- io_mode

       write(ADM_LOG_FID,*) 'restart out: z0_roughness_sea', &
                            maxval(diagvar1(:,:,:,I_ROUGHNESS_SEA)), &
                            minval(diagvar1(:,:,:,I_ROUGHNESS_SEA))
    endif

    ! [mod] 11/08/16 M.Satoh
    call diagvar_timeinfo ( cdate )

      ! del below; 11/08/16 M.Satoh 
      !
      ! 07/12/05 [Add] T.Mitsui
!!$ call diagvar_timeinfo(DIAG_VMAX,diag_cname(1:DIAG_VMAX))
    ! 08/03/10 [Mod] T.Mitsui
    ! 10/04/28 [Mod] M.Satoh

    ! 10/06/09 A.T.Noda
!    call diagvar_timeinfo ( &
!         DIAG_VMAX,         &
!         diag_cname(1:DIAG_VMAX)//trim(cdate), &
!         DIAG_VMAX_1layer,  &
!         diag_cname_1layer(1:DIAG_VMAX_1layer)//trim(cdate) &
!         )
!!$    ! [fix] 10/08/03 T.Mitsui, bug(?) by gfortran
!!$    cwork(:)=""
!!$    cwork_1layer(:)=""
!!$    cwork_nlayer(:)="" ! 11/08/16 M.Satoh
!!$    do iv=1, DIAG_VMAX
!!$       cwork(iv)=diag_cname(iv)//trim(cdate)
!!$    end do
!!$    if( DIAG_VMAX_1layer>0 )then
!!$       do iv=1, DIAG_VMAX_1layer
!!$          cwork_1layer(iv)=diag_cname_1layer(iv)//trim(cdate)
!!$       end do
!!$    end do
!!$    if( DIAG_VMAX_nlayer>0 )then ! 11/08/16 M.Satoh
!!$       do iv=1, DIAG_VMAX_nlayer
!!$          cwork_nlayer(iv)=diag_cname_nlayer(iv)//trim(cdate)
!!$       end do
!!$    end do

!!$       call diagvar_timeinfo ( &
!!$            DIAG_VMAX,         &
!!$            cwork(1:DIAG_VMAX),&
!!$            DIAG_VMAX_1layer,  &
!!$            cwork_1layer(1:DIAG_VMAX_1layer) &
!!$            )
!!$    else
!!$       call diagvar_timeinfo ( &
!!$            DIAG_VMAX,         &
!!$            cwork(1:DIAG_VMAX),&
!!$            DIAG_VMAX_1layer   &
!!$            )
!!$    endif
    !
  end subroutine diagvar_restart_output
  !-----------------------------------------------------------------------------
  ! 11/08/16 [Mod] M.Satoh
  subroutine diagvar_timeinfo &
       ( cdate ) !--- IN
    !
    use mod_misc, only :        &
         MISC_get_available_fid
    use mod_adm, only :         &
         ADM_NSYS,              &
         ADM_prc_me,            &
         ADM_kmin,ADM_kmax,     &
         ADM_prc_run_master
    use mod_time, only :        &
         TIME_DTL, &
         TIME_LSTEP_MAX
    use mod_grd, only :         &
         GRD_gz
    !
    implicit none
    character(len=14), intent(in) :: cdate

    character(len=ADM_NSYS) :: cname(DIAG_VMAX_DEF)
    character(len=ADM_NSYS) :: cname_1layer(DIAG_VMAX_DEF)
    character(len=ADM_NSYS) :: cname_nlayer(DIAG_VMAX_DEF)

    integer :: num = 1
    integer :: fid
    integer :: k
    integer :: n
    integer :: iv
    character(len=3) :: wc       ! [Add] 2012/06/07 T.Seiki

    cname(:) = ""
    cname_1layer(:) = ""
    cname_nlayer(:) = ""

    do iv = 1, DIAG_VMAX
       cname(iv) = diag_cname(iv)//trim(cdate)
    end do
    if( DIAG_VMAX_1layer > 0 )then
       do iv = 1, DIAG_VMAX_1layer
          cname_1layer(iv) = diag_cname_1layer(iv)//trim(cdate)
       end do
    end if
    if( DIAG_VMAX_nlayer > 0 )then ! 11/08/16 M.Satoh
       do iv = 1, DIAG_VMAX_nlayer
          cname_nlayer(iv) = diag_cname_nlayer(iv)//trim(cdate)
       end do
    end if

    if(ADM_prc_me == ADM_prc_run_master) then
       !
       fid = MISC_get_available_fid()
          open(fid,file='diagvar.info',form='formatted',status='replace')
       do n = 1, DIAG_VMAX
          write(fid,'(I4,F16.2)') 1, TIME_DTL*TIME_LSTEP_MAX
          write(fid,'(I4)') ADM_kall
          do k=1, ADM_kall
             write(fid,'(F16.4)') GRD_gz(k)
          end do
          write(fid,'(I4)') num
          write(fid,'(a32)') trim(cname(n))
       end do
       if ( DIAG_VMAX_1layer > 0 ) then
          do n = 1, DIAG_VMAX_1layer
             write(fid,'(I4,F16.2)') 1, TIME_DTL*TIME_LSTEP_MAX
             write(fid,'(I4)') ADM_KNONE
             write(fid,'(F16.4)') 0.0
             write(fid,'(I4)') num
             write(fid,'(a32)') trim(cname_1layer(n))
          end do
       end if
       if ( DIAG_VMAX_nlayer > 0 ) then
          do n = 1, DIAG_VMAX_nlayer
             write(fid,'(I4,F16.2)') 1, TIME_DTL*TIME_LSTEP_MAX
             write(fid,'(I4)') diag_knum_nlayer(n)
             do k=1, diag_knum_nlayer(n)
                write(fid,'(F16.4)') float(k)
             end do
             write(fid,'(I4)') num
             write(fid,'(a32)') trim(cname_nlayer(n))
          end do
       end if
       close(fid)
       !
    end if
    !
    return
    !
  end subroutine diagvar_timeinfo

end module mod_diagvar
