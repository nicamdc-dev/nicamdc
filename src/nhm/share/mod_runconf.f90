!-------------------------------------------------------------------------------
!>
!! run configuration module
!!
!! @par Description
!!         admin modlue for 3D-model
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)   Imported from igdc-4.34
!! @li      2005-10-28 (M.Satoh)    add ISCCP parameter
!! @li      2005-12-27 (M.Satoh)    introduce RAD_DIV_NUM
!! @li      2006-04-19 (H.Tomita)   Add 'DRY' option.
!! @li      2006-05-06 (H.Tomita)   abolish TURB_DIV_NUM.
!! @li      2006-08-11 (H.Tomita)   Add TRC_ADV_TYPE, TRC_NEG_FIX
!! @li      2006-09-28 (S.Iga)      Add OUT_FILE_TYPE
!! @li      2007-01-26 (H.Tomita)   Add the 'SIMPLE2' EIN_TYPE.
!!                                  Control the EIN_TYPE in thie module. as CVW(I_Q?) and LHV,LHF, and LHS.
!! @li      2007-03-23 (Y.Niwa)     add FLAG_NUDGING
!! @li      2007-05-08 (H.Tomita)   1. Move physics type configuration from mod_physicsinit.f90 to here.
!!                                  2. Add TRC_ADD_MAX for future implementation of turbulence scheme.
!! @li      2007-06-26 (Y.Niwa)     move LAND_TYPE from mod_land_driver to here
!!                                  move OCEAN_TYPE from mod_ocean_driver to here
!! @li      2007-07-17 (A.T.Noda)   Add RAD_CLOUD_TYPE for use of the partial cloud in rd_driver
!! @li      2007-07-23 (K.Suzuki)   SPRINTARS aerosol model
!!                                  1. Add number of aerosol tracers
!!                                  2. Add KAPCL: number of aerosol for radiation
!!                                  3. Add AE_TYPE for aerosol type configuration
!! @li      2007-11-07 (T.Mitsui)   add "OPT_OUTPUT_ALL" to omit output_all
!! @li      2008-01-24 (Y.Niwa)     add MIURA2004OLD
!!                                     TRC_ADV_TYPE='DEFAULT', TRC_NEG_FIX='ON'
!!                                  => TRC_ADV_TYPE='MIURA2004', TRC_NEG_FIX='OFF'
!!                                  add TRC_SAVE_MEMORY
!! @li      2008-01-30 (Y.Niwa)     bug fixes
!! @li      2008-03-10 (T.Mitsui)   add intermediate output of restart file
!! @li      2008-04-12 (T.Mitsui)   add 2-moment hydrometeors, incloud aerosol and their labeling
!! @li      2008-04-23 (H.Tomita)   Add MP_DIVNUM for microphysics.
!! @li      2008-08-11 (H.Tomita)   Add SFC_DIV_NUM and TB_DIV_NUM.
!! @li      2008-10-05 (T.Mitsui)   add option : ALL_PHYSTEP_POST
!! @li      2009-01-28 (A.T.Noda)   Implement MYNN
!! @li      2009-04-14 (T.Mitsui)   add opt_carb, opt_dust, opt_salt, opt_sulf trivial changing incloud aerosols
!!                                  add opt_aerosol_forcing for ARF
!! @li      2009-07-28 (H.Tomita)   add PRCIP_TRN_ECORRECT
!!                                  for energy adjustment in the rain-sedimentation process
!! @li      2010-03-08 (C.Kodama)   Add overwrite_restart option
!! @li      2010-04-26 (M.Satoh)    add ROUGHNESS_SEA_TYPE
!! @li      2010-06-19 (A.T.Noda)   Allow to use a convection parameterization
!!                                  with an advanced microphysics schemes, such as G98, NSW?,
!! @li      2010-11-11 (A.T.Noda)   1. add CORIOLIS, RAD_FIX_LAT/LON for Giga-LES
!! @li      2011-06-30 (T.Seiki)    fill undefined indices
!! @li      2011-07-22 (T.Ohno)     add CORIOLIS_PARAM
!! @li      2011-09-03 (H.Yashiro)  add TRC_name for New I/O
!! @li      2012-02-01 (T.Seiki)    add incloud aerosol indices+initialization
!! @li      2012-07-23 (H.Yashiro)  [add] River         by K.Yoshimura
!! @li      2012-07-23 (H.Yashiro)  [add] Water Isotope by K.Yoshimura
!! @li      2012-11-05 (H.Yashiro)  NICAM milestone project (Phase I:cleanup of shared module)
!!
!<
module mod_runconf
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_NSYS
  use mod_cnst, only: & ! [add] 10/11/11 A.Noda
     CNST_UNDEF
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: runconf_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=16),       public, allocatable, save  :: TRC_name(:) ! short name  of tracer [add] H.Yashiro 20110819
  character(len=ADM_NSYS), public, allocatable, save  :: WLABEL  (:) ! description of tracer

  integer, public, save :: TRC_vmax =  0 ! total number of tracers
  integer, public, save :: NQW_MAX  =  0 ! subtotal number of water mass tracers
  integer, public, save :: NQW_STR  = -1 ! start index of water mass tracers
  integer, public, save :: NQW_END  = -1 ! end   index of water mass tracers

  integer, public, save :: NNW_MAX  =  0 ! subtotal number of water number tracers
  integer, public, save :: NNW_STR  = -1 ! start index of water number tracers
  integer, public, save :: NNW_END  = -1 ! end   index of water number tracers

  integer, public, save :: NTB_MAX  =  0 ! subtotal number of turbulent tracers

  integer, public, save :: NCHEM_MAX  =  0 ! subtotal number of chemical (or general purpose) tracers
  integer, public, save :: NCHEM_STR  = -1 ! start index of chemical (or general purpose) tracers
  integer, public, save :: NCHEM_END  = -1 ! end   index of chemical (or general purpose) tracers

  integer, public, save :: NQA_MAX  =  0 ! subtotal number of aerosol mass tracers
  integer, public, save :: NQA_STR  = -1 ! start index of aerosol mass tracers
  integer, public, save :: NQA_END  = -1 ! end   index of aerosol mass tracers

  integer, public, save :: NDU_MAX  =  0 ! soil dust
  integer, public, save :: NQDU_STR = -1
  integer, public, save :: NQDU_END = -1
  integer, public, save :: NCB_MAX  =  0 ! carbonaceous
  integer, public, save :: NQCB_STR = -1
  integer, public, save :: NQCB_END = -1
  integer, public, save :: NSU_MAX  =  0 ! sulfate
  integer, public, save :: NQSU_STR = -1
  integer, public, save :: NQSU_END = -1
  integer, public, save :: NSA_MAX  =  0 ! sea salt
  integer, public, save :: NQSA_STR = -1
  integer, public, save :: NQSA_END = -1

  ! [add] 12/02/01, T.Seiki incloud aerosols are predicted as a passive tracer
  integer, public, save :: NQAIN_MAX  =  0 ! subtotal number of incloud aerosol mass tracers
  integer, public, save :: NQAIN_STR  = -1 ! start index of incloud aerosol mass tracers
  integer, public, save :: NQAIN_END  = -1 ! end   index of incloud aerosol mass tracers

  integer, public, save :: NQDUIN_MAX =  0 ! incloud soil dust
  integer, public, save :: NQDUIN_STR = -1
  integer, public, save :: NQDUIN_END = -1
  integer, public, save :: NQCBIN_MAX =  0 ! incloud carbonaceous
  integer, public, save :: NQCBIN_STR = -1
  integer, public, save :: NQCBIN_END = -1
  integer, public, save :: NQSUIN_MAX =  0 ! incloud sulfate
  integer, public, save :: NQSUIN_STR = -1
  integer, public, save :: NQSUIN_END = -1
  integer, public, save :: NQSAIN_MAX =  0 ! incloud sea salt
  integer, public, save :: NQSAIN_STR = -1
  integer, public, save :: NQSAIN_END = -1
  ! [add] 12/02/01, T.Seiki <=

  !--- preserved index name for tracers
  integer, public, save :: I_QV     = -1  ! Water vapor
  integer, public, save :: I_QC     = -1  ! Cloud water
  integer, public, save :: I_QR     = -1  ! Rain
  integer, public, save :: I_QI     = -1  ! Ice
  integer, public, save :: I_QS     = -1  ! Snow
  integer, public, save :: I_QG     = -1  ! Graupel
  integer, public, save :: I_QH     = -1  ! Hail (not used)

  integer, public, save :: I_NC     = -1  ! Cloud water (number)
  integer, public, save :: I_NR     = -1  ! Rain        (number)
  integer, public, save :: I_NI     = -1  ! Ice         (number)
  integer, public, save :: I_NS     = -1  ! Snow        (number)
  integer, public, save :: I_NG     = -1  ! Graupel     (number)

  integer, public, save :: I_TKE    = -1 ! turbulence kinetic energy

  integer, public, save :: I_QKEp   = -1
  integer, public, save :: I_TSQp   = -1
  integer, public, save :: I_QSQp   = -1
  integer, public, save :: I_COVp   = -1

  ! 09/04/14 T.Mitsui [Add] =>
  integer, public, save :: NQSU_SO4 = -1 ! SO4(2-)
  integer, public, save :: NQSU_SO2 = -1 ! SO2(gas and liquid)
  integer, public, save :: NQSU_DMS = -1 ! Dimethyl Sulphide(water-insoluble liquid)
  ! 09/04/14 T.Mitsui [Add] <=
  integer, public, save :: NQSO4IN  = -1  ! incloud SO4(2-) [add] 12/02/12, T.Seiki

  !--- Nonhydrostatic/hydrostatic flag
  integer,                 public, save :: NON_HYDRO_ALPHA = 1
  !                                        1 : Non-hydrosataic
  !                                        0 : Hydrosataic

  !--- Location of numerical filtering
  character(len=ADM_NSYS), public, save :: NDIFF_LOCATION = 'IN_LARGE_STEP'
  !                                        'IN_LARGE_STEP'  : insert in  the large step

  !--- Number of division of numerical filtering at the out large step
  ! Note: this number is used only at a large diffusion
  integer, public, save :: NDIFF_DIVISION_NUM = 1 ! division number of num. diff.

  !--- Estimation scheme for the internal energy
  character(len=ADM_NSYS), public, save :: EIN_TYPE = 'EXACT'
  !                                     'SIMPLE': standard approximation CVD * T
  !                                     'EXACT': exact formulation
  !                       -> if warm rain
  !                           qd*CVD*T + qv*CVV*T + (qc+qr)*CPL*T
  !                       -> if cold rain
  !                           qd*CVD*T + qv*CVV*T + (qc+qr)*CPL*T
  !                                               + (qi+qs)*CPI*T

  !--- output file type
  character(len=ADM_NSYS), public, save :: OUT_FILE_TYPE = 'DEFAULT'

  !--- run type
  character(len=ADM_NSYS), public, save :: RUN_TYPE = 'DEFAULT'

  !--- rain type
  character(len=ADM_NSYS), public, save :: RAIN_TYPE = 'WARM'
  !                                                    'COLD'
  !                                                    'CLOUD_PARAM'
  !                                                    'DRY'
 
  !--- radiation type
  character(len=ADM_NSYS), public, save :: RAD_TYPE = 'ISCCP'

  character(len=ADM_NSYS), public, save :: RAD_CLOUD_TYPE = 'DEFAULT'
  !                                                         'PARTIAL_CLOUD': Use partial cloud computed in turbulent module

  !--- Tracer advection scheme
  character(len=ADM_NSYS), public, save :: TRC_ADV_TYPE = 'MIURA2004'    ! Y.Niwa 080124
  !                                                       'DEFAULT'

  !--- negative fixer of water contents
  character(len=ADM_NSYS), public, save :: TRC_NEG_FIX = 'OFF'  ! Y.Niwa 080124
  !                                                      'ON'

  !--- option for saving memory in tracer advection  ! Y.Niwa 080124
  logical, public, save :: TRC_SAVE_MEMORY = .true.

  !--- radiation dividion number
  integer, public, save :: RAD_DIV_NUM = 1
  integer, public, save :: TB_DIV_NUM  = 1
  integer, public, save :: SFC_DIV_NUM = 1
  !
  !
  logical, public, save :: FLAG_NUDGING = .false. ! 07/03/23 Y.Niwa [add]

  logical, public, save :: OPT_OUTPUT_ALL = .false. ! 07/11/07 T.Mitsui [add]
  !
  ! 08/03/10 T.Mitsui [Add]
  !--- number of intermediate restart file, and their date
  !    date should be described as "yyyymmddhhmmss", (ex. "20080307214700"
  integer,           public, save :: num_restart = -1   ! max num.= 64
  character(len=14), public, save :: cdate_restart(32)  !
  character(len=14), public, save :: overwrite_restart  ! [Add] 10/03/08 C.Kodama

  !--- Latent heat
  real(8), public, save :: LHV
  real(8), public, save :: LHF
  real(8), public, save :: LHS
  !
  integer, public, save :: MP_DIV_NUM = 1

  !--- Max number of hydrometeor-species
  integer, parameter,   public        :: HYDRO_MAX=7 ! 08/04/12 [Add] T.Mitsui

  !--- specific heat of water on const pressure
  real(8), allocatable, public, save :: CVW(:)
  real(8), allocatable, public, save :: CPW(:)


  character(len=ADM_NSYS), public, save :: TSAT_TYPE = 'LIST'
  !                               = 'LIST': list vector: default
  !                               = 'NOLIST': no list vecto

  character(len=ADM_NSYS), public, save :: TSAT_TYPE_ENT = 'LIST'
  !                               = 'LIST': list vector: default
  !                               = 'NOLIST': no list vecto

  !==== radiation parameter
  logical, public, save :: RAD_GLOBAL_TROPICS = .false.
  logical, public, save :: CLEAR_SKY_RAD = .false.
  ! [Add] 09/04/14 T.Mitsui, moved from rd_mstrnx_ar5
  ! option for estimating aerosol direct radiative forcing
  logical, public, save :: opt_aerosol_forcing=.false.
  ! [Add] 10/11/14 A.Noda for Giga-LES
  logical, public, save :: CORIOLIS    = .true.
  real(8), public, save :: RAD_FIX_LON = CNST_UNDEF
  real(8), public, save :: RAD_FIX_LAT = CNST_UNDEF
  !
  !---  No.of radiation calculation
  !  integer, parameter, public :: NCRF =  1
  integer, parameter, public :: NCRF =  2
  !
  !---  No. of band for rad.
  integer, parameter, public :: NRBND = 3
  !
  integer, parameter, public :: NRBND_VIS = 1
  integer, parameter, public :: NRBND_NIR = 2
  integer, parameter, public :: NRBND_IR = 3
  !
  !---  direct/diffuse
  integer, parameter, public :: NRDIR = 2
  integer, parameter, public :: NRDIR_DIRECT = 1
  integer, parameter, public :: NRDIR_DIFFUSE = 2
  !
  !--- ISCCP parameter 05/10/28 M.Satoh
  integer, parameter, public :: NTAU_ISCCP = 7
  integer, parameter, public :: NPRES_ISCCP = 7
  !
  !
  !==== roughness  parameter
  integer, public, parameter :: NTYPE_Z0 = 3
  integer, public, parameter :: N_Z0M = 1
  integer, public, parameter :: N_Z0H = 2
  integer, public, parameter :: N_Z0E = 3
  !
  !--- Physics type
  character(len=ADM_NSYS), public, save :: CP_TYPE = 'NONE'   ! Y.Niwa add save 080130
  !                                      = 'KUO'

  character(len=ADM_NSYS), public, save :: MP_TYPE = 'O2001'   ! Y.Niwa add save 080130
  !                                      = 'NONE'
  !                                      = 'KESSLER' Kessler(1969)
  !                                      = 'KW78' Klemp & Wilhemson(1978,JAS)
  !                                         tentatively same as O2000
  !                                      = 'O2001' Ooyama(2001,JAS)
  !                                         equivalent to KW78
  !                                         but more feasible than KW78
  !                                      = 'G98' Grabowski(1998,JAS,3283-3298)
  !                                      = 'G98W5' modified G98
  !                                         use of five categories of water
  !                                      = 'G99' Grabowski(1999,Atm.Res.,17-41)

  character(len=ADM_NSYS), public, save :: TB_TYPE = 'NONE'   ! Y.Niwa add save 080130
  !                                = 'CONST'
  !                                = 'SMAGORINSKY'
  !                                = 'DEARDORFF'
  !                                = 'MY2'
  !                                = 'MY2.5'
  !                                = 'MY2MOIST' ! 090128 A.T.Noda
  !                                = 'MYNN2.5'
  !                                = 'MYNN3'

  character(len=ADM_NSYS), public, save :: LAND_TYPE  = 'BUCKET'     ! Y.Niwa add save 080130

  character(len=ADM_NSYS), public, save :: RIV_TYPE   = 'NONE'       ! [add] K.Yoshimura 20110414

  character(len=ADM_NSYS), public, save :: OCEAN_TYPE = 'SST'        ! Y.Niwa add save 080130

  character(len=ADM_NSYS), public, save :: AE_TYPE = 'NONE'   ! Y.Niwa add save 080130
  !                                                = 'SPRINTARS'

  character(len=ADM_NSYS), public, save :: CHEM_TYPE = 'NONE'   ! Y.Niwa add save 080130
  !                                                  = 'PASSIVE'
  !                                                  = 'CHASER'

  character(len=ADM_NSYS), public, save :: RD_TYPE = 'NONE'     ! Y.Niwa add save 080130
  !                                 = 'CONST'     : specified heating
  !                                 = 'FIX' = 'CONST'
  !                                 = 'MSTRNX': Nakajima rad from AGCM5.6

  character(len=ADM_NSYS), public, save :: SV_TYPE_T = 'ISOTHERM'
  !                             = 'ISOTHERM'
  !                             = 'HELD-SUAREZ'  : Held & Suarez

  character(len=ADM_NSYS), public, save :: SV_TYPE_QV = 'DRY'  ! Y.Niwa add save 080130
  !                             = 'SATURATION'
  !                             = 'DRY'

  character(len=ADM_NSYS), public, save :: SF_TYPE = 'NO-FLUX'
  !                                = 'CONST'
  !                                = 'NO-FLUX'
  !                                = 'BULK-CONST'
  !                                = 'NEUTRAL'
  !                                = 'LOUIS'

  character(len=ADM_NSYS), public, save :: FR_TYPE = 'NONE'  ! Y.Niwa add save 080130
  !                             = 'NONE'     : nothing
  !                             = 'HELD-SUAREZ'  : Held & Suarez
  !                             = 'RIGID-ROTATION' : nudging to Rigid rotation


  !--- roughness type
  character(len=ADM_NSYS), public, save :: ROUGHNESS_SEA_TYPE &
                                   = 'MILLER' ! : default as CCSR/NIES 5.6
  !                                = 'YQW'    ! : Yuqing Wang [add] 2010.4.26 M.Satoh

  character(len=ADM_NSYS), public, save :: ND_TYPE = 'NONE'   ! Y.Niwa add save 080130

  character(len=ADM_NSYS), public, save :: AF_TYPE = 'NONE'   ! Y.Niwa add save 080130

  character(len=ADM_NSYS), public, save :: GWD_TYPE = 'NONE'    ! Y.Niwa add save 080130


  !--- number of aerosol species for radiation : 07/07/23 K.Suzuki for SPRINTARS
  integer, parameter, public :: KAPCL = 7
  !<= K.Suzuki add for SPRINTARS 07/07/23
  !
  ! 08/04/12 T.Mitsui add ==>
  ! for 2-moment bulk cloud microphysics
  logical, public, save :: opt_2moment_water = .false.
  logical, public, save, allocatable :: flag_diagnose_number(:)
  ! for accurate incloud scavenging
  logical, public, save :: opt_incloud_aerosol = .false.
  ! option for idealized or regional experiments
  logical, public, save :: opt_carb_on = .true.
  logical, public, save :: opt_dust_on = .true.
  logical, public, save :: opt_salt_on = .true.
  logical, public, save :: opt_sulf_on = .true.
  !
  character(len=ADM_NSYS), public, save ::  PRCIP_TRN_ECORRECT = 'KIN2KIN'
  !                                                          'KIN2EIN'
  real(8), public, save :: CORIOLIS_PARAM = CNST_UNDEF ! [add] T.Ohno 110722

  !=> [add] K.Yoshimura 20110414
  character(len=ADM_NSYS), public, save ::  ISOTOPE = 'OFF'
  character(len=ADM_NSYS), public, save ::  NOFRAC  = 'OFF'
  integer, public, save :: ISO_MAX  = 0
  integer, public, save :: ISO_STR  = -1
  integer, public, save :: ISO_STR2 = -1
  integer, public, save :: ISO_END  = -1

  integer, public, save :: ISO1_QV = -1
  integer, public, save :: ISO1_QC = -1
  integer, public, save :: ISO1_QR = -1
  integer, public, save :: ISO1_QI = -1
  integer, public, save :: ISO1_QS = -1
  integer, public, save :: ISO1_QG = -1
  integer, public, save :: ISO2_QV = -1
  integer, public, save :: ISO2_QC = -1
  integer, public, save :: ISO2_QR = -1
  integer, public, save :: ISO2_QI = -1
  integer, public, save :: ISO2_QS = -1
  integer, public, save :: ISO2_QG = -1
  !<= [add] K.Yoshimura 20110414

  !--- limiter switch
  logical, public, save :: THUBURN_LIM = .true.  ![add] 20130613 R.Yoshida

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
  !>
  !> Description of the subroutine RUNCONF_setup
  !>
  subroutine RUNCONF_setup
    use mod_adm, only :  &
         ADM_LOG_FID,    &
         ADM_CTL_FID,    &
         ADM_proc_stop
    use mod_cnst, only: &
         CNST_CV,    &
         CNST_CP,    &
         CNST_CVV,   &
         CNST_CPV,   &
         CNST_CL,    &
         CNST_CI,    &
         CNST_LH00,  &
         CNST_LH0,   &
         CNST_LHS00, &
         CNST_LHS0,  &
         CNST_LHF00, &
         CNST_LHF0
    use mod_chemvar, only: &
       CHEM_TRC_vmax, &
       CHEM_TRC_name, &
       CHEM_TRC_desc
    implicit none

    namelist /RUNCONFPARAM/ &
         RUN_TYPE,          &
         NON_HYDRO_ALPHA,   & !--- hydro/non-hydro flag
         RD_TYPE,           &
         ND_TYPE,           &
         AF_TYPE,           &
         FR_TYPE,           &
         TB_TYPE,           &
         SV_TYPE_T,         &
         SV_TYPE_QV,        &
         SF_TYPE,           &
         MP_TYPE,           &
         CP_TYPE,           &
         GWD_TYPE,          &
         LAND_TYPE,         &  ! Y.Niwa add 070627
         OCEAN_TYPE,        &  ! Y.Niwa add 070627
         ROUGHNESS_SEA_TYPE,&  ! M.Satoh add 2010.4.26
         AE_TYPE,           &  ! K.Suzuki add 07/07/23 [add]
         CHEM_TYPE,         &
         RAIN_TYPE,         &
         RAD_TYPE ,         & !--- radiation type with/without ISCCP 05/10/28
         TRC_ADV_TYPE,      &
         NDIFF_LOCATION,    & !--- location of numerical diff.
         NDIFF_DIVISION_NUM,& !--- division number of numerical diffusion
         RAD_CLOUD_TYPE ,   & !--- cloud type for radiative computation 07/07/05 A.T.Noda
         EIN_TYPE,          & !--- scheme for temperature estimate
         RAD_GLOBAL_TROPICS,&
         CORIOLIS,          & ! [add] 10/11/14 A.Noda
         CORIOLIS_PARAM,    & ! [add] 11/07/22 T.Ohno
         RAD_FIX_LON,       & ! [add] 10/11/14 A.Noda
         RAD_FIX_LAT,       & ! [add] 10/11/14 A.Noda
         CLEAR_SKY_RAD,     &
         RAD_DIV_NUM,       &
         TRC_NEG_FIX,       &
         OUT_FILE_TYPE,     & ! add Iga(060927) determine the way of output (e.g. 'cfmip')
         OPT_OUTPUT_ALL,    & ! add option to omit output_all
         num_restart,       & ! add 07/03/10 T.Mitsui
         cdate_restart,     & ! add 07/03/10 T.Mitsui
         overwrite_restart, & ! [Add] 10/03/08 C.Kodama
         FLAG_NUDGING,      & !--- nudging flag 07/03/23 Y.Niwa [add]
         opt_2moment_water,  & ! 08/04/12 T.Mitsui add
         opt_incloud_aerosol,& ! 08/04/12 T.Mitsui add
         opt_carb_on, opt_dust_on, & ! 09/04/14 T.Mitsui
         opt_salt_on, opt_sulf_on, & ! 09/04/14 T.Mitsui
         opt_aerosol_forcing, & ! [Add] 09/04/14 T.Mitsui
         MP_DIV_NUM,          &  ! 08/04/23 [add] H.Tomita
         TB_DIV_NUM,          &
         SFC_DIV_NUM,         &
         PRCIP_TRN_ECORRECT,  &
         THUBURN_LIM             ! R.Yoshida 13/06/13 [add]

    character(len=2) :: is ! [add] H.Yashiro 20110819

    integer :: ierr
    integer :: nq, i
    !---------------------------------------------------------------------------

    cdate_restart(:)="YYYYMMDDHHMMSS" ! add 07/03/10 T.Mitsui
    overwrite_restart=""              ! [Add] 10/03/08 C.Kodama

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[runconf]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=RUNCONFPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** RUNCONFPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist RUNCONFPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist RUNCONFPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,RUNCONFPARAM)

    !--- counting tracer
    TRC_vmax = 0

    !--- Mass tracer for Water
    if    ( RAIN_TYPE == 'DRY' ) then
       NQW_MAX = 1
       I_QV    = TRC_vmax + 1
    elseif( RAIN_TYPE == 'CLOUD_PARAM' ) then
       NQW_MAX = 2
       I_QV    = TRC_vmax + 1
       I_QC    = TRC_vmax + 2
    elseif( RAIN_TYPE == 'WARM' ) then
       NQW_MAX = 3
       I_QV    = TRC_vmax + 1
       I_QC    = TRC_vmax + 2
       I_QR    = TRC_vmax + 3
    elseif( RAIN_TYPE == 'COLD' ) then
       NQW_MAX = 6
       I_QV    = TRC_vmax + 1
       I_QC    = TRC_vmax + 2
       I_QR    = TRC_vmax + 3
       I_QI    = TRC_vmax + 4
       I_QS    = TRC_vmax + 5
       I_QG    = TRC_vmax + 6
    else
       write(*,          *) 'xxx You must set RAIN_TYPE to DRY,CLOUD_PARAM,WARM or COLD. STOP.'
       write(ADM_LOG_FID,*) 'xxx You must set RAIN_TYPE to DRY,CLOUD_PARAM,WARM or COLD. STOP.'
       call ADM_proc_stop
    endif
    !--- index range
    NQW_STR = TRC_vmax + 1
    NQW_END = TRC_vmax + NQW_MAX

    !--- update total number
    TRC_vmax = TRC_vmax + NQW_MAX

    !--- Number tracer for Water
    ! 08/04/12 T.Mitsui =>
    if ( opt_2moment_water ) then
       if    ( RAIN_TYPE == 'DRY' ) then
          NNW_MAX = 0
       elseif( RAIN_TYPE == 'CLOUD_PARAM' ) then
          NNW_MAX = 1
          I_NC    = TRC_vmax + 1
       elseif( RAIN_TYPE == 'WARM' ) then
          NNW_MAX = 2
          I_NC    = TRC_vmax + 1
          I_NR    = TRC_vmax + 2
       elseif( RAIN_TYPE == 'COLD' ) then
          NNW_MAX = 5
          I_NC    = TRC_vmax + 1
          I_NR    = TRC_vmax + 2
          I_NI    = TRC_vmax + 3
          I_NS    = TRC_vmax + 4
          I_NG    = TRC_vmax + 5
       endif
       !--- index range
       NNW_STR = TRC_vmax + min(1,NNW_MAX)
       NNW_END = TRC_vmax + NNW_MAX

       !--- update total number
       TRC_vmax = TRC_vmax + NNW_MAX
    endif
    ! 08/04/12 T.Mitsui <=

    !--- Tracer for turbulence
    if    ( trim(TB_TYPE) == 'MY2.5' ) then
       NTB_MAX = 1
       I_TKE   = TRC_vmax + 1
    elseif( trim(TB_TYPE) == 'MYNN2.5' ) then !=> 09/01/28 A.T.Noda
       NTB_MAX = 1
       I_QKEp  = TRC_vmax + 1
    elseif( trim(TB_TYPE) == 'MYNN3' ) then
       NTB_MAX = 4
       I_QKEp  = TRC_vmax + 1
       I_TSQp  = TRC_vmax + 2
       I_QSQp  = TRC_vmax + 3
       I_COVp  = TRC_vmax + 4
    endif                                     !<= 09/01/28 A.T.Noda
    !--- update total number
    TRC_vmax = TRC_vmax + NTB_MAX

    !--- Tracer for chemistry
    if (      trim(CHEM_TYPE) == 'PASSIVE' &
         .OR. trim(CHEM_TYPE) == 'CHASER'  )then
       NCHEM_MAX = CHEM_TRC_vmax
       NCHEM_STR = TRC_vmax + min(1,NCHEM_MAX)
       NCHEM_END = TRC_vmax + NCHEM_MAX
    endif
    TRC_vmax = TRC_vmax + NCHEM_MAX

    !--- Aerosol tracer
    if (      trim(AE_TYPE) == 'SPRINTARS'     &
         .OR. trim(AE_TYPE) == 'SPRINTARS_CRM' ) then

       NQA_MAX = 0

       if ( opt_dust_on ) then ! soil dust
          NDU_MAX  = 10
          NQDU_STR = TRC_vmax + NQA_MAX + 1
          NQDU_END = TRC_vmax + NQA_MAX + NDU_MAX
       endif
       NQA_MAX = NQA_MAX + NDU_MAX

       if ( opt_carb_on ) then ! carbonaceous
          NCB_MAX  = 7
          NQCB_STR = TRC_vmax + NQA_MAX + 1
          NQCB_END = TRC_vmax + NQA_MAX + NCB_MAX
       endif
       NQA_MAX = NQA_MAX + NCB_MAX

       if ( opt_sulf_on ) then ! sulfate
          NSU_MAX  = 3
          NQSU_STR = TRC_vmax + NQA_MAX + 1
          NQSU_END = TRC_vmax + NQA_MAX + NSU_MAX

          ! 09/04/14 T.Mitsui [Add]
          NQSU_SO4 = TRC_vmax + NQA_MAX + 1
          NQSU_SO2 = TRC_vmax + NQA_MAX + 2
          NQSU_DMS = TRC_vmax + NQA_MAX + 3
       endif
       NQA_MAX = NQA_MAX + NSU_MAX

       if ( opt_salt_on ) then ! sea salt
          NSA_MAX  = 4
          NQSA_STR = TRC_vmax + NQA_MAX + 1
          NQSA_END = TRC_vmax + NQA_MAX + NSA_MAX
       endif
       NQA_MAX = NQA_MAX + NSA_MAX

       !--- index range
       NQA_STR = TRC_vmax + min(1,NQA_MAX)
       NQA_END = TRC_vmax + NQA_MAX

       !--- update total number
       TRC_vmax = TRC_vmax + NQA_MAX


       !--- Incloud aerosol tracer
       if( opt_2moment_water ) opt_incloud_aerosol = .true.

       if ( opt_incloud_aerosol ) then
          ! 09/04/14 [Mod] T.Mitsui, option works

          NQAIN_MAX = 0

          if ( opt_dust_on ) then ! soil dust
             NQDUIN_MAX = NDU_MAX
             NQDUIN_STR = TRC_vmax + NQAIN_MAX + 1
             NQDUIN_END = TRC_vmax + NQAIN_MAX + NQDUIN_MAX
          endif
          NQAIN_MAX = NQAIN_MAX + NQDUIN_MAX

          if ( opt_carb_on ) then ! carbonaceous
             NQCBIN_MAX = NCB_MAX
             NQCBIN_STR = TRC_vmax + NQAIN_MAX + 1
             NQCBIN_END = TRC_vmax + NQAIN_MAX + NQCBIN_MAX
          endif
          NQAIN_MAX = NQAIN_MAX + NQCBIN_MAX

          if ( opt_sulf_on ) then ! sulfate
             NQSUIN_MAX = 1
             NQSUIN_STR = TRC_vmax + NQAIN_MAX + 1
             NQSUIN_END = TRC_vmax + NQAIN_MAX + NQSUIN_MAX

             NQSO4IN = TRC_vmax + NQAIN_MAX + 1
          endif
          NQAIN_MAX = NQAIN_MAX + NQSUIN_MAX

          if ( opt_salt_on ) then ! sea salt
             NQSAIN_MAX = NSA_MAX
             NQSAIN_STR = TRC_vmax + NQAIN_MAX + 1
             NQSAIN_END = TRC_vmax + NQAIN_MAX + NQSAIN_MAX
          endif
          NQAIN_MAX = NQAIN_MAX + NQSAIN_MAX

          !--- index range
          NQAIN_STR = TRC_vmax + min(1,NQAIN_MAX)
          NQAIN_END = TRC_vmax + NQAIN_MAX

          !--- update total number
          TRC_vmax = TRC_vmax + NQAIN_MAX
       endif
    endif

    !=> [add] K.Yoshimura 20110414
    if(trim(ISOTOPE)=='ON') then
       if(RAIN_TYPE=='WARM') then
          ISO_MAX = 6
          ISO_STR = TRC_VMAX + 1
          ISO_STR2 = TRC_VMAX + 4
          ISO_END = TRC_VMAX + ISO_MAX
          ISO1_QV = ISO_STR
          ISO1_QC = ISO_STR + 1
          ISO1_QR = ISO_STR + 2
          ISO2_QV = ISO_STR2 
          ISO2_QC = ISO_STR2 + 1
          ISO2_QR = ISO_STR2 + 2
       else if(RAIN_TYPE=='COLD') then
          ISO_MAX = 12
          ISO_STR = TRC_VMAX + 1
          ISO_STR2 = TRC_VMAX + 7
          ISO_END = TRC_VMAX + ISO_MAX
          ISO1_QV = ISO_STR
          ISO1_QC = ISO_STR + 1
          ISO1_QR = ISO_STR + 2
          ISO1_QI = ISO_STR + 3
          ISO1_QS = ISO_STR + 4
          ISO1_QG = ISO_STR + 5
          ISO2_QV = ISO_STR2 
          ISO2_QC = ISO_STR2 + 1
          ISO2_QR = ISO_STR2 + 2
          ISO2_QI = ISO_STR2 + 3
          ISO2_QS = ISO_STR2 + 4
          ISO2_QG = ISO_STR2 + 5
       else if(RAIN_TYPE=='CLOUD_PARAM') then
          ISO_MAX = 4
          ISO_STR = TRC_VMAX + 1
          ISO_STR2 = TRC_VMAX + 3
          ISO_END = TRC_VMAX + ISO_MAX
          ISO1_QV = ISO_STR
          ISO1_QC = ISO_STR + 1
          ISO2_QV = ISO_STR2 
          ISO2_QC = ISO_STR2 + 1
       else if(RAIN_TYPE=='DRY') then
          ISO_MAX = 2
          ISO_STR = TRC_VMAX + 1
          ISO_STR2 = TRC_VMAX + 2
          ISO_END = TRC_VMAX + ISO_MAX
          ISO1_QV = ISO_STR
          ISO2_QV = ISO_STR2 
       end if
       TRC_VMAX = ISO_END
    end if
    !<= [add] K.Yoshimura 20110414

    allocate( TRC_name(TRC_vmax) ) ! [add] H.Yashiro 20110819
    allocate( WLABEL  (TRC_vmax) ) ! 08/04/12 [Add] T.Mitsui
    TRC_name(:) = ""
    WLABEL  (:) = ""

    allocate( flag_diagnose_number(TRC_vmax) )
    flag_diagnose_number(:) = .false.

    !--- Labeling
    do nq = 1, TRC_vmax
       if    ( nq == I_QV ) then
          TRC_name(nq) = 'qv'
          WLABEL  (nq) = 'VAPOR'
       elseif( nq == I_QC ) then
          TRC_name(nq) = 'qc'
          WLABEL  (nq) = 'CLOUD'
       elseif( nq == I_QR ) then
          TRC_name(nq) = 'qr'
          WLABEL  (nq) = 'RAIN'
       elseif( nq == I_QI ) then
          TRC_name(nq) = 'qi'
          WLABEL  (nq) = 'ICE'
       elseif( nq == I_QS ) then
          TRC_name(nq) = 'qs'
          WLABEL  (nq) = 'SNOW'
       elseif( nq == I_QG ) then
          TRC_name(nq) = 'qg'
          WLABEL  (nq) = 'GRAUPEL'

       elseif( nq == I_NC )then
          TRC_name(nq) = 'nc'
          WLABEL  (nq) = 'CLOUD_NUM'
       elseif( nq == I_NR )then
          TRC_name(nq) = 'nr'
          WLABEL  (nq) = 'RAIN_NUM'
       elseif( nq == I_NI )then
          TRC_name(nq) = 'ni'
          WLABEL  (nq) = 'ICE_NUM'
       elseif( nq == I_NS )then
          TRC_name(nq) = 'ns'
          WLABEL  (nq) = 'SNOW_NUM'
       elseif( nq == I_NG )then
          TRC_name(nq) = 'ng'
          WLABEL  (nq) = 'GRAUPEL_NUM'

       elseif( nq == NQDU_STR )then
          do i = 1, NDU_MAX
             write(is,'(I2.2)') i
             TRC_name(nq+i-1) = 'dust'//is
             WLABEL  (nq+i-1) = 'DUST'
          enddo
       elseif( nq == NQCB_STR ) then
          do i = 1, NCB_MAX
             write(is,'(I2.2)') i
             TRC_name(nq+i-1) = 'carbon'//is
             WLABEL  (nq+i-1) = 'CARBON'
          enddo
       elseif( nq == NQSU_STR ) then
          do i = 1, NSU_MAX
             write(is,'(I2.2)') i
             TRC_name(nq+i-1) = 'sulfate'//is
             WLABEL  (nq+i-1) = 'SULFATE'
          enddo
       elseif( nq == NQSA_STR ) then
          do i = 1, NSA_MAX
             write(is,'(I2.2)') i
             TRC_name(nq+i-1) = 'seasalt'//is
             WLABEL  (nq+i-1) = 'SALT'
          enddo

       elseif( nq == NQDUIN_STR )then
          do i = 1, NQDUIN_MAX
             write(is,'(I2.2)') i
             TRC_name(nq+i-1) = 'qcdust'//is
             WLABEL  (nq+i-1) = 'DUST_INCLOUD'
          enddo
       elseif( nq == NQCBIN_STR ) then
          do i = 1, NQCBIN_MAX
             write(is,'(I2.2)') i
             TRC_name(nq+i-1) = 'qccb'//is
             WLABEL  (nq+i-1) = 'CB_INCLOUD'
          enddo
       elseif( nq == NQSUIN_STR ) then
          TRC_name(nq) = 'qcso4'
          WLABEL  (nq) = 'SO4_INCLOUD'
       elseif( nq == NQSAIN_STR ) then
          do i = 1, NQSAIN_MAX
             write(is,'(I2.2)') i
             TRC_name(nq+i-1) = 'qcsalt'//is
             WLABEL  (nq+i-1) = 'SALT_INCLOUD'
          enddo
       !=> [add] K.Yoshimura 20110414
       elseif(nq == ISO1_QV)then
          TRC_name(nq) = 'qvo18'
          WLABEL  (nq) = 'O18_VAPOR'
       elseif(nq == ISO1_QC)then
          TRC_name(nq) = 'qco18'
          WLABEL  (nq) = 'O18_CLOUD'
       elseif(nq == ISO1_QR)then
          TRC_name(nq) = 'qro18'
          WLABEL  (nq) = 'O18_RAIN'
       elseif(nq == ISO1_QI)then
          TRC_name(nq) = 'qio18'
          WLABEL  (nq) = 'O18_ICE'
       elseif(nq == ISO1_QS)then
          TRC_name(nq) = 'qso18'
          WLABEL  (nq) = 'O18_SNOW'
       elseif(nq == ISO1_QG)then
          TRC_name(nq) = 'qgo18'
          WLABEL  (nq) = 'O18_GRAUPEL'
       elseif(nq == ISO2_QV)then
          TRC_name(nq) = 'qvhdo'
          WLABEL  (nq) = 'HDO_VAPOR'
       elseif(nq == ISO2_QC)then
          TRC_name(nq) = 'qchdo'
          WLABEL  (nq) = 'HDO_CLOUD'
       elseif(nq == ISO2_QR)then
          TRC_name(nq) = 'qrhdo'
          WLABEL  (nq) = 'HDO_RAIN'
       elseif(nq == ISO2_QI)then
          TRC_name(nq) = 'qihdo'
          WLABEL  (nq) = 'HDO_ICE'
       elseif(nq == ISO2_QS)then
          TRC_name(nq) = 'qshdo'
          WLABEL  (nq) = 'HDO_SNOW'
       elseif(nq == ISO2_QG)then
          TRC_name(nq) = 'qghdo'
          WLABEL  (nq) = 'HDO_GRAUPEL'
       !<= [add] K.Yoshimura 20110414
       elseif( nq == NCHEM_STR ) then
          do i = 1, NCHEM_MAX
             TRC_name(nq+i-1) = CHEM_TRC_name(i)
             WLABEL  (nq+i-1) = CHEM_TRC_desc(i)
          enddo
       endif
    enddo

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '--- Prognostic Tracers'
    write(ADM_LOG_FID,*) '|=========================================================|'
    write(ADM_LOG_FID,*) '|       :varname         :description                     |'
    do nq = 1, TRC_vmax
       write(ADM_LOG_FID,'(1x,A,I4,A,A16,A,A,A)') '|ID=', nq, ':', TRC_name(nq), ':', WLABEL(nq),'|'
    enddo
    write(ADM_LOG_FID,*) '|=========================================================|'

    !--- Heat capacity for thermodynamics
    allocate( CVW(NQW_STR:NQW_END) )
    allocate( CPW(NQW_STR:NQW_END) )

    if ( EIN_TYPE == 'EXACT' ) then
       LHV = CNST_LH00
       LHS = CNST_LHS00
       LHF = CNST_LHF00
       do nq = NQW_STR, NQW_END
          if    ( nq == I_QV ) then ! vapor
             CVW(nq) = CNST_CVV
             CPW(nq) = CNST_CPV
          elseif( nq == I_QC ) then ! cloud
             CVW(nq) = CNST_CL
             CPW(nq) = CNST_CL
          elseif( nq == I_QR ) then ! rain
             CVW(nq) = CNST_CL
             CPW(nq) = CNST_CL
          elseif( nq == I_QI ) then ! ice
             CVW(nq) = CNST_CI
             CPW(nq) = CNST_CI
          elseif( nq == I_QS ) then ! snow
             CVW(nq) = CNST_CI
             CPW(nq) = CNST_CI
          elseif( nq == I_QG ) then ! graupel
             CVW(nq) = CNST_CI
             CPW(nq) = CNST_CI
          endif
       enddo
    elseif( EIN_TYPE == 'SIMPLE' ) then
       LHV = CNST_LH0
       LHS = CNST_LHS0
       LHF = CNST_LHF0
       do nq = NQW_STR, NQW_END
          if    ( nq == I_QV ) then ! vapor
             CVW(nq) = CNST_CV
             CPW(nq) = CNST_CP
          elseif( nq == I_QC ) then ! cloud
             CVW(nq) = CNST_CV
             CPW(nq) = CNST_CV
          elseif( nq == I_QR ) then ! rain
             CVW(nq) = CNST_CV
             CPW(nq) = CNST_CV
          elseif( nq == I_QI ) then ! ice
             CVW(nq) = CNST_CV
             CPW(nq) = CNST_CV
          elseif( nq == I_QS ) then ! snow
             CVW(nq) = CNST_CV
             CPW(nq) = CNST_CV
          elseif( nq == I_QG ) then ! graupel
             CVW(nq) = CNST_CV
             CPW(nq) = CNST_CV
          endif
       enddo
    elseif( EIN_TYPE == 'SIMPLE2' ) then
       LHV = CNST_LH0
       LHS = CNST_LHS0
       LHF = CNST_LHF0
       do nq = NQW_STR, NQW_END
          if    ( nq == I_QV ) then ! vapor
             CVW(nq) = CNST_CVV
             CPW(nq) = CNST_CPV
          elseif( nq == I_QC ) then ! cloud
             CVW(nq) = CNST_CPV
             CPW(nq) = CNST_CPV
          elseif( nq == I_QR ) then ! rain
             CVW(nq) = CNST_CPV
             CPW(nq) = CNST_CPV
          elseif( nq == I_QI ) then ! ice
             CVW(nq) = CNST_CPV
             CPW(nq) = CNST_CPV
          elseif( nq == I_QS ) then ! snow
             CVW(nq) = CNST_CPV
             CPW(nq) = CNST_CPV
          elseif( nq == I_QG ) then ! graupel
             CVW(nq) = CNST_CPV
             CPW(nq) = CNST_CPV
          endif
       enddo
    endif

    if( mod(TB_DIV_NUM,SFC_DIV_NUM)/=0 ) then
       TB_DIV_NUM = SFC_DIV_NUM
       write(ADM_LOG_FID,*) '### FIX to TB_DIV_NUM =',TB_DIV_NUM
    endif

    if( THUBURN_LIM ) then  ![add] 20130613 R.Yoshida
       write(ADM_LOG_FID,*) 'Run with \"Thuburn Limiter\" in MIURA2004 Advection'
    else
       write(ADM_LOG_FID,*) '### Without \"Thuburn Limiter\" in MIURA2004 Advection'
    endif

    return
  end subroutine RUNCONF_setup

end module mod_runconf
!-------------------------------------------------------------------------------
