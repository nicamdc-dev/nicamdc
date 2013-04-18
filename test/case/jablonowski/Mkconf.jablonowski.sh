#! /bin/bash -x

##### set parameters from ARGV
GLEV=${1:-5}
RLEV=${2:-0}
NMPI=${3:-5}
GL=`printf %02d ${GLEV}`
RL=`printf %02d ${RLEV}`
if   [ ${NMPI} -ge 10000 ]; then
	NP=`printf %05d ${NMPI}`
elif [ ${NMPI} -ge 1000 ]; then
	NP=`printf %04d ${NMPI}`
elif [ ${NMPI} -ge 100 ]; then
	NP=`printf %03d ${NMPI}`
else
	NP=`printf %02d ${NMPI}`
fi
ZL=${4:-40}
vgrid=${5:-vgrid40_24000-600m.dat}
LSMAX=${6:-0}
DTL=${7:-1200}
DIFCF=${8:-1.29D16}
NHIST=${9:-0}

if [ ${LSMAX} == 0 ]; then
   # 11day
   let LSMAX=" 11 * 24 * 60 * 60 / ${DTL} "
fi

if [ ${NHIST} == 0 ]; then
   # 1day
   let NHIST="  1 * 24 * 60 * 60 / ${DTL} "
fi

dir2d=gl${GL}rl${RL}pe${NP}
res2d=GL${GL}RL${RL}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res3d=GL${GL}RL${RL}z${ZL}

let ZLp2="${ZL} + 2"

outdir=${dir3d}
cd ${outdir}

##### Generate nhm_driver.cnf
cat << EOFNHM > nhm_driver.cnf
################################################################################
#
# NICAM driver config (Jablonowski Baroclinic Wave)
#
################################################################################

###--- ADMIN SETUP---###
&ADMPARAM
    glevel      = ${GLEV},
    rlevel      = ${RLEV},
    vlayer      = ${ZL},
    rgnmngfname = "rl${RL}-prc${NP}.info",
/

###--- SET Grid system, topography, vegetation type
&GRDPARAM
    hgrid_io_mode   = "ADVANCED",
    hgrid_fname     = "boundary_${res2d}",
    vgrid_fname     = "${vgrid}",
    topo_fname      = "Jablonowski",
/

###--- for FULL RUN (11 days): LSTEP_MAX = 792
&TIMEPARAM
    DTL         = ${DTL}.D0,
    INTEG_TYPE  = "RK2",
    LSTEP_MAX   = ${LSMAX},
    SSTEP_MAX   = 4,
    SPLIT       = .true.,
    start_year  = 1000,
    start_month = 1,
    start_day   = 1,
/

###--- SET RUN TYPE
&RUNCONFPARAM 
    RUN_TYPE       = 'Jablonowski',
    EIN_TYPE       = 'SIMPLE',
    NDIFF_LOCATION = 'IN_LARGE_STEP',
    AF_TYPE        = 'NONE',
/

###--- SET BASIC STATE
&BSSTATEPARAM
    ref_type = 'NOBASE',
/

&NUMFILTERPARAM
    dep_hgrid         = .false.,
    lap_order_divdamp = 2,
    hdiff_type        = 'DIRECT',
    Kh_coef_maxlim    = 9.9D99,
    Kh_coef_minlim    = 0.D0,
    divdamp_type      = 'DIRECT',
    alpha_d           = ${DIFCF},
    alpha_dv          = 0.D0,
    gamma_h           = ${DIFCF},
    gamma_v           = 0.D0,
    hdiff_fact_rho    = 0.01D0,
    hdiff_fact_q      = 0.D0,
    ZD                = 99999.9D0,
    alpha_r           = 0.D0,
    DEEP_EFFECT       = .false.,
/

&RESTARTPARAM
    input_io_mode     = 'IDEAL',
    output_io_mode    = 'ADVANCED',
    output_basename   = 'restart_all_${res3d}',
    restart_layername = 'ZSALL${ZLp2}',
/

&DYCORETESTPARAM
    init_type = 'Jablonowski',
    test_case = '1'
/

&EMBUDGETPARAM MNT_ON = .true., MNT_INTV = ${NHIST} /

&NMHISD
    output_io_mode    = 'ADVANCED' ,
    histall_fname     = 'history'  ,
    hist3D_layername  = 'ZSDEF${ZL}',
    DIRECT_ACCESS     = .true.     ,
    NO_VINTRPL        = .false.    ,
    output_type       = 'SNAPSHOT' ,
    step              = ${NHIST}   ,
    doout_step0       = .true.     ,
/

&NMHIST item='ml_u',     file='u',   ktype='3D' /
&NMHIST item='ml_v',     file='v',   ktype='3D' /
&NMHIST item='ml_w',     file='w',   ktype='3D' /
&NMHIST item='ml_pres',  file='prs', ktype='3D' /
&NMHIST item='ml_tem',   file='t',   ktype='3D' /
&NMHIST item='sl_ps',    file='ps',  ktype='2D' /

################################################################################
EOFNHM
