#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

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

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=GL${GL}RL${RL}
res3d=GL${GL}RL${RL}z${ZL}

MNGINFO=rl${RL}-prc${NP}.info

LSMAX=${8:-0}
DTL=${9:-1200}
DIFCF=${10:-1.29D16}
NHIST=${11:-0}

if [ ${LSMAX} == 0 ]; then
   # 11day
   let LSMAX=" 11 * 24 * 60 * 60 / ${DTL} "
fi

if [ ${NHIST} == 0 ]; then
   # 1day
   let NHIST="  1 * 24 * 60 * 60 / ${DTL} "
fi

let ZLp2="${ZL} + 2"

outdir=${dir3d}
cd ${outdir}

##### Generate config file
cat << EOFNHM > ${BINNAME}.cnf
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
    VGRID_fname     = "${VGRID}",
    topo_fname      = "Jablonowski",
/

###--- for FULL RUN (11 days): LSTEP_MAX = 792
&TIMEPARAM
    DTL         = ${DTL}.D0,
    INTEG_TYPE  = "RK3",
    LSTEP_MAX   = ${LSMAX},
    SSTEP_MAX   = 6,
    SPLIT       = .true.,
    start_year  = 1000,
    start_month = 1,
    start_day   = 1,
/

###--- SET RUN TYPE
&RUNCONFPARAM 
    RUN_TYPE       = 'Jablonowski',
    EIN_TYPE       = 'SIMPLE',
    NDIFF_LOCATION = 'IN_LARGE_STEP2',
    AF_TYPE        = 'NONE',
/

###--- SET BASIC STATE
&BSSTATEPARAM
    ref_type = 'NOBASE',
/

&NUMFILTERPARAM
    hdiff_type        = 'DIRECT',
    lap_order_hdiff   = 2,
    gamma_h           = ${DIFCF},
    Kh_coef_minlim    = 0.D0,
    Kh_coef_maxlim    = 9.9D99,
    divdamp_type      = 'DIRECT',
    lap_order_divdamp = 2,
    alpha_d           = ${DIFCF},
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

&NMHIST item='ml_u',     file='u',    ktype='3D' /
&NMHIST item='ml_v',     file='v',    ktype='3D' /
&NMHIST item='ml_w',     file='w',    ktype='3D' /
&NMHIST item='ml_pres',  file='prs',  ktype='3D' /
&NMHIST item='ml_tem',   file='t',    ktype='3D' /
&NMHIST item='sl_ps',    file='ps',   ktype='2D' /
&NMHIST item='sl_u850',  file='u850', ktype='2D' /
&NMHIST item='sl_v850',  file='v850', ktype='2D' /
&NMHIST item='sl_w850',  file='w850', ktype='2D' /
&NMHIST item='sl_t850',  file='t850', ktype='2D' /


################################################################################
EOFNHM
