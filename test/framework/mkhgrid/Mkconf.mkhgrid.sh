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

outdir=${dir3d}
cd ${outdir}

##### Generate config file
cat << EOFNHM > ${BINNAME}.cnf
################################################################################
#
# NICAM tool config
#
################################################################################

&ADMPARAM
    glevel      = ${GLEV},
    rlevel      = ${RLEV},
    vlayer      = 1,
    rgnmngfname = "${MNGINFO}",
!    ADM_HGRID_SYSTEM = "1DMD-ON-SPHERE",
/

&PARAM_MKGRD
    MKGRD_DOPREROTATE      = .true.,
    MKGRD_DOSTRETCH        = .false.,
    MKGRD_DOSHRINK         = .true.,
    MKGRD_DOROTATE         = .true.,
    MKGRD_IN_BASENAME      = "rawgrid_${res2d}",
    MKGRD_OUT_BASENAME     = "hgrid_${res2d}",
    MKGRD_prerotation_tilt =  90.D0,
    MKGRD_stretch_alpha    =   4.D0,
    MKGRD_shrink_level     =      1,
    MKGRD_rotation_lon     = 140.D0,
    MKGRD_rotation_lat     =  40.D0,
/

################################################################################
EOFNHM
