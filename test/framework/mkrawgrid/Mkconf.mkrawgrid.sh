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
/

&PARAM_MKGRD
    MKGRD_DOSPRING         = .true.,
    MKGRD_OUT_BASENAME     = "rawgrid_${res2d}",
    MKGRD_spring_beta      = 1.15D0,
/

################################################################################
EOFNHM
