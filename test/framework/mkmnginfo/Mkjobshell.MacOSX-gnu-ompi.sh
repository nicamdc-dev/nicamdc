#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ For MacOSX & gfortran7.3 & OpenMPI3.0 -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

ln -svf ${TOPDIR}/bin/${BINNAME} .

# run
./${BINNAME} mkmnginfo.cnf || exit
mkdir -p     ${TOPDIR}/data/mnginfo
mv -f *.info ${TOPDIR}/data/mnginfo

################################################################################
EOF1
