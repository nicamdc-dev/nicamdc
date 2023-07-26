#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=1
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

NNODE=`expr \( $NMPI - 1 \) / 64 + 1`
NPROC=`expr $NMPI / $NNODE`
NPIND=`expr \( 255 \) / $NPROC + 1`

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Oakforest-PACS -----
#
################################################################################
#PJM -g jh180023
#PJM -L rscgrp=regular-cache
#PJM -L node=${NNODE}
#PJM -L elapse=00:30:00
#PJM -N NICAMDC
#PJM -j
#PJM -s
#
module load hdf5_szip
module load hdf5
module load netcdf
module load netcdf-fortran

export FORT_FMT_RECL=400

ln -svf ${TOPDIR}/bin/${BINNAME} .

# run
./${BINNAME} mkmnginfo.cnf || exit
mkdir -p     ${TOPDIR}/data/mnginfo
mv -f *.info ${TOPDIR}/data/mnginfo

################################################################################
EOF1
