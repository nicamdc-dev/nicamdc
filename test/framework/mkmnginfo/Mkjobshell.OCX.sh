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
# ------ FOR Oakbridge-CX -----
#
################################################################################
#PJM -g gn11
#PJM -L rscgrp=regular
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=00:30:00
#PJM -N NICAMDC
#PJM -X
#PJM -j
#PJM -s
#
export FORT_FMT_RECL=400
export OMP_NUM_THREADS=1

module load hdf5
module load netcdf
module load netcdf-fortran


ln -svf ${TOPDIR}/bin/${BINNAME} .

# run
./${BINNAME} mkmnginfo.cnf || exit
mkdir -p     ${TOPDIR}/data/mnginfo
mv -f *.info ${TOPDIR}/data/mnginfo

################################################################################
EOF1
