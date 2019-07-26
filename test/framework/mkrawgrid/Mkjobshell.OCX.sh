#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpiexec.hydra -n \${PJM_MPI_PROC}"

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

MNGINFO=rl${RL}-prc${NP}.info

NNODE=`expr \( $NMPI - 1 \) / 28 + 1`

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Oakbridge-CX -----
#
################################################################################
#PJM -g gn11
#PJM -L rscgrp=regular
#PJM -L node=${NNODE}
#PJM --mpi proc=${NMPI}
#PJM -L elapse=06:00:00
#PJM -N NICAMDC
#PJM -j
#PJM -s
#
export FORT_FMT_RECL=400
export OMP_NUM_THREADS=1

module load hdf5
module load netcdf
module load netcdf-fortran


ln -svf ${TOPDIR}/bin/${BINNAME} .
EOF1

if   [ -f ${TOPDIR}/data/mnginfo/${MNGINFO} ]; then
   echo "mnginfo file is found in default database"
   echo "ln -svf ${TOPDIR}/data/mnginfo/${MNGINFO} ." >> run.sh
elif [ -f ../../mkmnginfo/rl${RL}pe${NP}/${MNGINFO} ]; then
   echo "mnginfo file is found in test directory"
   echo "ln -svf ../../mkmnginfo/rl${RL}pe${NP}/${MNGINFO} ." >> run.sh
else
   echo "mnginfo file is not found!"
   exit 1
fi

cat << EOF2 >> run.sh

# run
${MPIEXEC} ./${BINNAME} mkrawgrid.cnf || exit

################################################################################
EOF2
