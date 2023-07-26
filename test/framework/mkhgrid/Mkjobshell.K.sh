#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpiexec"

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
res2d=GL${GL}RL${RL}
BNDDIR="${TOPDIR}/data/grid/boundary"

MNGINFO=rl${RL}-prc${NP}.info

# for K computer
if [ ${NMPI} -gt 36864 ]; then
   rscgrp="huge"
elif [ ${NMPI} -gt 384 ]; then
   rscgrp="large"
else
   rscgrp="small"
fi

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ For K computer
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${NMPI}"
#PJM --rsc-list "elapse=00:30:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ${TOPDIR}/bin/${BINNAME} %r:./"
#PJM --stgin  "rank=* ./mkhgrid.cnf            %r:./"
EOF1

if   [ -f ${TOPDIR}/data/mnginfo/${MNGINFO} ]; then
   echo "mnginfo file is found in default database"
   echo "#PJM --stgin  \"rank=* ${TOPDIR}/data/mnginfo/${MNGINFO} %r:./\"" >> run.sh
elif [ -f ../../mkmnginfo/rl${RL}pe${NP}/${MNGINFO} ]; then
   echo "mnginfo file is found in test directory"
   echo "#PJM --stgin  \"rank=* ../../mkmnginfo/rl${RL}pe${NP}/${MNGINFO} %r:./\"" >> run.sh
else
   echo "mnginfo file is not found!"
   exit 1
fi

if ls ../../mkrawgrid/${dir2d}/rawgrid_${res2d}.pe* > /dev/null 2>&1
then
   echo "#PJM --stgin  \"rank=* ../../mkrawgrid/${dir2d}/rawgrid_${res2d}.pe%06r %r:./\"" >> run.sh
else
   echo "rawgrid file is not found!"
   exit 1
fi

cat << EOF2 >> run.sh
#PJM --stgout "rank=* %r:./boundary*.pe%06r       ./${BNDDIR}/${dir2d}/"
#PJM --stgout "rank=* %r:./*                      ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2

# run
${MPIEXEC} ./${BINNAME} mkhgrid.cnf || exit

################################################################################
EOF2
