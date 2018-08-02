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
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=GL${GL}RL${RL}
res3d=GL${GL}RL${RL}z${ZL}

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
# for K computer
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${NMPI}"
#PJM --rsc-list "elapse=01:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ${TOPDIR}/bin/${BINNAME}                            %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/data/mnginfo/${MNGINFO}                   %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/data/zaxis/*                              %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/data/grid/llmap/gl${GL}rl${RL}/llmap.*    %r:./"
#PJM --stgin  "rank=* ../../mkllmap/gl${GL}rl${RL}/sample_${res2d}.pe%06r %r:./"
#PJM --stgout "rank=* %r:./llmap.*     ${TOPDIR}/data/grid/llmap/gl${GL}rl${RL}/"
#PJM --stgout "rank=* %r:./*           ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2

# run
${MPIEXEC} ./${BINNAME} sample_${res2d} \
           glevel=${GLEV} \
           rlevel=${RLEV} \
           mnginfo="./${MNGINFO}" \
           layerfile_dir="./zaxis" \
           llmap_base="./llmap" \
           -lon_swap -output_netcdf || exit

################################################################################
EOF1
