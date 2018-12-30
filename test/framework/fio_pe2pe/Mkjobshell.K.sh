#! /bin/bash -x

GLEV=${1}
RLEV_I=${2}
NMPI_I=${3}
RLEV_O=${4}
NMPI_O=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpiexec"

GL=`printf %02d ${GLEV}`
RL_I=`printf %02d ${RLEV_I}`
RL_O=`printf %02d ${RLEV_O}`
if   [ ${NMPI_I} -ge 10000 ]; then
	NP_I=`printf %05d ${NMPI_I}`
elif [ ${NMPI_I} -ge 1000 ]; then
	NP_I=`printf %04d ${NMPI_I}`
elif [ ${NMPI_I} -ge 100 ]; then
	NP_I=`printf %03d ${NMPI_I}`
else
	NP_I=`printf %02d ${NMPI_I}`
fi
if   [ ${NMPI_O} -ge 10000 ]; then
	NP_O=`printf %05d ${NMPI_O}`
elif [ ${NMPI_O} -ge 1000 ]; then
	NP_O=`printf %04d ${NMPI_O}`
elif [ ${NMPI_O} -ge 100 ]; then
	NP_O=`printf %03d ${NMPI_O}`
else
	NP_O=`printf %02d ${NMPI_O}`
fi

dir2d_I=gl${GL}rl${RL_I}pe${NP_I}
res2d_I=GL${GL}RL${RL_I}
dir2d_O=gl${GL}rl${RL_O}pe${NP_O}
res2d_O=GL${GL}RL${RL_O}

MNGINFO_I=rl${RL_I}-prc${NP_I}.info
MNGINFO_O=rl${RL_O}-prc${NP_O}.info

# for K computer
if [ ${NMPI_I} -gt 36864 ]; then
   rscgrp="huge"
elif [ ${NMPI_I} -gt 384 ]; then
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
#PJM --rsc-list "node=${NMPI_I}"
#PJM --rsc-list "elapse=01:00:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ${TOPDIR}/bin/${BINNAME}                            %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/data/mnginfo/${MNGINFO_I}                 %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/data/mnginfo/${MNGINFO_O}                 %r:./"
#PJM --stgin  "rank=* ${TOPDIR}/data/grid/boundary/${dir2d}/boundary_${res2d}.pe%06r %r:./"
#PJM --stgout "rank=* %r:./* ./"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2

# run
mkdir -p ./${dir2d_O}
rm -f dump*
${MPIEXEC} ./${BINNAME} boundary_${res2d_I} outfile="./${dir2d_O}/boundary_${res2d_O}" \
           glevel=${GLEV} \
           rlevel_in=${RLEV_I}  mnginfo_in="./${MNGINFO_I}"  \
           rlevel_out=${RLEV_O} mnginfo_out="./${MNGINFO_O}" || exit
rm -f dump*

################################################################################
EOF1
