#! /bin/bash -x

GLEV=${1:-5}
RLEV=${2:-0}
NMPI=${3:-5}
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
ZL=${4:-40}
vgrid=${5:-vgrid40_24000.dat}
case=${6:-"heldsuarez"}

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=gl${GL}rl${RL}

outdir=${NICAMRUNDIR}/${case}/${dir3d}
mkdir -p ${outdir}

cat << EOF1 > ${outdir}/run.sh
#! /bin/bash -x
################################################################################
#
# for K Computer
#
################################################################################
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=${NP}"
#PJM --rsc-list "elapse=00:10:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ${NICAMBINDIR}/nhm_driver %r:./"
#PJM --stgin  "rank=* ${outdir}/nhm_driver.cnf %r:./"
#PJM --stgin  "rank=* ${NICAMDIR}/data/mnginfo/rl${RL}-prc${NP}.info %r:./"
#PJM --stgin  "rank=* ${NICAMDIR}/data/grid/vgrid/${vgrid} %r:./"
#PJM --stgin  "rank=* ${NICAMDIR}/data/grid/boundary/${dir2d}/boundary_${res2d}.pe%06r %r:./"
#PJM --stgout "rank=* %r:./* ${outdir}/"
#PJM --stgout "rank=* %r:./prof/* ${outdir}/prof/"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=\${PARALLEL}
export LPG="/opt/FJSVxosmmm/sbin/lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export FLIB_FASTOMP=TRUE
export THREAD_STACK_SIZE=8192
export XOS_MMM_L_ARENA_FREE=2

# run
export fprof="fipp -C -Ihwm -d prof"
rm -rf ./prof
# run
\${fprof} ${MPIEXEC} ./nhm_driver

################################################################################
EOF1


cat << EOFICO2LL1 > ${outdir}/ico2ll.sh
#! /bin/bash -x
################################################################################
#
# for K Computer (frontend)
#
################################################################################
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:10:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ${NICAMBINDIR}/ico2ll %r:./"
#PJM --stgin  "rank=* ${outdir}/ico2ll.cnf %r:./"
#PJM --stgin  "rank=* ${outdir}/history.info %r:./"
#PJM --stgin  "rank=* ${outdir}/*.rgn* %r:./"
#PJM --stgin  "rank=* ${NICAMDIR}/data/grid/llmap/gl${GL}/rl${RL}/* %r:./"
#PJM --stgout "rank=* %r:./* ${outdir}/"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=\${PARALLEL}
export LPG="/opt/FJSVxosmmm/sbin/lpgparm -s 32MB -d 32MB -h 32MB -t 32MB -p 32MB"
export FLIB_FASTOMP=TRUE
export THREAD_STACK_SIZE=8192
export XOS_MMM_L_ARENA_FREE=2

# run
./ico2ll

################################################################################
EOFICO2LL1
