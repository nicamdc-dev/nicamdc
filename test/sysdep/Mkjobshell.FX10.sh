#! /bin/bash -x

##### set parameters from ARGV
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
vgrid=${5:-vgrid40_24000-600m.dat}
LSMAX=${6:-792}
DTL=${7:-1200}
DIFCF=${8:-1.29D16}
NHIST=${9:-72}

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}

MPIEXEC="mpiexec"

outdir=${dir3d}
cd ${outdir}

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# for OAKLEAF-FX
#
################################################################################
#PJM --rsc-list "rscgrp=short"
#PJM --rsc-list "node=${NMPI}"
#PJM --rsc-list "elapse=01:00:00"
#PJM -j
#PJM -s
export PARALLEL=8
export OMP_NUM_THREADS=8
export fu30bf=10

ln -sv ../../../../bin/nhm_driver .
ln -sv ../../../../data/mnginfo/rl${RL}-prc${NP}.info .
ln -sv ../../../../data/grid/vgrid/${vgrid} .
EOF1

for f in $( ls ../../../../data/grid/boundary/${dir2d} )
do
   echo "ln -sv ../../../../data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

fprof="fipp -C -Srange -Ihwm -d prof"
rm -rf ./prof

# run
\${fprof} ${MPIEXEC} ./nhm_driver || exit
# ${MPIEXEC} ./nhm_driver || exit

################################################################################
EOF2


cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# for OAKLEAF-FX
#
################################################################################
#PJM --rsc-list "rscgrp=short"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:30:00"
#PJM -j
#PJM -s
export PARALLEL=8
export OMP_NUM_THREADS=8
export fu08bf=10

ln -sv ../../../../bin/fio_ico2ll_mpi .
ln -sv ../../../../data/mnginfo/rl${RL}-prc${NP}.info .
ln -sv ../../../../data/zaxis .
EOFICO2LL1

for f in $( ls ../../../../data/grid/llmap/gl${GL}/rl${RL}/ )
do
   echo "ln -sv ../../../../data/grid/llmap/gl${GL}/rl${RL}/${f} ." >> ico2ll.sh
done

cat << EOFICO2LL2 >> ico2ll.sh

# run
${MPIEXEC} ./fio_ico2ll_mpi \
history \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="./rl${RL}-prc${NP}.info" \
layerfile_dir="./zaxis" \
llmap_base="./llmap" \
-lon_swap \
-comm_smallchunk

################################################################################
EOFICO2LL2
