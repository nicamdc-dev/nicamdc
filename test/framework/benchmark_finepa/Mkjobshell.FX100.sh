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

NNODE=`expr $NMPI / 2`

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# for FX100
#
################################################################################
#PJM -L "rscgrp=all"
#PJM -L node=${NNODE}
#PJM --mpi proc=${NMPI}
#PJM -L elapse=01:00:00
#PJM -j
#PJM -s
#
. /fefs/home/system/Env_base
#
export PARALLEL=16
export OMP_NUM_THREADS=16
export XOS_MMM_L_ARENA_FREE=2

module load TCSuite
module load HDF5
module load NetCDF-C
module load NetCDF-Fortran

ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d} )
do
   echo "ln -sv ${TOPDIR}/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

# run
for PA_ID in `seq 1 11`; do
   echo "rm -rf pa${PA_ID}" >> run.sh
   echo "fapp -C -d pa${PA_ID} -Hpa=${PA_ID} -Hmethod=raw ${MPIEXEC} -n ${NMPI} ./${BINNAME} || exit" >> run.sh
done

cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# for FX100
#
################################################################################
#PJM -L "rscgrp=all"
#PJM -L node=${NNODE}
#PJM --mpi proc=${NMPI}
#PJM -L elapse=01:00:00
#PJM -j
#PJM -s
#
. /fefs/home/system/Env_base
#
export PARALLEL=16
export OMP_NUM_THREADS=16
export XOS_MMM_L_ARENA_FREE=2

module load TCSuite
module load HDF5
module load NetCDF-C
module load NetCDF-Fortran

ln -sv ${TOPDIR}/bin/fio_ico2ll_mpi .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/zaxis .
EOFICO2LL1

for f in $( ls ${TOPDIR}/data/grid/llmap/gl${GL}rl${RL}/ )
do
   echo "ln -sv ${TOPDIR}/data/grid/llmap/gl${GL}rl${RL}/${f} ." >> ico2ll.sh
done

cat << EOFICO2LL2 >> ico2ll.sh

# run
${MPIEXEC} ./fio_ico2ll_mpi \
history \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="./${MNGINFO}" \
layerfile_dir="./zaxis" \
llmap_base="./llmap" \
-lon_swap \
-comm_smallchunk

################################################################################
EOFICO2LL2
