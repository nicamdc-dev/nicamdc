#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpiexec.hydra -n ${NMPI}"

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

NNODE=`expr \( $NMPI - 1 \) / 32 + 1`
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
#PJM --mpi proc=${NMPI}
#PJM --omp thread=2
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

export HFI_NO_CPUAFFINITY=1
export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_POLICY=hbw_bind,,
export I_MPI_FABRICS_LIST=tmi
unset KMP_AFFINITY
#export KMP_AFFINITY=verbose
#export I_MPI_DEBUG=5

export OMP_NUM_THREADS=2
export I_MPI_PIN_DOMAIN=${NPIND}
export I_MPI_PERHOST=${NPROC}
export KMP_HW_SUBSET=1t
export I_MPI_FABRICS=shm:tmi
export I_MPI_HARD_FINALIZE=1

ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d} )
do
   echo "ln -sv ${TOPDIR}/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

# run
${MPIEXEC} ./${BINNAME} || exit

################################################################################
EOF2


cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Oakforest-PACS -----
#
################################################################################
export FORT_FMT_RECL=400


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
