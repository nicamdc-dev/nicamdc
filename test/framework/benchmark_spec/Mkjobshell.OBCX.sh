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

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=GL${GL}RL${RL}
res3d=GL${GL}RL${RL}z${ZL}

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
#PJM -L elapse=00:30:00
#PJM -N NICAMDC
#PJM -j
#PJM -s
#
export FORT_FMT_RECL=400
export OMP_NUM_THREADS=1

module load hdf5
module load netcdf
module load netcdf-fortran
module unload impi/2019.4.243
module unload intel/2019.4.243
module load intel/2018.3.222
module load impi/2018.3.222


ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d} )
do
   echo "ln -sv ${TOPDIR}/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

for f in $( ls ${TOPDIR}/data/reference/benchmark_spec/${dir2d} )
do
   echo "ln -sv ${TOPDIR}/data/reference/benchmark_spec/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

# run
${MPIEXEC} ./${BINNAME} || exit

################################################################################
EOF2
