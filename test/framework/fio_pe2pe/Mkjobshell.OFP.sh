#! /bin/bash -x

GLEV=${1}
RLEV_I=${2}
NMPI_I=${3}
RLEV_O=${4}
NMPI_O=${5}
TOPDIR=${6}
BINNAME=${7}

NMPI=1 # serial execution

# System specific
MPIEXEC="mpiexec.hydra -n ${NMPI}"

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
res3d_I=GL${GL}RL${RL_I}z94
dir2d_O=gl${GL}rl${RL_O}pe${NP_O}
res2d_O=GL${GL}RL${RL_O}
res3d_O=GL${GL}RL${RL_O}z94

MNGINFO_I=rl${RL_I}-prc${NP_I}.info
MNGINFO_O=rl${RL_O}-prc${NP_O}.info

NNODE=`expr \( $NMPI - 1 \) / 64 + 1`
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
#PJM --omp thread=1
#PJM -L elapse=05:00:00
#PJM -N pe2pe_gl${GL}
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

export OMP_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=${NPIND}
export I_MPI_PERHOST=${NPROC}
export KMP_HW_SUBSET=1t
export I_MPI_FABRICS=shm:tmi
export I_MPI_HARD_FINALIZE=1


ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO_I} .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO_O} .
EOF1

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d_I} )
do
   echo "ln -sv ${TOPDIR}/data/grid/boundary/${dir2d_I}/${f} ." >> run.sh
done

for f in $( ls ${TOPDIR}/data/reference/benchmark_spec/${dir2d_I} )
do
   echo "ln -sv ${TOPDIR}/data/reference/benchmark_spec/${dir2d_I}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

# run
mkdir -p ./${dir2d_O}
rm -f dump*
${MPIEXEC} ./${BINNAME} boundary_${res2d_I} outfile="./${dir2d_O}/boundary_${res2d_O}" \
           glevel=${GLEV} \
           rlevel_in=${RLEV_I}  mnginfo_in="./${MNGINFO_I}"  \
           rlevel_out=${RLEV_O} mnginfo_out="./${MNGINFO_O}" || exit
rm -f dump*
${MPIEXEC} ./${BINNAME} restart_all_${res3d_I} outfile="./${dir2d_O}/restart_all_${res3d_O}" \
           glevel=${GLEV} \
           rlevel_in=${RLEV_I}  mnginfo_in="./${MNGINFO_I}"  \
           rlevel_out=${RLEV_O} mnginfo_out="./${MNGINFO_O}" || exit
rm -f dump*
mkdir -p ${TOPDIR}/data/grid/boundary/${dir2d_O}
mv -f ./${dir2d_O}/boundary_*    ${TOPDIR}/data/grid/boundary/${dir2d_O}/
mkdir -p ${TOPDIR}/data/reference/benchmark_spec/${dir2d_O}
mv -f ./${dir2d_O}/restart_all_* ${TOPDIR}/data/reference/benchmark_spec/${dir2d_O}/

################################################################################
EOF2
