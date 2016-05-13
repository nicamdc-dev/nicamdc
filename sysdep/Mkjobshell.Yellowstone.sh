#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpirun.lsf"

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

# for DCMIP2016
prjcode="SCIS0006"

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# for NCAR Yellowstone (IBM iDataPlex Sandybridge)
#
################################################################################
#BSUB -a poe                  # set parallel operating environment
#BSUB -P ${prjcode}           # project code
#BSUB -J nicamdc              # job name
#BSUB -W 00:10                # wall-clock time (hrs:mins)
#BSUB -n ${NMPI}              # number of tasks in job
#BSUB -R "span[ptile=4]"      # run four MPI tasks per node
#BSUB -q regular              # queue
#BSUB -e errors.%J.nicamdc    # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.nicamdc    # output file name in which %J is replaced by the job ID
 
export OMP_NUM_THREADS=4
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS
 
${MPIEXEC} ./${BINNAME}

################################################################################
EOF1


cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# for NCAR Yellowstone (IBM iDataPlex Sandybridge)
#
################################################################################
#BSUB -a poe                  # set parallel operating environment
#BSUB -P ${prjcode}           # project code
#BSUB -J ico2ll               # job name
#BSUB -W 00:10                # wall-clock time (hrs:mins)
#BSUB -n ${NMPI}              # number of tasks in job
#BSUB -R "span[ptile=4]"      # run four MPI tasks per node
#BSUB -q regular              # queue
#BSUB -e errors.%J.ico2ll     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.ico2ll     # output file name in which %J is replaced by the job ID
 
export OMP_NUM_THREADS=4
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

# run
${MPIEXEC} ./fio_ico2ll_mpi \
history \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="./${MNGINFO}" \
layerfile_dir="./." \
llmap_base="./llmap" \
outfile_dir="../" \
-lon_swap \
-comm_smallchunk

################################################################################
EOFICO2LL1
