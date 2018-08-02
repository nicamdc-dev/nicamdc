#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpirun --mca btl openib,self"

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

NNODE=`expr \( $NMPI - 1 \) / 20 + 1`
NPROC=`expr $NMPI / $NNODE`

if [ ${NNODE} -gt 16 ]; then
   rscgrp="l"
elif [ ${NNODE} -gt 3 ]; then
   rscgrp="m"
else
   rscgrp="s"
fi

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & mpt & torque -----
#
################################################################################
#PBS -q ${rscgrp}
#PBS -l nodes=${NNODE}:ppn=${NPROC}
#PBS -N ico2ll_gl${GL}rl${RL}
#PBS -l walltime=05:00:00
#PBS -o STDOUT
#PBS -e STDERR
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

cd \$PBS_O_WORKDIR

ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/zaxis .
EOF1

if ls ../../mkllmap/gl${GL}rl${RL}/sample_${res2d}.pe* > /dev/null 2>&1
then
   for f in $( ls ../../mkllmap/gl${GL}rl${RL}/sample_${res2d}.pe* )
   do
      echo "ln -sv ../../mkllmap/gl${GL}rl${RL}/${f} ." >> run.sh
   done
else
   echo "sample data in mkllmap is not found!"
   exit 1
fi

if ls ${TOPDIR}/data/grid/llmap/gl${GL}rl${RL}/llmap.info > /dev/null 2>&1
then
   for f in $( ls ${TOPDIR}/data/grid/llmap/gl${GL}rl${RL}/ )
   do
      echo "ln -sv ${TOPDIR}/data/grid/llmap/gl${GL}rl${RL}/${f} ." >> run.sh
   done
else
   echo "llmap data is not found!"
   exit 1
fi

cat << EOF2 >> run.sh

# run
${MPIEXEC} ./${BINNAME} sample_${res2d} \
           glevel=${GLEV} \
           rlevel=${RLEV} \
           mnginfo="./${MNGINFO}" \
           layerfile_dir="./zaxis" \
           llmap_base="./llmap" \
           -lon_swap -output_netcdf || exit

################################################################################
EOF2
