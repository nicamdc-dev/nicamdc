#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="mpirun -np ${NMPI}"

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

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & intel mpi -----
#
################################################################################
export FORT_FMT_RECL=400



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
