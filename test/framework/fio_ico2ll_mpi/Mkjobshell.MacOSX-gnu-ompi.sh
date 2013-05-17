#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}

# System specific
MPIEXEC="openmpirun -np ${NMPI}"

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

outdir=${dir3d}
cd ${outdir}

cat << EOFICO2LL1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR MacOSX & gfortran4.6 & OpenMPI1.6 -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ../../mkmnginfo/${dir3d}/${MNGINFO} .
ln -sv ${TOPDIR}/data/zaxis .
EOFICO2LL1

for f in $( ls ../../mkllmap/${dir3d}/sample_${res2d}.pe* )
do
   echo "ln -sv ${f} ." >> run.sh
done

for f in $( ls ../../mkllmap/${dir3d}/llmap.* )
do
   echo "ln -sv ${f} ." >> run.sh
done

cat << EOFICO2LL2 >> run.sh

# run
${MPIEXEC} ./fio_ico2ll_mpi \
sample_${res2d} \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="./${MNGINFO}" \
layerfile_dir="./zaxis" \
llmap_base="./llmap" \
nlim_llgrid=1000000 \
-lon_swap \
-comm_smallchunk

################################################################################
EOFICO2LL2
