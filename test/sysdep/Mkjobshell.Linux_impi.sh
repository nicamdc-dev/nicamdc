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
LSMAX=${6:-792}
DTL=${7:-1200}
DIFCF=${8:-1.29D16}
NHIST=${9:-72}
case=${10:-"heldsuarez"}

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}

outdir=${NICAMRUNDIR}/${case}/${dir3d}
mkdir -p ${outdir}

cat << EOF1 > ${outdir}/run.sh
#! /bin/bash -x
################################################################################
#
# for Linux(Westmere) & ifort & intelmpi
#
################################################################################
ulimit -s 54975581388
export KMP_NUM_THREADS=1
export OMP_NUM_THREADS=1
export KMP_STACKSIZE=16m
#
echo "####################"
echo "### MPI NUM PROCESSES:" \${MPI_NUM_PROCESS}
echo "### Stack Size:" '`ulimit -s`'
echo "### OMP NUM THREADS:"   \${KMP_NUM_THREADS}
echo "### OMP STACK SIZE:"    \${KMP_STACKSIZE}
echo "####################"

ln -sv ${NICAMBINDIR}/nhm_driver .
ln -sv ${NICAMDIR}/data/mnginfo/rl${RL}-prc${NP}.info .
ln -sv ${NICAMDIR}/data/grid/vgrid/${vgrid} .
EOF1

for f in $( ls ${NICAMDIR}/data/grid/boundary/${dir2d} )
do
   echo "ln -sv ${NICAMDIR}/data/grid/boundary/${dir2d}/${f} ." >> ${outdir}/run.sh
done

cat << EOF2 >> ${outdir}/run.sh

# run
echo "### JOB start time:"
date
${MPIEXEC} -np ${NMPI} ./nhm_driver
echo "### JOB end time:"
date
echo "####################"

################################################################################
EOF2


cat << EOFICO2LL1 > ${outdir}/ico2ll.sh
#! /bin/bash -x
################################################################################
#
# for Linux(Westmere) & ifort & intelmpi
#
################################################################################
ulimit -s 54975581388
export KMP_NUM_THREADS=1
export OMP_NUM_THREADS=1
export KMP_STACKSIZE=16m
#
echo "####################"
echo "### MPI NUM PROCESSES:" \${MPI_NUM_PROCESS}
echo "### Stack Size:" '`ulimit -s`'
echo "### OMP NUM THREADS:"   \${KMP_NUM_THREADS}
echo "### OMP STACK SIZE:"    \${KMP_STACKSIZE}
echo "####################"

ln -sv ${NICAMBINDIR}/fio_ico2ll_mpi .
ln -sv ${NICAMDIR}/data/mnginfo/rl${RL}-prc${NP}.info .
ln -sv ${NICAMDIR}/data/zaxis .
EOFICO2LL1

for f in $( ls ${NICAMDIR}/data/grid/llmap/gl${GL}/rl${RL}/ )
do
   echo "ln -sv ${NICAMDIR}/data/grid/llmap/gl${GL}/rl${RL}/${f} ." >> ${outdir}/ico2ll.sh
done

cat << EOFICO2LL2 >> ${outdir}/ico2ll.sh

# run
${MPIEXEC} -np ${NMPI} ./fio_ico2ll_mpi \
history \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="rl${RL}-prc${NP}.info" \
layerfile_dir="./zaxis" \
llmap_base="./llmap" \
nlim_llgrid=20000 \
-lon_swap \
-comm_smallchunk

# draw figures
ln -s ../../../test/grads/jab_ps.gs ./
ln -s ../../../test/grads/jab_t.gs ./
ln -s ../../../test/grads/mul.gs ./
ln -s ../../../test/grads/cbarn.gs ./
../../../test/grads/bin/grads -blc "run jab_ps.gs"
../../../test/grads/bin/grads -blc "run jab_t.gs"

################################################################################
EOFICO2LL2
