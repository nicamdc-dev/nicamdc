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

MPIEXEC="mpirun -np ${NMPI}"

outdir=${dir3d}
cd ${outdir}

cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Linux & intel fortran & intelmpi -----
#
################################################################################
ulimit -s 54975581388
export OMP_NUM_THREADS=1
export KMP_NUM_THREADS=1
export KMP_STACKSIZE=16m
#
echo "####################"
echo "### MPI NUM PROCESSES:" \${MPI_NUM_PROCESS}
echo "### Stack Size:" '`ulimit -s`'
echo "### OMP NUM THREADS:"   \${KMP_NUM_THREADS}
echo "### OMP STACK SIZE:"    \${KMP_STACKSIZE}
echo "####################"

ln -sv ../../../../bin/nhm_driver .
ln -sv ../../../../data/mnginfo/rl${RL}-prc${NP}.info .
ln -sv ../../../../data/grid/vgrid/${vgrid} .
EOF1

for f in $( ls ../../../../data/grid/boundary/${dir2d} )
do
   echo "ln -sv ../../../../data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

# run
${MPIEXEC} ./nhm_driver || exit

################################################################################
EOF2


cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR MacOSX & gfortran4.6 & OpenMPI1.6 -----
#
################################################################################
ulimit -s 54975581388
export OMP_NUM_THREADS=1
export KMP_NUM_THREADS=1
export KMP_STACKSIZE=16m
#
echo "####################"
echo "### MPI NUM PROCESSES:" \${MPI_NUM_PROCESS}
echo "### Stack Size:" '`ulimit -s`'
echo "### OMP NUM THREADS:"   \${KMP_NUM_THREADS}
echo "### OMP STACK SIZE:"    \${KMP_STACKSIZE}
echo "####################"

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
