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

OMP_NUM_THREADS=1

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
#BSUB -R "span[ptile=10]"     # run four MPI tasks per node
#BSUB -q regular              # queue
#BSUB -e errors.%J.nicamdc    # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.nicamdc    # output file name in which %J is replaced by the job ID
 
export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export MP_TASK_AFFINITY=core:${OMP_NUM_THREADS}
module load mkl/10.3.11

ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d} )
do
   echo "ln -sv ${TOPDIR}/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

${MPIEXEC} ./${BINNAME}

################################################################################
EOF2


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
#BSUB -R "span[ptile=10]"     # run four MPI tasks per node
#BSUB -q regular              # queue
#BSUB -e errors.%J.ico2ll     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.ico2ll     # output file name in which %J is replaced by the job ID
 
export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export MP_TASK_AFFINITY=core:${OMP_NUM_THREADS}
module load mkl/10.3.11

ln -sv ${TOPDIR}/bin/fio_ico2ll_mpi .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/zaxis .
EOFICO2LL1

for f in $( ls ${TOPDIR}/data/grid/llmap/gl${GL}/rl${RL}/ )
do
   echo "ln -sv ${TOPDIR}/data/grid/llmap/gl${GL}/rl${RL}/${f} ." >> ico2ll.sh
done

cat << EOFICO2LL2 >> ico2ll.sh

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
EOFICO2LL2


cat << EOFICO2LLNC1 > ico2ll_netcdf.sh
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
#BSUB -R "span[ptile=10]"     # run four MPI tasks per node
#BSUB -q regular              # queue
#BSUB -e errors.%J.ico2ll     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.ico2ll     # output file name in which %J is replaced by the job ID

export OMP_NUM_THREADS=1
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS
module load mkl/10.3.11
module load cdo/1.6.3

# User Settings
# ------------------------------------------------------------------------------

glev=5          # g-level of original grid
case=163        # test case number
out_intev='1hr' # output interval (format: "1hr", "6hr", "day", "100s")

# ------------------------------------------------------------------------------

ln -sv ${TOPDIR}/bin/fio_ico2ll_mpi .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/zaxis .
EOFICO2LLNC1

for f in $( ls ${TOPDIR}/data/grid/llmap/gl${GL}/rl${RL}/ )
do
   echo "ln -sv ${TOPDIR}/data/grid/llmap/gl${GL}/rl${RL}/${f} ." >> ico2ll_netcdf.sh
done

cat << EOFICO2LLNC2 >> ico2ll_netcdf.sh

# run
${MPIEXEC} ./fio_ico2ll_mpi \
history \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="./${MNGINFO}" \
layerfile_dir="./zaxis" \
llmap_base="./llmap" \
output_netcdf=.true. \
dcmip2016=.true. \
-lon_swap \
-comm_smallchunk

# netcdf combine
case \${glev} in
5)
 reso='r200'
;;
6)
 reso='r100'
;;
7)
 reso='r50'
;;
8)
 reso='r25'
;;
esac

case \${case} in
161) # Moist baroclinic wave
levs='L30'
target=( u.nc v.nc w.nc prs.nc t.nc qv.nc qc.nc qr.nc pasv1.nc pasv2.nc \
u850.nc v850.nc w850.nc t850.nc w500.nc t500.nc ps.nc prcp.nc \
forcing_vx.nc forcing_vy.nc forcing_vz.nc forcing_e.nc \
forcing_qv.nc forcing_qc.nc forcing_qr.nc forcing_cl.nc forcing_cl2.nc \
cl_column.nc cl2_column.nc cly_column.nc )
;;
162) # Idealized tropical cyclone
levs='L30'
target=( u.nc v.nc w.nc prs.nc t.nc qv.nc qc.nc qr.nc \
u850.nc v850.nc w850.nc t850.nc w500.nc t500.nc ps.nc prcp.nc \
forcing_vx.nc forcing_vy.nc forcing_vz.nc forcing_e.nc \
forcing_qv.nc forcing_qc.nc forcing_qr.nc )
;;
163) # Supercell
levs='L40'
target=( u.nc v.nc w.nc prs.nc t.nc qv.nc qc.nc qr.nc \
u850.nc v850.nc w850.nc t850.nc w500.nc t500.nc ps.nc prcp.nc \
forcing_vx.nc forcing_vy.nc forcing_vz.nc forcing_e.nc \
forcing_qv.nc forcing_qc.nc forcing_qr.nc )
;;
esac

model='nicam'
org_grid='hex'
dat_grid='interp_latlon'
equation='nonhydro'
description='g-level_'\${glev}

fname=\${model}.\${case}.\${reso}.\${levs}.\${dat_grid}.\${equation}.nc

input=""
for (( i=0; i <\${#target[@]}; ++i ))
do
input=\$input" "\${target[\$i]}
done

echo "start cdo process"
cdo -s merge \${input}                                           temporary_A.nc
cdo -s setgatt,model,\${model}                                   temporary_A.nc temporary_B.nc
rm temporary_A.nc; cdo -s setgatt,test_case,\${case}             temporary_B.nc temporary_A.nc
rm temporary_B.nc; cdo -s setgatt,horizontal_resolution,\${reso} temporary_A.nc temporary_B.nc
rm temporary_A.nc; cdo -s setgatt,levels,\${levs}                temporary_B.nc temporary_A.nc
rm temporary_B.nc; cdo -s setgatt,grid,\${dat_grid}              temporary_A.nc temporary_B.nc
rm temporary_A.nc; cdo -s setgatt,equation,\${equation}          temporary_B.nc temporary_A.nc
rm temporary_B.nc; cdo -s setgatt,time_frequency,\${out_intev}   temporary_A.nc temporary_B.nc
rm temporary_A.nc; cdo -s setgatt,description,\${description}    temporary_B.nc temporary_A.nc
cdo -s setgatt,history,dcmip2016                                temporary_A.nc \${fname}

for (( i=0; i <\${#target[@]}; ++i ))
do
 rm \${target[\$i]}
done
rm temporary_A.nc temporary_B.nc

################################################################################
EOFICO2LLNC2
