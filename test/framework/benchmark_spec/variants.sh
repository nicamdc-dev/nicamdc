#! /bin/bash -x

# resolution dependent parameters
#              GL04    GL05    GL06    GL07    GL08    GL09
DTLLIST=(   2400.D0 1200.D0  600.D0  300.D0  180.D0   30.D0 )
LSTEPLIST=(     252     504    1008    2016    3360    1440 )
NDIFFLIST=( 1.0D+17 1.3D+16 1.6D+15 2.0D+14 2.5D+13 3.2D+12 )

for glevel in 5 6 7
do

   rlevel=0
   for nmpi in 1 2 5 10
   do
      GL=`printf %02d ${glevel}`
      RL=`printf %02d ${rlevel}`
      if   [ ${nmpi} -ge 10000 ]; then
         PE=`printf %05d ${nmpi}`
      elif [ ${nmpi} -ge 1000 ]; then
         PE=`printf %04d ${nmpi}`
      elif [ ${nmpi} -ge 100 ]; then
         PE=`printf %03d ${nmpi}`
      else
         PE=`printf %02d ${nmpi}`
      fi

      let idx="${glevel} - 4"
      DTL=${DTLLIST[${idx}]}
      LSTEP=${LSTEPLIST[${idx}]}
      NDIFF=${NDIFFLIST[${idx}]}

      outdir=gl${GL}rl${RL}z94pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile       |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/nhm_driver.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  |
      sed -e "s/#DTL#/${DTL}/g"                                  |
      sed -e "s/#LSTEP#/${LSTEP}/g"                                  |
      sed -e "s/#NDIFF#/${NDIFF}/g"                              > ${outdir}/nhm_driver.cnf
   done

   rlevel=1
   for nmpi in 1 2 4 5 8 10 20 40
   do
      GL=`printf %02d ${glevel}`
      RL=`printf %02d ${rlevel}`
      if   [ ${nmpi} -ge 10000 ]; then
         PE=`printf %05d ${nmpi}`
      elif [ ${nmpi} -ge 1000 ]; then
         PE=`printf %04d ${nmpi}`
      elif [ ${nmpi} -ge 100 ]; then
         PE=`printf %03d ${nmpi}`
      else
         PE=`printf %02d ${nmpi}`
      fi

      let idx="${glevel} - 4"
      DTL=${DTLLIST[${idx}]}
      LSTEP=${LSTEPLIST[${idx}]}
      NDIFF=${NDIFFLIST[${idx}]}

      outdir=gl${GL}rl${RL}z94pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile       |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/nhm_driver.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  |
      sed -e "s/#DTL#/${DTL}/g"                                  |
      sed -e "s/#LSTEP#/${LSTEP}/g"                                  |
      sed -e "s/#NDIFF#/${NDIFF}/g"                              > ${outdir}/nhm_driver.cnf
   done

   rlevel=2
   for nmpi in 1 2 4 5 8 10 16 20 32 40 80 160
   do
      GL=`printf %02d ${glevel}`
      RL=`printf %02d ${rlevel}`
      if   [ ${nmpi} -ge 10000 ]; then
         PE=`printf %05d ${nmpi}`
      elif [ ${nmpi} -ge 1000 ]; then
         PE=`printf %04d ${nmpi}`
      elif [ ${nmpi} -ge 100 ]; then
         PE=`printf %03d ${nmpi}`
      else
         PE=`printf %02d ${nmpi}`
      fi

      let idx="${glevel} - 4"
      DTL=${DTLLIST[${idx}]}
      LSTEP=${LSTEPLIST[${idx}]}
      NDIFF=${NDIFFLIST[${idx}]}

      outdir=gl${GL}rl${RL}z94pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile       |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/nhm_driver.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  |
      sed -e "s/#DTL#/${DTL}/g"                                  |
      sed -e "s/#LSTEP#/${LSTEP}/g"                                  |
      sed -e "s/#NDIFF#/${NDIFF}/g"                              > ${outdir}/nhm_driver.cnf
   done

done

for glevel in 9
do

   rlevel=2
   for nmpi in 32
   do
      GL=`printf %02d ${glevel}`
      RL=`printf %02d ${rlevel}`
      if   [ ${nmpi} -ge 10000 ]; then
         PE=`printf %05d ${nmpi}`
      elif [ ${nmpi} -ge 1000 ]; then
         PE=`printf %04d ${nmpi}`
      elif [ ${nmpi} -ge 100 ]; then
         PE=`printf %03d ${nmpi}`
      else
         PE=`printf %02d ${nmpi}`
      fi

      let idx="${glevel} - 4"
      DTL=${DTLLIST[${idx}]}
      LSTEP=${LSTEPLIST[${idx}]}
      NDIFF=${NDIFFLIST[${idx}]}

      outdir=gl${GL}rl${RL}z94pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile       |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/nhm_driver.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  |
      sed -e "s/#DTL#/${DTL}/g"                                  |
      sed -e "s/#LSTEP#/${LSTEP}/g"                                  |
      sed -e "s/#NDIFF#/${NDIFF}/g"                              > ${outdir}/nhm_driver.cnf
   done

   rlevel=5
   for nmpi in 2048
   do
      GL=`printf %02d ${glevel}`
      RL=`printf %02d ${rlevel}`
      if   [ ${nmpi} -ge 10000 ]; then
         PE=`printf %05d ${nmpi}`
      elif [ ${nmpi} -ge 1000 ]; then
         PE=`printf %04d ${nmpi}`
      elif [ ${nmpi} -ge 100 ]; then
         PE=`printf %03d ${nmpi}`
      else
         PE=`printf %02d ${nmpi}`
      fi

      let idx="${glevel} - 4"
      DTL=${DTLLIST[${idx}]}
      LSTEP=${LSTEPLIST[${idx}]}
      NDIFF=${NDIFFLIST[${idx}]}

      outdir=gl${GL}rl${RL}z94pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile       |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/nhm_driver.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  |
      sed -e "s/#DTL#/${DTL}/g"                                  |
      sed -e "s/#LSTEP#/${LSTEP}/g"                                  |
      sed -e "s/#NDIFF#/${NDIFF}/g"                              > ${outdir}/nhm_driver.cnf
   done

   rlevel=6
   for nmpi in 8192
   do
      GL=`printf %02d ${glevel}`
      RL=`printf %02d ${rlevel}`
      if   [ ${nmpi} -ge 10000 ]; then
         PE=`printf %05d ${nmpi}`
      elif [ ${nmpi} -ge 1000 ]; then
         PE=`printf %04d ${nmpi}`
      elif [ ${nmpi} -ge 100 ]; then
         PE=`printf %03d ${nmpi}`
      else
         PE=`printf %02d ${nmpi}`
      fi

      let idx="${glevel} - 4"
      DTL=${DTLLIST[${idx}]}
      LSTEP=${LSTEPLIST[${idx}]}
      NDIFF=${NDIFFLIST[${idx}]}

      outdir=gl${GL}rl${RL}z94pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile       |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/nhm_driver.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                            |
      sed -e "s/#nmpi#/${nmpi}/g"                                |
      sed -e "s/#rl#/rl${RL}/g"                                  |
      sed -e "s/#prc#/prc${PE}/g"                                |
      sed -e "s/#GL#/GL${GL}/g"                                  |
      sed -e "s/#RL#/RL${RL}/g"                                  |
      sed -e "s/#DTL#/${DTL}/g"                                  |
      sed -e "s/#LSTEP#/${LSTEP}/g"                                  |
      sed -e "s/#NDIFF#/${NDIFF}/g"                              > ${outdir}/nhm_driver.cnf
   done

done
