#! /bin/bash -x

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

      outdir=gl${GL}rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile      |
      sed -e "s/#rlevel#/${rlevel}/g"                           |
      sed -e "s/#nmpi#/${nmpi}/g"                               |
      sed -e "s/#rl#/rl${RL}/g"                                 |
      sed -e "s/#prc#/prc${PE}/g"                               |
      sed -e "s/#GL#/GL${GL}/g"                                 |
      sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/mkrawgrid.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                           |
      sed -e "s/#nmpi#/${nmpi}/g"                               |
      sed -e "s/#rl#/rl${RL}/g"                                 |
      sed -e "s/#prc#/prc${PE}/g"                               |
      sed -e "s/#GL#/GL${GL}/g"                                 |
      sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/mkrawgrid.cnf
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

      outdir=gl${GL}rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile      |
      sed -e "s/#rlevel#/${rlevel}/g"                           |
      sed -e "s/#nmpi#/${nmpi}/g"                               |
      sed -e "s/#rl#/rl${RL}/g"                                 |
      sed -e "s/#prc#/prc${PE}/g"                               |
      sed -e "s/#GL#/GL${GL}/g"                                 |
      sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/mkrawgrid.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                           |
      sed -e "s/#nmpi#/${nmpi}/g"                               |
      sed -e "s/#rl#/rl${RL}/g"                                 |
      sed -e "s/#prc#/prc${PE}/g"                               |
      sed -e "s/#GL#/GL${GL}/g"                                 |
      sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/mkrawgrid.cnf
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

      outdir=gl${GL}rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile      |
      sed -e "s/#rlevel#/${rlevel}/g"                           |
      sed -e "s/#nmpi#/${nmpi}/g"                               |
      sed -e "s/#rl#/rl${RL}/g"                                 |
      sed -e "s/#prc#/prc${PE}/g"                               |
      sed -e "s/#GL#/GL${GL}/g"                                 |
      sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/Makefile

      sed -e "s/#glevel#/${glevel}/g" ./templates/mkrawgrid.cnf |
      sed -e "s/#rlevel#/${rlevel}/g"                           |
      sed -e "s/#nmpi#/${nmpi}/g"                               |
      sed -e "s/#rl#/rl${RL}/g"                                 |
      sed -e "s/#prc#/prc${PE}/g"                               |
      sed -e "s/#GL#/GL${GL}/g"                                 |
      sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/mkrawgrid.cnf
   done

done
