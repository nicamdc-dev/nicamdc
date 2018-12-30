#! /bin/bash -x

for glevel in 5 6 7
do

   rlevel=0
   for nmpi in 2 5 10
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

      outdir=gl${GL}rl00pe01_rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile |
      sed -e "s/#rlevel_out#/${rlevel}/g"                  |
      sed -e "s/#nmpi_out#/${nmpi}/g"                      > ${outdir}/Makefile
   done

   rlevel=1
#   for nmpi in 1 2 4 5 8 10 20 40
   for nmpi in 4 8 20 40
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

      outdir=gl${GL}rl00pe01_rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile |
      sed -e "s/#rlevel_out#/${rlevel}/g"                  |
      sed -e "s/#nmpi_out#/${nmpi}/g"                      > ${outdir}/Makefile
   done

   rlevel=2
#   for nmpi in 1 2 4 5 8 10 16 20 32 40 80 160
   for nmpi in 1 2 4 16 20 32 40 80 160
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

      outdir=gl${GL}rl00pe01_rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile |
      sed -e "s/#rlevel_out#/${rlevel}/g"                  |
      sed -e "s/#nmpi_out#/${nmpi}/g"                      > ${outdir}/Makefile
   done

done
