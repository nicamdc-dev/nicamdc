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
   for nmpi in 4 8 10 20 40
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
   for nmpi in 16 32 80 160
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

for glevel in 9
do

   rlevel=2
   for nmpi in 40 80 160
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

      outdir=gl${GL}rl02pe32_rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile_pe32 |
      sed -e "s/#rlevel_out#/${rlevel}/g"                       |
      sed -e "s/#nmpi_out#/${nmpi}/g"                           > ${outdir}/Makefile
   done

   rlevel=3
   for nmpi in 128 256 320 640
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

      outdir=gl${GL}rl02pe32_rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile_pe32 |
      sed -e "s/#rlevel_out#/${rlevel}/g"                       |
      sed -e "s/#nmpi_out#/${nmpi}/g"                           > ${outdir}/Makefile
   done

   rlevel=4
   for nmpi in 512 1280 2560
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

      outdir=gl${GL}rl02pe32_rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile_pe32 |
      sed -e "s/#rlevel_out#/${rlevel}/g"                       |
      sed -e "s/#nmpi_out#/${nmpi}/g"                           > ${outdir}/Makefile
   done

   rlevel=5
   for nmpi in 1024 2048 5120 10240
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

      outdir=gl${GL}rl02pe32_rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile_pe32 |
      sed -e "s/#rlevel_out#/${rlevel}/g"                       |
      sed -e "s/#nmpi_out#/${nmpi}/g"                           > ${outdir}/Makefile
   done

   rlevel=6
   for nmpi in 4096 8192 20480 40960
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

      outdir=gl${GL}rl02pe32_rl${RL}pe${PE}
      mkdir -p ${outdir}

      sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile_pe32 |
      sed -e "s/#rlevel_out#/${rlevel}/g"                       |
      sed -e "s/#nmpi_out#/${nmpi}/g"                           > ${outdir}/Makefile
   done

done
