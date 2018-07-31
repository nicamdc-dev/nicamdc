#! /bin/bash -x

rlevel=0
for nmpi in 1 2 5 10
do
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

   outdir=rl${RL}pe${PE}
   mkdir -p ${outdir}

   cp ./templates/Makefile ${outdir}/Makefile

   sed -e "s/#rlevel#/${rlevel}/g" ./templates/mkmnginfo.cnf |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               > ${outdir}/mkmnginfo.cnf
done

rlevel=1
for nmpi in 1 2 4 5 8 10 20 40
do
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

   outdir=rl${RL}pe${PE}
   mkdir -p ${outdir}

   cp ./templates/Makefile ${outdir}/Makefile

   sed -e "s/#rlevel#/${rlevel}/g" ./templates/mkmnginfo.cnf |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               > ${outdir}/mkmnginfo.cnf
done

rlevel=2
for nmpi in 1 2 4 5 8 10 16 20 32 40 80 160
do
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

   outdir=rl${RL}pe${PE}
   mkdir -p ${outdir}

   cp ./templates/Makefile ${outdir}/Makefile

   sed -e "s/#rlevel#/${rlevel}/g" ./templates/mkmnginfo.cnf |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               > ${outdir}/mkmnginfo.cnf
done
