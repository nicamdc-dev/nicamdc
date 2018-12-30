#! /bin/bash -x

for glevel in 5 6 7
do

   rlevel=0
   nmpi=10
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

   outdir=gl${GL}rl${RL}
   mkdir -p ${outdir}

   sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile      |
   sed -e "s/#rlevel#/${rlevel}/g"                           |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               |
   sed -e "s/#GL#/GL${GL}/g"                                 |
   sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/Makefile

   sed -e "s/#glevel#/${glevel}/g" ./templates/mkllmap.cnf |
   sed -e "s/#rlevel#/${rlevel}/g"                           |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               |
   sed -e "s/#GL#/GL${GL}/g"                                 |
   sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/mkllmap.cnf

   rlevel=1
   nmpi=40
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

   outdir=gl${GL}rl${RL}
   mkdir -p ${outdir}

   sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile      |
   sed -e "s/#rlevel#/${rlevel}/g"                           |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               |
   sed -e "s/#GL#/GL${GL}/g"                                 |
   sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/Makefile

   sed -e "s/#glevel#/${glevel}/g" ./templates/mkllmap.cnf |
   sed -e "s/#rlevel#/${rlevel}/g"                           |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               |
   sed -e "s/#GL#/GL${GL}/g"                                 |
   sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/mkllmap.cnf

   rlevel=2
   nmpi=160
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

   outdir=gl${GL}rl${RL}
   mkdir -p ${outdir}

   sed -e "s/#glevel#/${glevel}/g" ./templates/Makefile      |
   sed -e "s/#rlevel#/${rlevel}/g"                           |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               |
   sed -e "s/#GL#/GL${GL}/g"                                 |
   sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/Makefile

   sed -e "s/#glevel#/${glevel}/g" ./templates/mkllmap.cnf |
   sed -e "s/#rlevel#/${rlevel}/g"                           |
   sed -e "s/#nmpi#/${nmpi}/g"                               |
   sed -e "s/#rl#/rl${RL}/g"                                 |
   sed -e "s/#prc#/prc${PE}/g"                               |
   sed -e "s/#GL#/GL${GL}/g"                                 |
   sed -e "s/#RL#/RL${RL}/g"                                 > ${outdir}/mkllmap.cnf

done
