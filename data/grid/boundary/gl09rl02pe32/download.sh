#! /bin/bash

# please get from http://r-ccs-climate.riken.jp/members/yashiro/spec_nicamdc/data/grid/boundary/gl09rl02pe32/

dir=data/grid/boundary/gl09rl02pe32

for prc in `seq 0 31`
do
   PE=`printf %06d ${prc}`

   file=boundary_GL09RL02.pe${PE}

   rm -f ./$file
   wget http://r-ccs-climate.riken.jp/members/yashiro/spec_nicamdc/$dir/$file || exit 1

   echo "Checking md5sum:"
   md5sum -c ${file}.md5 || exit 1
done
