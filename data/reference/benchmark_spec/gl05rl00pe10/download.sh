#! /bin/bash

# please get from http://r-ccs-climate.riken.jp/members/yashiro/spec_nicamdc/data/reference/benchmark_spec/gl05rl00pe10/restart_all_GL05RL00z94.pe000000

dir=data/reference/benchmark_spec/gl05rl00pe10

for prc in `seq 0 9`
do
   PE=`printf %06d ${prc}`

   file=restart_all_GL05RL00z94.pe${PE}

   rm -f ./$file
   wget http://r-ccs-climate.riken.jp/members/yashiro/spec_nicamdc/$dir/$file || exit 1
done

echo "Checking md5sum:"
md5sum -c restart_all_GL05RL00z94.md5 || exit 1

