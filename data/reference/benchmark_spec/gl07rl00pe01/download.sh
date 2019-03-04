#! /bin/bash

# please get from http://r-ccs-climate.riken.jp/members/yashiro/spec_nicamdc/data/reference/benchmark_spec/gl07rl00pe01/restart_all_GL05RL00z94.pe000000

dir=data/reference/benchmark_spec/gl07rl00pe01
file=restart_all_GL07RL00z94.pe000000

rm -f ./$file
wget http://r-ccs-climate.riken.jp/members/yashiro/spec_nicamdc/$dir/$file || exit 1

echo "Checking md5sum:"
md5sum -c ${file}.md5 || exit 1

