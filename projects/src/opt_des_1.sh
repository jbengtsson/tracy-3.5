#!/bin/sh

prm1=$1

#queue="ap-medium.q"
queue="ap-high.q"
#queue="test-medium.q"

t1="00:10:00"
t2="00:30:00"

#dir=`pwd`
dir=$HOME/git_repos/tracy-3.5_temp/projects/src

\rm opt_des_1.cmd.o*

qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue -v lat_file=$prm1 $dir/opt_des_1.cmd
