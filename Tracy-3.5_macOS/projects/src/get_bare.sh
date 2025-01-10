#!/bin/sh

prm1=$1

#queue="ap-medium.q"
queue="ap-high.q"
#queue="test-medium.q"

t1="00:10:00"
t2="00:30:00"

#dir=`pwd`
dir=$HOME/git_repos/tracy-3.5/projects/src

\rm dnu.cmd.o* fmap*.cmd.o*

qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue -v flat_file=$prm1 $dir/dnu.cmd
qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue -v flat_file=$prm1 $dir/fmap.cmd
qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue -v flat_file=$prm1 $dir/fmap_dp.cmd
