#!/bin/sh

queue="prime_bdl.q"

t1="24:00:00"
t2="96:00:00"

#dir=`pwd`
dir=$HOME/git_repos/projects/src

\rm leac*.cmd.o*

qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue $dir/leac_1.cmd
qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue $dir/leac_2.cmd
