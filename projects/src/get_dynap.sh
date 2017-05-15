#!/bin/sh

queue="prime_bdl.q"

t1="96:00:00"
t2="96:00:00"

#dir=`pwd`
dir=$HOME/git_repos/projects/src

\rm dynap.cmd.o*

qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue $dir/dynap.cmd
