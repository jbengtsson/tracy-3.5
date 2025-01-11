#!/bin/sh

prm1=$1

queue="ap-high.q"
#queue="ap-medium.q"
#queue="test-medium.q"

t1="00:10:00"
t2="12:00:00"

dir=$HOME/git_repos/tracy-3.5_temp/projects/src

\rm bare_da.cmd.o*

qsub -l h_rt=$t2 -q $queue -v flat_file=$prm1 $dir/bare_da.cmd
