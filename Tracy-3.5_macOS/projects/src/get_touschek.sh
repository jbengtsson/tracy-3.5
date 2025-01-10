#!/bin/sh

queue="ap-high.q"

t1="96:00:00"
t2="96:00:00"

dir=$TRACY_LIB/projects/src

\rm touschek.cmd.o*

qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue $dir/touschek.cmd
