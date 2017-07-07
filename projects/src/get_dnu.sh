#!/bin/sh

#queue="ap-medium.q"
#queue="ap-high.q"
queue="test-medium.q"

t1="00:10:00"
t2="00:30:00"

#dir=`pwd`
dir=$HOME/git_repos/tracy-3.5/projects/src

#~/projects/src/main /home/bengtsson/projects/in/lattice/sls-2
#wait

\rm dnu.cmd.o* fmap*.cmd.o*

qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue $dir/dnu.cmd
qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue $dir/fmap.cmd
qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue $dir/fmap_dp.cmd
