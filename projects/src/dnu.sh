#!/bin/sh

\rm nohup.out

DIR=$HOME/git_repos/tracy-3.5/projects/src

nohup $DIR/dnu      flat_file.dat &
nohup $DIR/ptc/fmap flat_file.dat 1 &
nohup $DIR/ptc/fmap flat_file.dat 2 &

# FILE_DIR=$HOME/git_repos/thor-2.0/thor/wrk

# nohup $DIR/dnu      $FILE_DIR/flat_file.fit &
# nohup $DIR/ptc/fmap $FILE_DIR/flat_file.fit 1 &
# nohup $DIR/ptc/fmap $FILE_DIR/flat_file.fit 2 &

wait
