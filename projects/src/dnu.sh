#!/bin/sh

DIR=$HOME/git_repos/projects/src

FILE_DIR=$HOME/git_repos/thor-2.0/thor/wrk

\rm nohup.out

# nohup $DIR/dnu      flat_file.dat &
# nohup $DIR/ptc/fmap flat_file.dat 1 &
# nohup $DIR/ptc/fmap flat_file.dat 2 &

nohup $DIR/dnu      $FILE_DIR/flat_file.fit &
nohup $DIR/ptc/fmap $FILE_DIR/flat_file.fit 1 &
nohup $DIR/ptc/fmap $FILE_DIR/flat_file.fit 2 &

wait
