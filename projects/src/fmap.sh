#!/bin/sh

DIR=$HOME/git_repos/projects/src

nohup $DIR/leac param.dat 1 >& leac_1.log &
nohup $DIR/leac param.dat 2 >& leac_2.log &
