#!/bin/sh

NUM_REC=$HOME/num_rec
TRACY_LIB=$HOME/git_repos/tracy-3.5

g++ -c -Wno-non-template-friend -fPIC \
    -I$TRACY_LIB/tracy/inc \
    -lstdc++ -lm -lgfortran \
    -o tracy_2_3.o tracy_2_3.cc

g++ -Wall -Wno-non-template-friend -fPIC \
    -o tracy_2_3 tracy_2_3.o \
    -L$TRACY_LIB/tracy/lib -L$NUM_REC/lib \
    -ltracy_ptc -lTPSALib -lLieLib -lnum_rec \
    -lstdc++ -lm -lgfortran

