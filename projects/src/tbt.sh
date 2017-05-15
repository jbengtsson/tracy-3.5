#!/bin/sh

gcc tbt.cc -o tbt \
    -L$NUM_REC/lib -lnum_rec \
    -lstdc++ \
    -lm
