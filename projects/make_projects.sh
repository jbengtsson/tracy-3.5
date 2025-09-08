#!/bin/sh

dir=`pwd`

rm -rf autom4te.cache
rm -rf aclocal.m4

# Create config dir – if it doesn't exist.
mkdir -p config

make distclean

./bootstrap
./configure --prefix=$dir/projects

make
