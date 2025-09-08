#!/bin/sh

dir=`pwd`

rm -rf autom4te.cache
rm -rf aclocal.m4
rm -rf tracy/lib/*

make distclean

# Configure libtool (for shared libraries).
# Linux.
# libtoolize
# Macbook.
glibtoolize

# Create config dir – if it doesn't exist.
mkdir -p config

./bootstrap
./configure --prefix=$dir/tracy

make install
