#!/bin/sh

dir=`pwd`

TRACY="tracy-3.5"

cd "$TRACY_LIB"

rm -rf autom4te.cache
rm -rf aclocal.m4
rm -rf tracy/lib/*

make distclean

./bootstrap
./configure --prefix=$TRACY_LIB/tracy

make install
