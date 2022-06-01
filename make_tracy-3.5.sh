#!/bin/sh

#echo $PATH
#PATH=$PATH:/net/home/lli/OPAL/OPAL-2.0.1/bin
#echo $PATH
#export PATH

dir=`pwd`

cd "$TRACY_LIB"

rm -rf autom4te.cache
rm -rf aclocal.m4
rm -rf tracy/lib/*

make distclean

./bootstrap
./configure --prefix=$TRACY_LIB/tracy

make install
