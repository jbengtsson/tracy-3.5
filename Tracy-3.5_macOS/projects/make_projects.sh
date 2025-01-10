#!/bin/sh

dir=`pwd`

rm -rf autom4te.cache
rm -rf aclocal.m4

make distclean

./bootstrap
./configure --prefix=$dir/projects

make
