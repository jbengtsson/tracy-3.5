#!/bin/sh

dir=`pwd`

rm -rf autom4te.cache
rm -rf aclocal.m4

./bootstrap
./configure --prefix=$dir

make
