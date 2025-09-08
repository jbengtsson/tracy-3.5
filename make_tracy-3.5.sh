#!/bin/sh

dir=`pwd`

rm -rf autom4te.cache
rm -rf aclocal.m4
rm -rf tracy/lib/*

# Create config dir – if it doesn't exist.
mkdir -p config

make distclean

# Configure libtool (for shared libraries).
case "$(uname -s)" in
    Linux*)
        echo "Running on Linux"
        libtoolize
        ;;
    Darwin*)
	echo "Running on macOS"
        glibtoolize
	;;
esac

./bootstrap
./configure --prefix=$dir/tracy

make install
