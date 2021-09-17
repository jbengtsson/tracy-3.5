#!/bin/sh
#------------------------------------------------------------------------------
# Set up of Johan's MAC environment
# example for anyone that needs a dedicated setup

# Set up of required programs
CELLAR_DIR=/usr/local/Cellar/gcc/11.1.0_1/bin/
export PATH=$CELLAR_DIR:$PATH
CC=gcc-11
CXX=g++-11
FC=gfortran-11
F77=${FC}
LD=${CXX}

# Libraries depemdencies
# arma
ARMA=/usr/local/Cellar/armadillo/10.6.1

# clags
ARMA_CFLAGS=-I/${ARMA}/include
GSL_CFLAGS=`gsl-config --cflags`
CFLAGS="${GSL_CFLAGS} ${ARMA_CFLAGS}"
CXXFLAGS="${CFLAGS} -I/usr/local/carma/include"

# library
LDLIBS=`gsl-config --libs`
LDFLAGS="$LDLIBS -L${ARMA}/lib/"
LIBS=
#------------------------------------------------------------------------------
# Todo:
#     check if LDLFAGS or LDLIBS are required 
export CC CXX FC CFLAGS CXXFLAGS LD LIBS LDFLAGS LDLIBS
cat <<EOF
CFLAGS:      $CFLAGS
CXXFLAGS:    $CXXFLAGS
FFLAGS:      $FFLAGS
LIBS:        $LIBS
LDFLAGS:     $LDFLAGS
LDLIBS:      $LDLIBS
EOF

# include a check that if TRACY LIB is set
if [ -z "$TRACY_LIB" ]
then
    echo '$TRACY_LIB has to be defined!'
    exit 1
fi

dir=`pwd`

cd "$TRACY_LIB"

# clean up build up
rm -rf autom4te.cache
rm -rf aclocal.m4
rm -rf tracy/lib/*

make distclean

./bootstrap
./configure --prefix=$TRACY_LIB/tracy
make install


