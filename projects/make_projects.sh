#!/bin/sh
# Build the Tracy scripts in src/ and src/ptc/.
#
# Requires TRACY_LIB, NUM_REC and THOR_LIB to be set, and the Tracy libraries
# to have been built first (../make_tracy-3.5.sh). See ../README.rst.
#
# Any arguments are passed through to ./configure.

dir=`pwd`

rm -rf autom4te.cache
rm -rf aclocal.m4

# Create config dir – if it doesn't exist.
mkdir -p config

# Only meaningful once a Makefile exists, and a failure is recoverable: the
# bootstrap below regenerates everything anyway.
if [ -f Makefile ]; then
    make distclean || echo "make distclean failed; continuing with a clean bootstrap"
fi

# Past this point every step is required, so stop on the first failure.
set -e

# Configure libtool. This generates the m4/ macros that aclocal needs. m4/ is
# not tracked in git, so it must be regenerated in every fresh checkout.
case "$(uname -s)" in
    Linux*)
        echo "Running on Linux"
        libtoolize
        ;;
    Darwin*)
        echo "Running on macOS"
        glibtoolize
        ;;
    *)
        echo "Unrecognised platform $(uname -s); trying libtoolize"
        libtoolize
        ;;
esac

./bootstrap
./configure --prefix="$dir/projects" "$@"

make
