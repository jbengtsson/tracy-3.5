# Process this file with autoconf to produce a configure script.
AC_INIT([Tracy],[2.0])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/main.cc)
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS(config.h)
# Checks for programs.
AC_PROG_AWK
# Use GNU vs. Clang/LLVM compilers.
#AC_PROG_CC(gcc-14)
AC_PROG_CC(cc)
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CPP
# Obsolete.
# AC_PROG_RANLIB
LT_INIT
# Use GNU vs. Clang/LLVM compilers.
#AC_PROG_CXX(g++-14)
AC_PROG_CXX(c++)
AC_PROG_F77([gfortran gfortran-15 gfortran-14 f77])
# Sets FLIBS: the FORTRAN runtime libraries. Needed because these programs
# link Fortran objects but are linked by the C++ driver, which does not add
# them itself.
AC_F77_LIBRARY_LDFLAGS
#AC_PROG_MAKE_SET(make)

# For macOS.
AC_CANONICAL_HOST
case "$host_os" in
  darwin*)
    MAC_LDFLAGS="-F/System/Library/Frameworks -framework CoreServices"
    ;;
  *)
    MAC_LDFLAGS=""
    ;;
esac
AC_SUBST(MAC_LDFLAGS)


# Checks for libraries.

# GSL. gsl-config ships with GSL and reports its own include and library
# paths, which differ across Linux, Intel macOS and Apple Silicon.
AC_ARG_WITH([gsl],
  [AS_HELP_STRING([--with-gsl=DIR],
                  [prefix of a GSL install, if gsl-config is not on PATH])])
AS_IF([test -n "$with_gsl"],
      [GSL_CONFIG="$with_gsl/bin/gsl-config"
       AS_IF([test ! -x "$GSL_CONFIG"],
             [AC_MSG_ERROR([no gsl-config at $GSL_CONFIG -- check --with-gsl])])],
      [AC_PATH_PROG([GSL_CONFIG], [gsl-config], [no])])
AS_IF([test "x$GSL_CONFIG" = xno],
      [AC_MSG_ERROR([GSL not found. Install it (brew install gsl, or
       apt install libgsl-dev) or pass --with-gsl=DIR.])])
GSL_CFLAGS=`$GSL_CONFIG --cflags`
GSL_LIBS=`$GSL_CONFIG --libs`
GSL_VERSION=`$GSL_CONFIG --version`
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

# The scripts link against three trees located by environment variable. When
# unset these expand to nonsense paths such as /thor/lib/libTPSALib.a and the
# failure only surfaces at link time, so check them here instead.
AS_IF([test ! -d "$TRACY_LIB/tracy/lib"],
      [AC_MSG_ERROR([TRACY_LIB is unset or wrong: '$TRACY_LIB'.
       Set it to the tracy-3.5 tree root and build the library there first
       (./make_tracy-3.5.sh). See README.rst.])])
AS_IF([test ! -d "$NUM_REC/lib"],
      [AC_MSG_ERROR([NUM_REC is unset or wrong: '$NUM_REC'.
       Set it to your Numerical Recipes tree. See README.rst.])])
AS_IF([test ! -d "$THOR_LIB/thor/lib"],
      [AC_MSG_ERROR([THOR_LIB is unset or wrong: '$THOR_LIB'.
       Set it to your thor tree; it is needed by projects/src/ptc.
       See README.rst.])])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CONFIG_FILES([Makefile src/Makefile src/ptc/Makefile])
AC_OUTPUT

# Report what was actually found, so build reports from different machines
# can be compared directly. Versions are resolved first: command substitution
# is not performed inside the AC_MSG_NOTICE argument.
tracy_cc_version=`$CC --version 2>/dev/null | head -n 1`
tracy_cxx_version=`$CXX --version 2>/dev/null | head -n 1`
tracy_f77_version=`$F77 --version 2>/dev/null | head -n 1`
AC_MSG_NOTICE([

  Tracy $PACKAGE_VERSION -- scripts build configured

    host:            $host ($host_cpu / $host_os)

    C compiler:      $CC
      version:       $tracy_cc_version
    C++ compiler:    $CXX
      version:       $tracy_cxx_version
    FORTRAN 77:      $F77
      version:       $tracy_f77_version

    GSL:             $GSL_VERSION
      gsl-config:    $GSL_CONFIG
      cflags:        $GSL_CFLAGS
      libs:          $GSL_LIBS

    TRACY_LIB:       $TRACY_LIB
    NUM_REC:         $NUM_REC
    THOR_LIB:        $THOR_LIB

    extra LDFLAGS:   ${MAC_LDFLAGS:-(none)}

  Next: make -- the scripts land in src/ and src/ptc/.
])
