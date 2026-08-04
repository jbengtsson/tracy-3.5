# Process this file with autoconf to produce a configure script.
AC_CONFIG_SRCDIR([tracy/src/t2elem.cc])
AC_INIT([Tracy],[2.0])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(tracy/src/t2elem.cc)
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS([config.h])
# Checks for programs.
AC_PROG_AWK
# Use GNU vs. Clang/LLVM compilers.
#AC_PROG_CC(gcc-14)
AC_PROG_CC(cc)
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CPP
LT_INIT
# Use GNU vs. Clang/LLVM compilers.
#AC_PROG_CXX(g++-14)
AC_PROG_CXX(c++)
AC_PROG_F77([gfortran gfortran-15 gfortran-14 f77])
#AC_PROG_MAKE_SET(gmake)

AC_CANONICAL_HOST

# -mcmodel=medium is an x86-64 option; it is rejected outright on arm64.
case "$host_cpu" in
  x86_64|amd64)
    MCMODEL_FFLAGS="-mcmodel=medium"
    ;;
  *)
    MCMODEL_FFLAGS=""
    ;;
esac
AC_SUBST(MCMODEL_FFLAGS)

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

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CONFIG_FILES([Makefile
          TPSA/Makefile
	  tracy/Makefile tracy/src/Makefile])
AC_OUTPUT

# Report what was actually found, so build reports from different machines
# can be compared directly. Versions are resolved first: command substitution
# is not performed inside the AC_MSG_NOTICE argument.
tracy_cc_version=`$CC --version 2>/dev/null | head -n 1`
tracy_cxx_version=`$CXX --version 2>/dev/null | head -n 1`
tracy_f77_version=`$F77 --version 2>/dev/null | head -n 1`
AC_MSG_NOTICE([

  Tracy $PACKAGE_VERSION -- library build configured

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

    FFLAGS extra:    ${MCMODEL_FFLAGS:-(none -- not x86-64)}
    install prefix:  $prefix

  Next: make install, then build the scripts in projects/ -- see README.rst.
])
