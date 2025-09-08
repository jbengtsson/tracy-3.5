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
AC_PROG_F77(gfortran-14)
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

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CONFIG_FILES([Makefile src/Makefile src/ptc/Makefile])
AC_OUTPUT
