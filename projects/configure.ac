# Process this file with autoconf to produce a configure script.
AC_INIT(Tracy, 2.0)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/main.cc)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)
# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AM_PROG_CC_STDC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CPP
AC_PROG_RANLIB
# Use GNU gcc; vs. LLVM.
AC_PROG_CXX(gcc-11)
AC_PROG_F77(gfortran)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_OUTPUT(Makefile src/Makefile src/ptc/Makefile)
