# C macro processor.
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([tracy/src/t2elem.cc])
AC_INIT([Tracy],[2.0])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(tracy/src/t2elem.cc)


# search for GSL see https://www.gnu.org/software/gsl/doc/html/autoconf.html
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])


AM_INIT_AUTOMAKE
# Shared library: libtool.
LT_CONFIG_LTDL_DIR()
# Remark: Execute glibtoolize before aclocal & automake.
#LT_INIT([disable-static])
#LT_INIT([disable-shared])
LT_INIT
AC_CONFIG_HEADERS(config.h)
# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CPP
# For MacBook: use GNU gcc vs. LLVM.
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
AC_CONFIG_FILES([Makefile
          TPSA/Makefile
	  tracy/Makefile tracy/src/Makefile])
AC_OUTPUT
