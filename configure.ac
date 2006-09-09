dnl Process this file with autoconf to produce a configure script.
AC_INIT([JAGS], [0.97.1], [jags@iarc.fr], [JAGS])
AC_CONFIG_SRCDIR(src/lib/Console.cc)
AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)

AC_PROG_CXX
AM_PROG_LEX
AC_PROG_YACC

case "${host_os}" in
  darwin*)
    ## Allow linking of code compiled with Apple GCC and code
    ## compiled with GNU GCC.
    AC_CHECK_LIB(cc_dynamic, main)
    ;;
  *)
    ;;
esac

dnl libtool dynamic loading
AC_LIBTOOL_DLOPEN
AC_LIBLTDL_CONVENIENCE
AC_PROG_LIBTOOL
AC_CONFIG_SUBDIRS(libltdl)
AC_SUBST(INCLTDL)
AC_SUBST(LIBLTDL)

dnl fortran stuff
ACX_F77_WRAPPERS
ACX_F77_LIBRARY_LDFLAGS
ACX_LAPACK([],[
	AC_MSG_ERROR("You need to install the LAPACK library")
])

dnl Rmath stuff

AC_DEFINE(MATHLIB_STANDALONE, 1, [Define if you have standalone R math library])

dnl Relevant parts lifted from R's configure.ac, version 2.3.0
AC_CHECK_LIB(m, sin)

## IEEE 754. We rely on this in e.g. the working log test.
R_IEEE_754

R_CHECK_FUNCS([expm1 log1p], [#include <math.h>])
## Used to build src/include/Rmath.h.
## <NOTE>
## we don't use AC_CONFIG_HEADERS on Rmath.h.in because
## a) that would comment out #undef statements in Rmath.h.in and
## b) Rmath.h should be a self-contained file for standalone Rmath use.
## </NOTE>
#if test "${ac_cv_func_expm1}" = yes; then
if test "${ac_cv_have_decl_expm1}" = yes; then
  RMATH_HAVE_EXPM1="# define HAVE_EXPM1 1"
else
  RMATH_HAVE_EXPM1="# undef HAVE_EXPM1"
fi
AC_SUBST(RMATH_HAVE_EXPM1)
#if test "${ac_cv_func_log1p}" = yes; then
if test "${ac_cv_have_decl_expm1}" = yes; then
  RMATH_HAVE_LOG1P="# define HAVE_LOG1P 1"
else
  RMATH_HAVE_LOG1P="# undef HAVE_LOG1P"
fi
AC_SUBST(RMATH_HAVE_LOG1P)

R_FUNC_LOG
R_FUNC_LOG1P

dnl End Rmath stuff

jagshomedir=${libdir}/jags
AC_SUBST(jagshomedir)
jagslibdir=${jagshomedir}/lib
AC_SUBST(jagslibdir)
jagsincludedir=${jagshomedir}/include
AC_SUBST(jagsincludedir)
jagsmoddir=${jagshomedir}/modules
AC_SUBST(jagsmoddir)

AC_CONFIG_HEADERS([src/include/config.h])
AC_CONFIG_FILES([
  Makefile
  m4/Makefile
  etc/Makefile
  etc/Makeconf
  src/Makefile
  src/lib/Makefile
  src/lib/sarray/Makefile
  src/lib/function/Makefile
  src/lib/distribution/Makefile
  src/lib/graph/Makefile
  src/lib/sampler/Makefile
  src/lib/model/Makefile
  src/lib/compiler/Makefile
  src/lib/rng/Makefile
  src/include/Makefile
  src/include/sarray/Makefile
  src/include/function/Makefile
  src/include/distribution/Makefile
  src/include/graph/Makefile
  src/include/sampler/Makefile
  src/include/model/Makefile
  src/include/compiler/Makefile
  src/include/rng/Makefile
  src/modules/Makefile
  src/modules/basefunctions/Makefile
  src/modules/baserngs/Makefile
  src/modules/basesamplers/Makefile
  src/modules/bugs/Makefile	
  src/modules/bugs/functions/Makefile	
  src/modules/bugs/distributions/Makefile	
  src/modules/bugs/samplers/Makefile	
  src/modules/bugs/rmath/Makefile	
  src/modules/bugs/matrix/Makefile	
  src/terminal/Makefile
])
AC_OUTPUT
