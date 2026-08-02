AC_DEFUN([JAGS_GCD],
[
AC_REQUIRE([AC_CANONICAL_HOST])

use_gcd="no"
have_dispatch_h="no"
use_gcd_libs="no"

AC_ARG_ENABLE([gcd],
    AS_HELP_STRING([--enable-gcd], [Enable Grand Central Dispatch]))

if test "x$enable_gcd" = "xyes" ; then
  AC_MSG_NOTICE([GCD support requested via --enable-gcd])
  case "${host_os}" in
  darwin*)
    AC_CHECK_FUNC([dispatch_get_global_queue], [use_gcd=yes])
    ;;
  *)
    AC_CHECK_HEADER([dispatch/dispatch.h], [have_dispatch_h="yes"])
    if test "x${have_dispatch_h}" = "xyes" ; then
      AC_CHECK_LIB([dispatch], [dispatch_get_global_queue], [use_gcd_libs="yes"])
    fi
    if test "x$use_gcd_libs" = "xyes"; then
       use_gcd="yes"
    fi
    ;;
  esac
  
  if test "x${use_gcd}" = "xyes"; then
    AC_DEFINE([HAVE_GCD], [1],
              [Define if you want to use Grand Central Dispatch for thread management])
    if test "x$use_gcd_libs" = "xyes"; then
      GCD_LIBS="-ldispatch"
      AC_SUBST([GCD_LIBS])
    fi
  else
    AC_MSG_ERROR([GCD is not available])
  fi

fi
])# JAGS_GCD

