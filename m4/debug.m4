AC_DEFUN([JAGS_DEBUG],
[
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_ENABLE([debug],
    AS_HELP_STRING([--enable-debug], [Enable additional debugging code]))

if test "x$enable_debug" = "xyes" ; then
  AC_DEFINE([HAVE_GCD], [1],
            [Define if you want to use Grand Central Dispatch for thread management])
fi
])# JAGS_DEBUG

