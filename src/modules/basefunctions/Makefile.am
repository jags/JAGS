dnl AM_CPPFLAGS = -I @top_srcdir@/include
INCLUDES = -I $(top_srcdir)/src/include

jagsmod_LTLIBRARIES = basefunctions.la

basefunctions_la_LDFLAGS = -no-undefined -module -avoid-version

basefunctions_la_LIBADD = $(top_builddir)/src/lib/libjags.la

basefunctions_la_SOURCES = Add.cc Divide.cc Equals.cc Multiply.cc Neg.cc Pow.cc \
Subtract.cc module.cc 

noinst_HEADERS = Add.h Divide.h Equals.h Multiply.h Neg.h Pow.h Subtract.h 
