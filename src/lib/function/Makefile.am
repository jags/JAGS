noinst_LTLIBRARIES = libfunc.la

libfunc_la_CPPFLAGS = -I$(top_srcdir)/src/include

libfunc_la_SOURCES = Function.cc ScalarFunction.cc LinkFunction.cc	\
VectorFunction.cc ArrayFunction.cc FunctionPtr.cc FuncTab.cc FuncError.cc \
DPQFunction.cc DFunction.cc PFunction.cc QFunction.cc \
ScalarLogDensity.cc VectorLogDensity.cc ArrayLogDensity.cc \
ScalarVectorFunction.cc

if CANCHECK
check_LTLIBRARIES = libfuntest.la
libfuntest_la_SOURCES = testfun.cc
libfuntest_la_CPPFLAGS = -I$(top_srcdir)/src/include
libfuntest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
libfuntest_la_LDFLAGS = $(CPPUNIT_LIBS) 
endif
