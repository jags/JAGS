noinst_LTLIBRARIES = libglmdist.la

libglmdist_la_CPPFLAGS = -I$(top_srcdir)/src/include

libglmdist_la_SOURCES = DScaledGamma.cc DScaledWishart.cc

noinst_HEADERS = DScaledGamma.h DScaledWishart.h

## Test library 

if CANCHECK
check_LTLIBRARIES = libglmdisttest.la
libglmdisttest_la_SOURCES = testglmdist.cc testglmdist.h
libglmdisttest_la_CPPFLAGS = -I$(top_srcdir)/src/include \
	-I$(top_srcdir)/src/modules/base/rngs
libglmdisttest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
endif
