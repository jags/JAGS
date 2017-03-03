noinst_LTLIBRARIES = libmixdist.la

libmixdist_la_CPPFLAGS = -I$(top_srcdir)/src/include

libmixdist_la_LDFLAGS = -no-undefined -module -avoid-version

libmixdist_la_SOURCES = DNormMix.cc DBetaBin.cc DPick.cc

noinst_HEADERS = DNormMix.h DBetaBin.h DPick.h

### Test library 

if CANCHECK
check_LTLIBRARIES = libmixdisttest.la
libmixdisttest_la_SOURCES = testmixdist.cc testmixdist.h
libmixdisttest_la_CPPFLAGS = -I$(top_srcdir)/src/include \
	-I$(top_srcdir)/src/modules/base/rngs
libmixdisttest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
endif
