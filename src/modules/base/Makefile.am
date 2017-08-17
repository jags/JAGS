SUBDIRS = functions samplers monitors rngs

### Base module: note that the name is "basemod" not "base" to avoid
### conflict with the base package in R.

jagsmod_LTLIBRARIES = basemod.la

basemod_la_SOURCES = base.cc

basemod_la_CPPFLAGS = -I$(top_srcdir)/src/include

basemod_la_LIBADD = functions/libbasefunctions.la  \
        monitors/libbasemonitors.la                \
        samplers/libbasesamplers.la                \
        rngs/libbaserngs.la                        \
	$(top_builddir)/src/jrmath/libjrmath.la \
	$(top_builddir)/src/lib/libjags.la        

basemod_la_LDFLAGS = -module -avoid-version
if WINDOWS
basemod_la_LDFLAGS += -no-undefined
endif

### Test library 

if CANCHECK
check_LTLIBRARIES = libbasetest.la
libbasetest_la_SOURCES = testbase.cc testbase.h
libbasetest_la_CPPFLAGS = -I$(top_srcdir)/src/include
libbasetest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
libbasetest_la_LDFLAGS = $(CPPUNIT_LIBS)
libbasetest_la_LIBADD = functions/libbasefuntest.la	\
	functions/libbasefunctions.la			\
	$(top_builddir)/src/lib/libtest.la		\
	$(top_builddir)/src/lib/libjags.la
	$(top_builddir)/src/jrmath/libjrmath.la
endif

