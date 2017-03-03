SUBDIRS = distributions samplers

jagsmod_LTLIBRARIES = mix.la

mix_la_SOURCES = mix.cc

mix_la_CPPFLAGS = -I$(top_srcdir)/src/include

mix_la_LIBADD = distributions/libmixdist.la               \
        samplers/libmixsamp.la         			\
        $(top_builddir)/src/jrmath/libjrmath.la		\
        $(top_builddir)/src/lib/libjags.la 

mix_la_LDFLAGS = -module -avoid-version
if WINDOWS
mix_la_LDFLAGS += -no-undefined
endif

### Test library 

if CANCHECK
check_LTLIBRARIES = libmixtest.la
libmixtest_la_SOURCES = testmix.cc testmix.h
libmixtest_la_CPPFLAGS = -I$(top_srcdir)/src/include
libmixtest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
libmixtest_la_LDFLAGS = $(CPPUNIT_LDFLAGS)
libmixtest_la_LIBADD = 	distributions/libmixdisttest.la		\
	distributions/libmixdist.la				\
	$(top_builddir)/src/modules/base/rngs/libbaserngs.la	\
	$(top_builddir)/src/lib/libtest.la			\
	$(top_builddir)/src/lib/libjags.la 			\
	$(top_builddir)/src/jrmath/libjrmath.la 		\
	@LAPACK_LIBS@ @BLAS_LIBS@

if WINDOWS
libmixtest_la_LDFLAGS += -no-undefined
else
libmixtest_la_LIBADD += @FLIBS@
endif

endif
