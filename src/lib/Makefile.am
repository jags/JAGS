SUBDIRS = util sarray distribution function sampler graph model	\
	compiler rng module

extra_dist = version.cc.in

lib_LTLIBRARIES = libjags.la

libjags_la_SOURCES = version.cc Console.cc

libjags_la_CPPFLAGS = -I$(top_srcdir)/src/include -I$(top_builddir)/src/include

libjags_la_LIBADD = model/libmodel.la		\
        util/libutil.la				\
        sarray/libsarray.la			\
        function/libfunc.la			\
        sampler/libsampler.la			\
        distribution/libdist.la			\
	graph/libgraph.la			\
	compiler/libcompiler.la			\
        rng/librng.la				\
	module/libmodule.la

libjags_la_LDFLAGS = -version-info $(JAGS_MAJOR):$(JAGS_MINOR):0
if WINDOWS
libjags_la_LDFLAGS += -no-undefined
libjags_la_SOURCES += dllversion.rc
endif

### Test library 

if CANCHECK
check_LTLIBRARIES = libtest.la
libtest_la_SOURCES = testlib.cc
libtest_la_CPPFLAGS = -I$(top_srcdir)/src/include
libtest_la_LIBADD = function/libfuntest.la

if WINDOWS
libtest_la_LDFLAGS = -no-undefined
endif

endif
