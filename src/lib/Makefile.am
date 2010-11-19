extra_dist = version.cc.in

lib_LTLIBRARIES = libjags.la

libjags_la_SOURCES = version.cc Console.cc Module.cc

libjags_la_CPPFLAGS = -I$(top_srcdir)/src/include

libjags_la_LIBADD = model/libmodel.la		\
        util/libutil.la				\
        sarray/libsarray.la			\
        function/libfunc.la			\
        sampler/libsampler.la			\
        distribution/libdist.la			\
	graph/libgraph.la			\
	compiler/libcompiler.la			\
        rng/librng.la

libjags_la_LDFLAGS = -version-info 2:2:0
if WINDOWS
libjags_la_LDFLAGS += -no-undefined
endif

SUBDIRS = util sarray distribution function sampler graph model	\
	compiler rng 
