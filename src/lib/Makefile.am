lib_LTLIBRARIES = libjags.la

libjags_la_SOURCES = Console.cc Module.cc

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

libjags_la_LDFLAGS = -no-undefined -version-info 1:1:0

SUBDIRS = util sarray distribution function sampler graph model	\
	compiler rng 
