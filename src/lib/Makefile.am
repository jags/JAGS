INCLUDES = -I $(top_srcdir)/src/include

jagslib_LTLIBRARIES = libjags.la

libjags_la_SOURCES = Console.cc Module.cc

libjags_la_LIBADD = model/libmodel.la		\
        sarray/libsarray.la			\
        function/libfunc.la			\
        sampler/libsampler.la			\
        distribution/libdist.la			\
	graph/libgraph.la			\
	compiler/libcompiler.la			\
        rng/librng.la

libjags_la_LDFLAGS = -no-undefined

SUBDIRS = sarray distribution function sampler graph model	\
	compiler rng
