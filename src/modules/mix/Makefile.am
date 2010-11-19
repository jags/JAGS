SUBDIRS = distributions samplers

jagsmod_LTLIBRARIES = mix.la

mix_la_SOURCES = mix.cc

mix_la_CPPFLAGS = -I$(top_srcdir)/src/include

mix_la_LIBADD = distributions/mixdist.la               \
        samplers/mixsamp.la         			\
        $(top_builddir)/src/jrmath/libjrmath.la		\
        $(top_builddir)/src/lib/libjags.la 

mix_la_LDFLAGS = -module -avoid-version
if WINDOWS
mix_la_LDFLAGS += -no-undefined
endif