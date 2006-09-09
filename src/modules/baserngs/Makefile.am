dnl AM_CPPFLAGS = -I @top_srcdir@/include
INCLUDES = -I $(top_srcdir)/src/include

jagsmod_LTLIBRARIES = baserngs.la

baserngs_la_LDFLAGS = -no-undefined -module -avoid-version

baserngs_la_LIBADD = $(top_builddir)/src/lib/libjags.la

baserngs_la_SOURCES = MarsagliaRNG.cc WichmannHillRNG.cc SuperDuperRNG.cc \
MersenneTwisterRNG.cc BaseRNGFactory.cc module.cc

noinst_HEADERS = MarsagliaRNG.h WichmannHillRNG.h SuperDuperRNG.h \
MersenneTwisterRNG.h BaseRNGFactory.h
