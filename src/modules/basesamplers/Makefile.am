AM_CPPFLAGS = -I @top_srcdir@/src/include

jagsmod_LTLIBRARIES = basesamplers.la

basesamplers_la_LDFLAGS = -no-undefined -module -avoid-version

basesamplers_la_LIBADD = $(top_builddir)/src/lib/libjags.la

basesamplers_la_SOURCES = DiscreteSliceSampler.cc FiniteFactory.cc	\
FiniteSampler.cc RealSliceSampler.cc SliceFactory.cc module.cc

noinst_HEADERS = DiscreteSliceSampler.h FiniteFactory.h		\
FiniteSampler.h RealSliceSampler.h SliceFactory.h 

