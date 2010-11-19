jagsmod_LTLIBRARIES = lecuyer.la

lecuyer_la_SOURCES = lecuyer.cc RngStream.cc RngStreamFactory.cc

lecuyer_la_CPPFLAGS = -I$(top_srcdir)/src/include

lecuyer_la_LIBADD = $(top_builddir)/src/lib/libjags.la

noinst_HEADERS = RngStream.h RngStreamFactory.h

lecuyer_la_LDFLAGS = -module -avoid-version
if WINDOWS
lecuyer_la_LDFLAGS += -no-undefined
endif

