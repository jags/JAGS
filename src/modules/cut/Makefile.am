jagsmod_LTLIBRARIES = cut.la

cut_la_CPPFLAGS = -I$(top_srcdir)/src/include

cut_la_LIBADD = $(top_builddir)/src/lib/libjags.la 		\

cut_la_SOURCES = cut.cc DCut.cc

noinst_HEADERS = DCut.h

cut_la_LDFLAGS = -module -avoid-version
if WINDOWS
cut_la_LDFLAGS += -no-undefined
endif
