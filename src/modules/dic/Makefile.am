jagsmod_LTLIBRARIES = dic.la

dic_la_CPPFLAGS = -I$(top_srcdir)/src/include

dic_la_LDFLAGS = -no-undefined -module -avoid-version

dic_la_LIBADD = $(top_builddir)/src/lib/libjags.la

dic_la_SOURCES = dic.cc DevianceMonitor.cc DevianceMonitorFactory.cc	\
PDMonitor.cc PDMonitorFactory.cc PoptMonitor.cc PoptMonitorFactory.cc

noinst_HEADERS = DevianceMonitor.h DevianceMonitorFactory.h PDMonitor.h PDMonitorFactory.h PoptMonitor.h PoptMonitorFactory.h
