jagsmod_LTLIBRARIES = dic.la

dic_la_CPPFLAGS = -I$(top_srcdir)/src/include

dic_la_LDFLAGS = -module -avoid-version
if WINDOWS
dic_la_LDFLAGS += -no-undefined
endif

dic_la_LIBADD = $(top_builddir)/src/lib/libjags.la

dic_la_SOURCES = dic.cc DevianceMean.cc DevianceTrace.cc		\
DevianceMonitorFactory.cc PDMonitor.cc PoptMonitor.cc			\
PDMonitorFactory.cc PDTrace.cc PDTraceFactory.cc			\
WAICMonitorFactory.cc WAICMonitor.cc

noinst_HEADERS = DevianceMean.h DevianceTrace.h				\
DevianceMonitorFactory.h PDMonitor.h PoptMonitor.h PDMonitorFactory.h	\
PDTrace.h PDTraceFactory.h						\
WAICMonitorFactory.h WAICMonitor.h
