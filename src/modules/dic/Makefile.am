jagsmod_LTLIBRARIES = dic.la

dic_la_CPPFLAGS = -I$(top_srcdir)/src/include

dic_la_LDFLAGS = -no-undefined -module -avoid-version

dic_la_LIBADD = $(top_builddir)/src/lib/libjags.la

dic_la_SOURCES = dic.cc DevianceMonitor.cc DevianceMonitorFactory.cc	\
PDMonitor.cc PDMonitorFactory.cc PoptMonitor.cc PoptMonitorFactory.cc	\
DefaultPDMonitor.cc KLPDMonitor.cc KLNorm.cc KLBin.cc KLPois.cc KLTab.cc

noinst_HEADERS = DevianceMonitor.h DevianceMonitorFactory.h		\
PDMonitor.h PDMonitorFactory.h PoptMonitor.h PoptMonitorFactory.h	\
DefaultPDMonitor.h KLPDMonitor.h KLNorm.h KLBin.h KLPois.h KLTab.h KL.h
