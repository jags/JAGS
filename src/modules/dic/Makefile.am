jagsmod_LTLIBRARIES = dic.la

dic_la_CPPFLAGS = -I$(top_srcdir)/src/include

dic_la_LIBADD = $(top_builddir)/src/lib/libjags.la

dic_la_SOURCES = dic.cc DevianceMean.cc DevianceTrace.cc		\
DevianceMonitorFactory.cc PDMonitor.cc PoptMonitor.cc			\
PDMonitorFactory.cc KL.cc KLNorm.cc KLBin.cc KLPois.cc KLTab.cc		\
CalKL.cc CalKLExact.cc CalKLApprox.cc PDTrace.cc PDTraceFactory.cc

noinst_HEADERS = DevianceMean.h DevianceTrace.h				\
DevianceMonitorFactory.h PDMonitor.h PoptMonitor.h PDMonitorFactory.h	\
KL.h KLNorm.h KLBin.h KLPois.h KLTab.h CalKL.h CalKLExact.h		\
CalKLApprox.h PDTrace.h PDTraceFactory.h

dic_la_LDFLAGS = -module -avoid-version
if WINDOWS
dic_la_LDFLAGS += -no-undefined
endif
