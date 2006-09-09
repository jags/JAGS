INCLUDES = -I $(top_srcdir)/src/include -I$(top_srcdir)/src/modules/bugs/matrix \
           -I$(top_srcdir)/src/modules/bugs/rmath

noinst_LTLIBRARIES = bugsdist.la

bugsdist_la_LDFLAGS = -no-undefined -module -avoid-version

bugsdist_la_SOURCES = DBern.cc DCat.cc DDirch.cc DHyper.cc DLogis.cc	\
DMulti.cc DNormMix.cc DSum.cc DWeib.cc DBeta.cc DChisqr.cc DExp.cc	\
DInterval.cc DMNorm.cc DNegBin.cc DPar.cc DT.cc DWish.cc DBin.cc	\
DDexp.cc DGamma.cc DLnorm.cc DMState.cc DNorm.cc DPois.cc DUnif.cc 

noinst_HEADERS = DBern.h DCat.h DDirch.h DHyper.h DLogis.h DMulti.h	\
DNormMix.h DSum.h DWeib.h DBeta.h DChisqr.h DExp.h DInterval.h		\
DMNorm.h DNegBin.h DPar.h DT.h DWish.h DBin.h DDexp.h DGamma.h		\
DLnorm.h DMState.h DNorm.h DPois.h DUnif.h 

