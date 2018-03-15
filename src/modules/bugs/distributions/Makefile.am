noinst_LTLIBRARIES = libbugsdist.la

libbugsdist_la_CPPFLAGS = -I$(top_srcdir)/src/include	\
-I$(top_srcdir)/src/modules/bugs/matrix

libbugsdist_la_SOURCES = DBern.cc DCat.cc DDirch.cc DHyper.cc DLogis.cc	\
DMulti.cc DSum.cc DWeib.cc DBeta.cc DChisqr.cc DExp.cc DInterval.cc	\
DMNorm.cc DNegBin.cc DPar.cc DT.cc DWish.cc DBin.cc DDexp.cc DGamma.cc	\
DLnorm.cc DNorm.cc DPois.cc DUnif.cc DMT.cc DGenGamma.cc		\
DF.cc DNChisqr.cc DRound.cc DNT.cc SumDist.cc DSample.cc DRW1.cc 	\
DMNormVC.cc DGamPois.cc

noinst_HEADERS = DBern.h DCat.h DDirch.h DHyper.h DLogis.h DMulti.h	\
DSum.h DWeib.h DBeta.h DChisqr.h DExp.h DInterval.h DMNorm.h		\
DNegBin.h DPar.h DT.h DWish.h DBin.h DDexp.h DGamma.h DLnorm.h		\
DNorm.h DPois.h DUnif.h DMT.h DGenGamma.h DF.h DNChisqr.h DRound.h	\
DNT.h SumDist.h DSample.h DRW1.h DMNormVC.h DGamPois.h

### Test library 

if CANCHECK
check_LTLIBRARIES = libbugsdisttest.la
libbugsdisttest_la_SOURCES = testbugsdist.cc testbugsdist.h
libbugsdisttest_la_CPPFLAGS = -I$(top_srcdir)/src/include \
	-I$(top_srcdir)/src/modules/base/rngs
libbugsdisttest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
endif
