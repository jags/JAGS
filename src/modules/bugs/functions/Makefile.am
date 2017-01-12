noinst_LTLIBRARIES = libbugsfunc.la

libbugsfunc_la_CPPFLAGS = -I$(top_srcdir)/src/include	\
-I$(top_srcdir)/src/modules/bugs/matrix \
-I$(top_srcdir)/src/modules/bugs/distributions

libbugsfunc_la_LDFLAGS = -no-undefined -module -avoid-version

libbugsfunc_la_SOURCES = Max.cc Probit.cc InProd.cc Mean.cc SD.cc	\
ICLogLog.cc Min.cc Sqrt.cc CLogLog.cc ILogit.cc Step.cc Log.cc Sum.cc	\
Exp.cc LogFact.cc LogGam.cc Phi.cc Logit.cc Inverse.cc LogDet.cc	\
Sort.cc Rank.cc Abs.cc Cos.cc Sin.cc Tan.cc Cosh.cc Sinh.cc Tanh.cc	\
ArcCos.cc ArcSin.cc ArcTan.cc ArcCosh.cc ArcSinh.cc ArcTanh.cc		\
InterpLin.cc Round.cc Trunc.cc Equals.cc Prod.cc MatMult.cc		\
Transpose.cc DSumFunc.cc DIntervalFunc.cc DRoundFunc.cc IfElse.cc	\
Order.cc Combine.cc Rep.cc

noinst_HEADERS = LogGam.h Phi.h Sum.h Log.h Logit.h Probit.h		\
CLogLog.h InProd.h Max.h SD.h ICLogLog.h Mean.h Sqrt.h ILogit.h Min.h	\
Step.h Exp.h LogFact.h Inverse.h LogDet.h Sort.h Rank.h Abs.h Cos.h	\
Sin.h Tan.h Cosh.h Sinh.h Tanh.h ArcCos.h ArcSin.h ArcTan.h ArcCosh.h	\
ArcSinh.h ArcTanh.h InterpLin.h Round.h Trunc.h Equals.h Prod.h		\
MatMult.h Transpose.h DSumFunc.h DIntervalFunc.h DRoundFunc.h		\
IfElse.h Order.h Combine.h Rep.h

### Test library 

if CANCHECK
check_LTLIBRARIES = libbugsfuntest.la
libbugsfuntest_la_SOURCES = testbugsfun.cc testbugsfun.h
libbugsfuntest_la_CPPFLAGS = -I$(top_srcdir)/src/include
libbugsfuntest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
endif
