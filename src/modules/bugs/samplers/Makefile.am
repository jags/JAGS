# SUBDIRS = cj

noinst_LTLIBRARIES = libbugssampler.la

libbugssampler_la_CPPFLAGS = -I$(top_srcdir)/src/include	\
-I$(top_srcdir)/src/modules/bugs/distributions		\
-I$(top_srcdir)/src/modules/bugs/matrix

libbugssampler_la_SOURCES = Censored.cc CensoredFactory.cc		\
ConjugateGamma.cc ConjugateWishart.cc ConjugateBeta.cc			\
ConjugateMNormal.cc ConjugateDirichlet.cc ConjugateNormal.cc		\
DSumFactory.cc ConjugateFactory.cc RWDSum.cc RealDSum.cc		\
DiscreteDSum.cc MNormal.cc MNormalFactory.cc ConjugateMethod.cc		\
Dirichlet.cc DirichletFactory.cc TruncatedGamma.cc DMultiDSum.cc	\
ShiftedCount.cc ShiftedMultinomial.cc SumMethod.cc SumFactory.cc 	\
RW1.cc RW1Factory.cc BinomSlicer.cc BinomSliceFactory.cc

noinst_HEADERS = Censored.h CensoredFactory.h ConjugateFactory.h	\
ConjugateNormal.h ConjugateBeta.h ConjugateGamma.h DSumFactory.h	\
ConjugateDirichlet.h ConjugateMNormal.h ConjugateWishart.h RWDSum.h	\
RealDSum.h DiscreteDSum.h MNormal.h MNormalFactory.h			\
ConjugateMethod.h Dirichlet.h DirichletFactory.h TruncatedGamma.h	\
DMultiDSum.h ShiftedCount.h ShiftedMultinomial.h SumMethod.h		\
SumFactory.h RW1.h RW1Factory.h BinomSlicer.h BinomSliceFactory.h

