INCLUDES = -I $(top_srcdir)/src/include \
-I$(top_srcdir)/src/modules/bugs/distributions \
-I$(top_srcdir)/src/modules/bugs/matrix \
-I$(top_srcdir)/src/modules/bugs/rmath 

noinst_LTLIBRARIES = bugssampler.la

bugssampler_la_LDFLAGS = -no-undefined -module -avoid-version

bugssampler_la_SOURCES = Censored.cc ConjugateGamma.cc		\
ConjugateWishart.cc ConjugateBeta.cc		\
ConjugateMNormal.cc 	\
ConjugateDirichlet.cc ConjugateNormal.cc DSumFactory.cc		\
ConjugateFactory.cc ConjugateSampler.cc	\
DSumSampler.cc  

noinst_HEADERS = Censored.h ConjugateFactory.h ConjugateNormal.h	\
	\
ConjugateBeta.h ConjugateGamma.h ConjugateSampler.h DSumFactory.h	\
ConjugateDirichlet.h ConjugateMNormal.h	\
ConjugateWishart.h DSumSampler.h

