noinst_LTLIBRARIES = libglmsampler.la

libglmsampler_la_CPPFLAGS = -I$(top_srcdir)/src/include \
		-I$(top_srcdir)/src/modules/glm/CSparse \
		-I$(top_srcdir)/src/modules/glm/SSparse/config \
		-I$(top_srcdir)/src/modules/glm/SSparse/CHOLMOD/Include

libglmsampler_la_SOURCES = GLMFactory.cc GLMSampler.cc GLMMethod.cc	\
 KS.cc			\
 IWLSFactory.cc	\
 IWLS.cc LGMix.cc AuxMixPoisson.cc AuxMixBinomial.cc Outcome.cc		\
 NormalLinear.cc BinaryProbit.cc BinaryLogit.cc Classify.cc		\
 IWLSOutcome.cc HolmesHeld.cc HolmesHeldFactory.cc GLMBlock.cc		\
 GLMGibbs.cc GLMGenericFactory.cc HolmesHeldGibbs.cc PolyaGamma.cc PGcommon.cc \
 ScaledWishart.cc ScaledWishartFactory.cc SampleWishart.cc \
 ScaledGamma.cc ScaledGammaFactory.cc \
 LogisticLinear.cc TLinear.cc OrderedLogit.cc OrderedProbit.cc \
 LogNormalLinear.cc \
 RESampler.cc REMethod.cc REFactory.cc \
 REGamma.cc REGammaSlicer.cc REGammaFactory.cc \
 REScaledGamma.cc REScaledGammaFactory.cc \
 REScaledWishart.cc REScaledWishartFactory.cc \
 MNormalLinear.cc \
 REMethod2.cc REFactory2.cc \
 REGammaSlicer2.cc REGamma2.cc REGammaFactory2.cc \
 REScaledGamma2.cc REScaledGammaFactory2.cc \
 REScaledWishart2.cc REScaledWishartFactory2.cc

noinst_HEADERS = GLMFactory.h GLMSampler.h GLMMethod.h			\
  KS.h 	\
  IWLSFactory.h IWLS.h LGMix.h		\
  AuxMixPoisson.h AuxMixBinomial.h Outcome.h	\
  NormalLinear.h BinaryProbit.h BinaryLogit.h Classify.h		\
  IWLSOutcome.h HolmesHeld.h HolmesHeldFactory.h GLMBlock.h		\
  GLMGibbs.h GLMGenericFactory.h HolmesHeldGibbs.h PolyaGamma.h PG.h    \
  ScaledWishart.h ScaledWishartFactory.h SampleWishart.h \
  ScaledGamma.h ScaledGammaFactory.h \
  LogisticLinear.h TLinear.h OrderedLogit.h OrderedProbit.h \
  LogNormalLinear.h \
  RESampler.h REMethod.h REFactory.h \
  REGamma.h REGammaSlicer.h REGammaFactory.h \
  REScaledGamma.h REScaledGammaFactory.h \
  REScaledWishart.h REScaledWishartFactory.h \
  MNormalLinear.h \
  REMethod2.h REFactory2.h \
  REGammaSlicer2.h REGamma2.h REGammaFactory2.h \
  REScaledGamma2.h REScaledGammaFactory2.h \
  REScaledWishart2.h REScaledWishartFactory2.h

if CANCHECK
check_LTLIBRARIES = libglmsamptest.la
libglmsamptest_la_SOURCES = testglmsamp.cc testglmsamp.h
libglmsamptest_la_CPPFLAGS = -I$(top_srcdir)/src/include
libglmsamptest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
endif
