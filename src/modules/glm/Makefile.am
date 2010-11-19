SUBDIRS = CSparse SSparse samplers 

jagsmod_LTLIBRARIES = glm.la

glm_la_SOURCES = glm.cc

glm_la_CPPFLAGS = -I$(top_srcdir)/src/include	\
-I$(top_srcdir)/src/modules/glm/CSparse		\
-I$(top_srcdir)/src/modules/glm/SSparse/UFconfig \
-I$(top_srcdir)/src/modules/glm/SSparse/CHOLMOD/Include

glm_la_LIBADD = samplers/glmsampler.la CSparse/csparse.la	\
	SSparse/ssparse.la $(top_builddir)/src/lib/libjags.la	\
	$(top_builddir)/src/jrmath/libjrmath.la @LAPACK_LIBS@	\
	@BLAS_LIBS@

glm_la_LDFLAGS = -module -avoid-version
if WINDOWS
glm_la_LDFLAGS += -no-undefined
endif
