SUBDIRS = distributions functions matrix 

jagsmod_LTLIBRARIES = msm.la

msm_la_SOURCES = msm.cc

msm_la_CPPFLAGS = -I$(top_srcdir)/src/include

msm_la_LIBADD = distributions/msmdist.la              \
        functions/msmfunc.la                           \
        matrix/msmmatrix.la                            \
        $(top_builddir)/src/jrmath/libjrmath.la		\
        $(top_builddir)/src/lib/libjags.la 		\
        @LAPACK_LIBS@ @BLAS_LIBS@ 

msm_la_LDFLAGS = -module -avoid-version
if WINDOWS
msm_la_LDFLAGS += -no-undefined
endif
