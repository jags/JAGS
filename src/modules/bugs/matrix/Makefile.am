INCLUDES = -I $(top_srcdir)/src/include

noinst_LTLIBRARIES = bugsmatrix.la

bugsmatrix_la_SOURCES = matexp.cc matrix.cc

bugsmatrix_la_LDFLAGS = -no-undefined -module -avoid-version

bugsmatrix_la_LIBADD = @FLIBS@ @BLAS_LIBS@ @LAPACK_LIBS@

noinst_HEADERS = matexp.h lapack.h matrix.h
