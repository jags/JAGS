noinst_LTLIBRARIES = libbasefunctions.la

libbasefunctions_la_CPPFLAGS = -I$(top_srcdir)/src/include

libbasefunctions_la_SOURCES = Infix.cc Add.cc Subtract.cc Divide.cc	\
Multiply.cc Neg.cc Pow.cc And.cc Or.cc Not.cc Equal.cc NotEqual.cc	\
GreaterThan.cc GreaterOrEqual.cc LessThan.cc LessOrEqual.cc Seq.cc

noinst_HEADERS = Infix.h Add.h Subtract.h Divide.h Multiply.h Neg.h	\
Pow.h And.h Or.h Not.h Equal.h NotEqual.h GreaterThan.h			\
GreaterOrEqual.h LessThan.h LessOrEqual.h Seq.h

### Test library 

if CANCHECK
check_LTLIBRARIES = libbasefuntest.la
libbasefuntest_la_SOURCES = testbasefun.cc testbasefun.h
libbasefuntest_la_CPPFLAGS = -I$(top_srcdir)/src/include
libbasefuntest_la_CXXFLAGS = $(CPPUNIT_CFLAGS)
endif
