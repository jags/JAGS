noinst_LTLIBRARIES = libcompiler.la
AM_YFLAGS = -d
BUILT_SOURCES = parser.h

libcompiler_la_CPPFLAGS = -I$(top_srcdir)/src/include

libcompiler_la_SOURCES = parser.yy scanner.ll Compiler.cc		\
LogicalFactory.cc ParseTree.cc Counter.cc CounterTab.cc			\
ConstantFactory.cc MixtureFactory.cc MixCompiler.cc NodeFactory.cc

noinst_HEADERS = parser.h parser_extra.h remap.h MixCompiler.h
