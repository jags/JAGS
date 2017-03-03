if CANCHECK

# Rules for the test code (use `make check` to execute)
TESTS = base bugs mix glm
check_PROGRAMS = $(TESTS)


## Base module

base_SOURCES = base.cc 
base_CXXFLAGS = $(CPPUNIT_CFLAGS)
base_LDFLAGS = $(CPPUNIT_LIBS)

base_LDADD = $(top_builddir)/src/modules/base/libbasetest.la

base_CPPFLAGS = -I$(top_srcdir)/src/include	\
	-I$(top_srcdir)/src/modules

## Bugs module

bugs_SOURCES = bugs.cc 
bugs_CXXFLAGS = $(CPPUNIT_CFLAGS)
bugs_LDFLAGS = $(CPPUNIT_LIBS)

bugs_LDADD = $(top_builddir)/src/modules/bugs/libbugstest.la

bugs_CPPFLAGS = -I$(top_srcdir)/src/include	\
	-I$(top_srcdir)/src/modules

## Mix module

mix_SOURCES = mix.cc 
mix_CXXFLAGS = $(CPPUNIT_CFLAGS)
mix_LDFLAGS = $(CPPUNIT_LIBS)

mix_LDADD = $(top_builddir)/src/modules/mix/libmixtest.la

mix_CPPFLAGS = -I$(top_srcdir)/src/include	\
	-I$(top_srcdir)/src/modules

## Glm module

glm_SOURCES = glm.cc 
glm_CXXFLAGS = $(CPPUNIT_CFLAGS)
glm_LDFLAGS = $(CPPUNIT_LIBS)

glm_LDADD = $(top_builddir)/src/modules/glm/libglmtest.la

glm_CPPFLAGS = -I$(top_srcdir)/src/include	\
	-I$(top_srcdir)/src/modules

endif
