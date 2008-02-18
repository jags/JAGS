BUILT_SOURCES = parser.h 
AM_YFLAGS = -d

libexec_PROGRAMS = jags-terminal

jags_terminal_SOURCES = parser.yy scanner.ll ReadData.cc 
jags_terminal_LDADD = @LIBLTDL@ ${top_builddir}/src/lib/libjags.la
jags_terminal_CPPFLAGS= -I$(top_srcdir)/src/include @INCLTDL@

LEX_OUTPUT_ROOT=lex.zz

noinst_HEADERS = ReadData.h 

## The shell script is not required under Windows, so we do not
## build or install it. Instead, we install a batch file 

if WINDOWS

dist_bin_SCRIPTS = jags.bat

else

bin_SCRIPTS = jags

edit = sed					\
        -e 's,@jagsmoddir\@,$(jagsmoddir),g'	\
        -e 's,@libexecdir\@,$(libexecdir),g'	\
        -e 's,@libdir\@,$(libdir),g'

jags: Makefile $(srcdir)/jags.in
	rm -f jags jags.tmp
	$(edit) $(srcdir)/jags.in >jags.tmp
	chmod +x jags.tmp
	mv jags.tmp jags

endif

## We need to give Automake some additional hints for distribution
## and cleaning of the jags shell script and its source.

EXTRA_DIST = jags.in
CLEANFILES = jags
