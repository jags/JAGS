BUILT_SOURCES = parser.h 
AM_YFLAGS = -d

libexec_PROGRAMS = jags-terminal

jags_terminal_SOURCES = parser.yy scanner.ll ReadData.cc 
## We list the modules here so that ltdl can emulate dynamic loading
## on platforms that do not support it, or where the user has chosen
## static linking
jags_terminal_LDADD = @LIBLTDL@ ${top_builddir}/src/lib/libjags.la	\
-dlopen ${top_builddir}/src/modules/base/basemod.la \
-dlopen ${top_builddir}/src/modules/bugs/bugs.la \
-dlopen ${top_builddir}/src/modules/dic/dic.la \
-dlopen ${top_builddir}/src/modules/mix/mix.la \
-dlopen ${top_builddir}/src/modules/msm/msm.la 
jags_terminal_CPPFLAGS= -I$(top_srcdir)/src/include @INCLTDL@

LEX_OUTPUT_ROOT=lex.zz

noinst_HEADERS = ReadData.h 

## The shell script is not required under Windows, so we do not
## build or install it. Instead, we install a batch file 

if WINDOWS

dist_bin_SCRIPTS = jags.bat

else

bin_SCRIPTS = jags

# The edit command substitutes paths and filenames into the jags script.
# Note that jterm (the name of the installed version of jags-terminal)
# is the odd one out as it is a shell variable not
# a Makefile variable. In particular, this means that the sed script
# must be in double quotes so that it is expanded properly.

edit = sed					\
        -e 's,@jagsmoddir\@,$(jagsmoddir),g'	\
        -e 's,@libexecdir\@,$(libexecdir),g'	\
        -e 's,@libdir\@,$(libdir),g'		\
	-e "s,@jterm\@,$${jterm},g"

jags: Makefile $(srcdir)/jags.in
	rm -f jags jags.tmp
	jterm=`echo jags-terminal | sed -e '$(program_transform_name)'`; \
	$(edit) $(srcdir)/jags.in > jags.tmp
	chmod +x jags.tmp
	mv jags.tmp jags

endif


## We need to give Automake some additional hints for distribution
## and cleaning of the jags shell script and its source.

EXTRA_DIST = jags.in
CLEANFILES = jags
