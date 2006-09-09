BUILT_SOURCES = parser.h
AM_YFLAGS = -d

jagsbindir = @jagshomedir@/bin
jagsbin_PROGRAMS = jags.bin 
jags_bin_SOURCES = parser.yy scanner.ll ReadData.cc 

bin_PROGRAMS = jags
jags_SOURCES = jags.in

noinst_HEADERS = ReadData.h 

jags_bin_LDADD = @LIBLTDL@ ${top_builddir}/src/lib/libjags.la
jags_bin_CPPFLAGS= -I $(top_srcdir)/src/include @INCLTDL@

LEX_OUTPUT_ROOT=lex.zz

edit = sed					\
        -e 's,@libdir\@,$(libdir),g'		\
        -e 's,@jagsmoddir\@,$(jagsmoddir),g'	\
        -e 's,@jagsbindir\@,$(jagsbindir),g'	\
        -e 's,@jagshomedir\@,$(jagshomedir),g'

jags: Makefile $(srcdir)/jags.in
	rm -f jags jags.tmp
	$(edit) $(srcdir)/jags.in >jags.tmp
	chmod +x jags.tmp
	mv jags.tmp jags
