ACLOCAL_AMFLAGS = -I m4
SUBDIRS = m4 libltdl src etc win32
EXTRA_DIST = Doxyfile

win32-installer:
	make prefix=`pwd`/win32/inst
	make install prefix=`pwd`/win32/inst
	(cd win32; for f in `ls lapack/*.dll`; do \
	    ${INSTALL} $${f} inst/bin; \
	  done; \
	  for f in `ls lapack/*.dll.a`; do \
	    ${INSTALL} $${f} inst/lib; \
	  done;) 
	makensis -DVERSION=${PACKAGE_VERSION} win32/jags.nsi
# The edit command substitutes paths and filenames into the jags script.
# Note that jterm (the name of the installed version of jags-terminal)
# is the odd one out as it is a shell variable not
# a Makefile variable. In particular, this means that the sed script
# must be in double quotes so that it is expanded properly.


