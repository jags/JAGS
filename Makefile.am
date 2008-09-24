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
