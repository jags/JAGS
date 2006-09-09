INCLUDES = -I $(top_srcdir)/src/include 

noinst_LTLIBRARIES = rmath.la

rmath_la_SOURCES = callbacks.cc \
mlutils.c stirlerr.c bd0.c choose.c gamma.c lgamma.c toms708.c \
d1mach.c i1mach.c fsign.c imax2.c imin2.c \
fmin2.c fmax2.c ftrunc.c \
chebyshev.c lgammacor.c log1p.c lbeta.c \
rbeta.c dbeta.c pbeta.c qbeta.c \
rbinom.c dbinom.c pbinom.c qbinom.c \
rchisq.c dchisq.c pchisq.c qchisq.c \
rexp.c dexp.c pexp.c qexp.c \
rgamma.c dgamma.c pgamma.c qgamma.c \
rhyper.c dhyper.c phyper.c qhyper.c \
rlnorm.c dlnorm.c plnorm.c qlnorm.c \
rlogis.c dlogis.c plogis.c qlogis.c \
rnbinom.c dnbinom.c pnbinom.c qnbinom.c \
rnorm.c dnorm.c pnorm.c qnorm.c \
rpois.c dpois.c ppois.c qpois.c \
rt.c dt.c pt.c qt.c \
runif.c dunif.c punif.c qunif.c \
rweibull.c dweibull.c pweibull.c qweibull.c 

rmath_la_LDFLAGS = -no-undefined -module -avoid-version

noinst_HEADERS = nmath.h dpq.h Rmath.h

EXTRA_DIST = README
