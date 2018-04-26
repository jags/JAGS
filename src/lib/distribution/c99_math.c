/* Remap some C99 functions for conformance with C++98 */
/* FIXME: Not needed when we use C++11 */

#include <math.h>

double jags_expm1(double x) {
   return expm1(x);
}

double jags_log1p(double x) {
   return log1p(x);
}
