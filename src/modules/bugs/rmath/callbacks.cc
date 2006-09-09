#include <rng/RNG.h>
#include "nmath.h"

extern "C" {

  double unif_rand(RNG *rng)
  {
    return rng->uniform();
  }
  
  double exp_rand(RNG *rng)
  {
    return rng->exponential();
  }
  
  double norm_rand(RNG *rng)
  {
    return rng->normal();
  }
  
}
