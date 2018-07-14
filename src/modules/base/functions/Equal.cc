#include <config.h>
#include "Equal.h"

#include <cmath>
#include <cfloat>

using std::vector;

namespace jags {
namespace base {

Equal::Equal () : Infix ("==")
{
}

double Equal::evaluate(vector<double const *> const &args) const
{
  return *args[0] == *args[1];
}

bool Equal::isDiscreteValued(vector<bool> const &) const
{
  return true;
}

}}
