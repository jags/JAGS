#include <config.h>
#include "GreaterOrEqual.h"

#include <cmath>
#include <cfloat>

using std::vector;

namespace jags {
namespace base {

GreaterOrEqual::GreaterOrEqual () : Infix (">=")
{
}

double GreaterOrEqual::evaluate(vector<double const *> const &args) const
{
  return *args[0] >= *args[1];
}

bool GreaterOrEqual::isDiscreteValued(vector<bool> const &) const
{
  return true;
}

}}
