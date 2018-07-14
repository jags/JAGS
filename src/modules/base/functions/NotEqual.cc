#include <config.h>
#include "NotEqual.h"

#include <cmath>
#include <cfloat>

using std::vector;

namespace jags {
namespace base {

NotEqual::NotEqual () : Infix ("!=")
{
}

double NotEqual::evaluate(vector<double const *> const &args) const
{
  return *args[0] != *args[1];
}

bool NotEqual::isDiscreteValued(vector<bool> const &) const
{
  return true;
}

}}
