#include <config.h>
#include "GreaterThan.h"

#include <cmath>
#include <cfloat>

using std::vector;

namespace jags {
namespace base {

GreaterThan::GreaterThan () : Infix (">")
{
}

double GreaterThan::evaluate(vector<double const *> const &args) const
{
  return *args[0] > *args[1];
}

bool GreaterThan::isDiscreteValued(vector<bool> const &) const
{
  return true;
}

}}
