#include <config.h>
#include "LessOrEqual.h"

#include <cmath>
#include <cfloat>

using std::vector;

namespace jags {
namespace base {

LessOrEqual::LessOrEqual () : Infix ("<=")
{
}

double LessOrEqual::evaluate(vector<double const *> const &args) const
{
  return *args[0] <= *args[1];
}

bool LessOrEqual::isDiscreteValued(vector<bool> const &) const
{
  return true;
}

}}
