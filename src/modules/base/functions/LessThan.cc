#include <config.h>
#include "LessThan.h"

#include <cmath>
#include <cfloat>

using std::vector;

namespace jags {
namespace base {

LessThan::LessThan () : Infix ("<")
{
}

double LessThan::evaluate(vector<double const *> const &args) const
{
  return *args[0] < *args[1];
}

bool LessThan::isDiscreteValued(vector<bool> const &) const
{
  return true;
}

}}
