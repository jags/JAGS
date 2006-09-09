#include <config.h>
#include <sarray/SArray.h>
#include "Log.h"

#include <cmath>

using std::vector;
using std::log;


Log::Log ()
  : ScalarFunc ("log", 1)
{
}

double Log::eval (vector <SArray const *> const &args) const
{
  return log (*args[0]->value());
}

bool Log::checkParameterValue (vector <SArray const *> const &args) const
{
  return *args[0]->value() >= 0;
}
