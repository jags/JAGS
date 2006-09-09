#include <config.h>
#include <sarray/SArray.h>
#include "CLogLog.h"

#include <cmath>

using std::vector;
using std::log;

CLogLog::CLogLog ():ScalarFunc ("cloglog", 1)
{
}

double CLogLog::eval (vector <SArray const *> const &args) const
{
  return log (-log (1 - *(args[0]->value())));
}

bool CLogLog::checkParameterValue (vector <SArray const *> const &args) const
{
  double arg = *(args[0]->value ());
  return (arg > 0 && arg < 1);
}
