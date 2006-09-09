#include <config.h>
#include <sarray/SArray.h>
#include "LogGam.h"

#include <Rmath.h>

using std::vector;

LogGam::LogGam ()
  : ScalarFunc ("loggam", 1)
{
}

double LogGam::eval (vector <SArray const *> const &args) const
{
  return lgammafn (*args[0]->value ());
}

bool LogGam::checkParameterValue (vector <SArray const *> const &args) const
{
  return *args[0]->value() > 0;
}
