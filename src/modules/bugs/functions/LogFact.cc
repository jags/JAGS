#include <config.h>
#include <sarray/SArray.h>
#include "LogFact.h"

#include <Rmath.h>

using std::vector;

LogFact::LogFact ()
  : ScalarFunc ("logfact", 1)
{
}

double LogFact::eval (vector <SArray const *> const &args) const
{
  return lgammafn (*args[0]->value() + 1);
}

bool LogFact::checkParameterValue (vector <SArray const *> const &args) const
{
  return *args[0]->value()  > -1;
}
