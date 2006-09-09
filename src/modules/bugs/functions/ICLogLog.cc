#include <config.h>
#include <sarray/SArray.h>
#include "ICLogLog.h"

#include <cmath>

using std::vector;
using std::log;
using std::exp;

ICLogLog::ICLogLog(): InverseLinkFunc("icloglog", "cloglog")
{
}

double
ICLogLog::eval (vector <SArray const *> const &args) const
{
  return 1 - exp(-exp(*args[0]->value()));
}

double ICLogLog::link(double mu) const
{
  return log (-log (1 - mu));
}

double ICLogLog::gradLink(double mu) const
{
  return -1/((1-mu)*log(1-mu));
}

