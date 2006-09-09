#include <config.h>
#include <sarray/SArray.h>
#include "SD.h"

#include <cmath>

using std::vector;
using std::sqrt;

SD::SD ()
  : ScalarFunc ("sd", 1)
{
}

double
SD::eval (vector < SArray const *>const &args) const
{
  double const *arg0 = args[0]->value ();
  long len = args[0]->length ();

  double svalue = 0;
  if (len > 1)
    {
      long i;
      double sum = 0;
      for (i = 0; i < len; i++)
	{
	  sum += arg0[i];
	}
      double mean = sum / len;

      double var = 0;
      for (i = 0; i < len; i++)
	{
	  var += (arg0[i] - mean) * (arg0[i] - mean);
	}
      svalue = sqrt (var / len);
    }
  return svalue;
}

bool SD::checkParameterDim (vector <Index> const &dims) const
{
  return dims[0][0] >= 2;
}
