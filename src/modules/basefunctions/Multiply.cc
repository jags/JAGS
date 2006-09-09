#include <config.h>
#include <sarray/SArray.h>
#include <sarray/util.h>
#include "Multiply.h"

using std::vector;

Multiply::Multiply ()
  : ScalarFunc ("*", 0)
{
}

double
Multiply::eval (vector < SArray const *>const &args) const
{
  double val = *args[0]->value();
  /* We adopt the convention that any argument multiplied by zero
     is zero. This is required to do calculations like x*log(x)
     when x is zero.
  */
  if (val == 0) {
    return 0;
  }
  for (unsigned int i = 1; i < args.size(); ++i) {
    double vali = *args[i]->value();
    if (vali == 0) {
      return 0;
    }
    val *= vali;
  }
  return val;
}

bool Multiply::isDiscreteValued(vector<bool> const &mask) const
{
  return allTrue(mask);
}

bool Multiply::isLinear(vector<bool> const &mask,
                        vector<bool> const &fix) const
{
    //Test for quadratic or higher terms
    unsigned long nfactor = 0;
    for (unsigned long i = 0; i < mask.size(); ++i) {
	nfactor += mask[i];
    }
    if (nfactor > 1) {
	return false;
    }

    if (!fix.empty()) {
	for (unsigned int i = 0; i < fix.size(); ++i) {
	    if (!mask[i] && !fix[i]) {
		return false;
	    }
	}
    }

    return true;
}

bool Multiply::isScale(unsigned int index, vector<bool> const &fix) const
{
    if (fix.empty()) {
	return true;
    }
    else {
	for (unsigned int i = 0; i < fix.size(); ++i) {
	    if (i != index && !fix[i]) {
		return false;
	    }
	}
    }

    return true;
}
