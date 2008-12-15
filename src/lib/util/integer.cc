#include <config.h>
#include <stdexcept>
#include <cmath>
#include <cfloat>
#include <climits>
#include <string>

using std::runtime_error;
using std::fabs;
using std::string;

static const double eps = 16 * DBL_EPSILON;

int asInteger(double fval)
{
    if (fval >= INT_MAX || fval <= INT_MIN) {
	throw runtime_error(string("double value out of range for conversion to int"));
    }
    int ival;
    if (fval > 0) {
	ival = static_cast<int>(fval + eps);
    }
    else {
	ival = static_cast<int>(fval - eps);
    }
    return ival;
}

int checkInteger(double fval, bool &flag)
{
    if (fval >= INT_MAX) {
	flag = false;
	return INT_MAX;
    }
    else if (fval <= INT_MIN) {
	flag = false;
	return INT_MIN;
    }
    else {
	int ival = asInteger(fval);
	flag = fabs(fval - ival) < eps;
	return ival;
    }
}

