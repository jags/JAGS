#include <config.h>

#include "QFunction.h"

#include <distribution/DistScalarRmath.h>

using std::vector;
using std::string;

namespace bugs {

    QFunction::QFunction(DistScalarRmath *dist)
	: DPQFunction(string("q") + dist->name().substr(1), dist)
    {}
    
    double QFunction::evaluate(vector<double const *> const &args) const
    {
	double x = *args[0];
	vector<double const *> param(args.size() - 1);
	for (unsigned int i = 1; i < args.size(); ++i) {
	    param[i-1] = args[i];
	}
	
	return dist()->q(x, param, true, false);
    }

    bool QFunction::checkParameterValue(vector<double const*> const &args) const
    {
	double x = *args[0];
	if (x <= 0 || x >= 1)
	    return false;

	return checkArgs(args);
    }

}