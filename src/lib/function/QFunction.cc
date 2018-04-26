#include <config.h>
#include <function/QFunction.h>
#include <distribution/RScalarDist.h>

using std::vector;
using std::string;

namespace jags {

    QFunction::QFunction(RScalarDist const *dist)
	: DPQFunction(string("q") + dist->name().substr(1), dist)
    {}
    
    double QFunction::evaluate(vector<double const *> const &args) const
    {
	double x = *args[0];
	vector<double const *> param(args.size() - 1);
	for (unsigned long i = 1; i < args.size(); ++i) {
	    param[i-1] = args[i];
	}
	
	return dist()->q(x, param, true, false);
    }

    bool QFunction::checkParameterValue(vector<double const*> const &args) const
    {
	return checkArgs(args);
    }

}
