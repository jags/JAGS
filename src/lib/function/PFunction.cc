#include <config.h>
#include <function/PFunction.h>
#include <distribution/RScalarDist.h>

using std::vector;
using std::string;

namespace jags {

    PFunction::PFunction(RScalarDist const *dist)
	: DPQFunction(string("p") + dist->name().substr(1), dist)
    {}
    
    double PFunction::evaluate(vector<double const *> const &args) const
    {
	double x = *args[0];
	vector<double const *> param(args.size() - 1);
	for (unsigned long i = 1; i < args.size(); ++i) {
	    param[i-1] = args[i];
	}
	
	return dist()->p(x, param, true, false);
    }

    bool 
    PFunction::checkParameterValue(vector<double const *> const &args) const
    {
	if (dist()->discrete()) {
	    double x = *args[0];
	    if (x != static_cast<int>(x))
		return false; //FIXME: use checkInteger
	}

	return checkArgs(args);
    }

}
