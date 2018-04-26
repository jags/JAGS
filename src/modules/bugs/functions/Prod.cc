#include <config.h>
#include <util/logical.h>
#include "Prod.h"

using std::vector;

namespace jags {
namespace bugs {

    Prod::Prod () : ScalarVectorFunction("prod", 0)
    {
    }

    double Prod::scalarEval(vector <double const *> const &args,
			    vector<unsigned long> const &lengths) const
    {
	double value = 1;
	for (unsigned long j = 0; j < args.size(); ++j) {
	    for (unsigned long i = 0; i < lengths[j]; ++i) {
		value *= args[j][i];
	    }
	}
	return value;
    }

    bool Prod::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

}}
