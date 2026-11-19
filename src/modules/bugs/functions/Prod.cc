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

    bool Prod::hasGradient(unsigned long i) const
    {
	return true;
    }
    
    void Prod::gradient(double *grad,
			vector<double const *> const &args,
			vector<unsigned long> const &lengths,
			unsigned long i) const
    {
	if (i >= args.size()) return;
	
	double y = scalarEval(args, lengths);

	for (unsigned long j = 0; j < lengths[i]; ++j) {	    
	    grad[j] += y/args[i][j];
	}
    }

    bool Prod::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

}}
