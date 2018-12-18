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

    bool Prod::isDifferentiable(unsigned long i) const
    {
	return true;
    }
    
    void Prod::gradient(double *grad,
			vector<double const *> const &args,
			vector<unsigned long> const &lengths,
			unsigned long i) const
    {
	if (i >= args.size()) return;
	
	double y = 1;
	for (unsigned long j = 0; j < args.size(); ++j) {
	    if (j == i) continue;
	    for (unsigned long k = 0; k < lengths[j]; ++k) {
		y *= args[j][k];
	    }
	}

	unsigned long N = lengths[i];
	for (unsigned long j = 0; j < N; ++j) {	    
	    double gj = y;
	    for (unsigned long k = 0; k < N; ++k) {
		if (k == j) continue;
		gj *= args[j][k];
	    }
	    grad[j] += gj;
	}
    }

    bool Prod::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

}}
