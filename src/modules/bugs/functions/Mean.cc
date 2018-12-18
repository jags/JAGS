#include <config.h>
#include "Mean.h"

using std::vector;

namespace jags {
namespace bugs {

    Mean::Mean ()
	: ScalarVectorFunction ("mean", 1)
    {
    }

    double Mean::scalarEval (vector<double const*> const &args,
			     vector<unsigned long> const &lengths) const
    {
	double svalue = 0;
	for (unsigned long i = 0; i < lengths[0]; i++) {
	    svalue += args[0][i];
	}
	svalue /= lengths[0];
	return svalue;
    }

    bool Mean::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    void Mean::gradient(double *grad, vector<double const *> const &args,
			vector<unsigned long> const &lengths,
			unsigned long i) const
    {
	unsigned long N = lengths[0];
	double y = 1.0/N;
	for (unsigned long j = 0; j < N; ++j) {
	    grad[j] += y;
	}
    }

    bool Mean::isScale(vector<bool> const &, vector<bool> const &) const
    {
	return true;
    }

}}
