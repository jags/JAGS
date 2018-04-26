#include <config.h>
#include "SD.h"

#include <cmath>

using std::vector;
using std::sqrt;

namespace jags {
namespace bugs {

    SD::SD ()
	: ScalarVectorFunction ("sd", 1)
    {
    }

    double SD::scalarEval (vector<double const *>const &args,
			   vector<unsigned long> const &lengths) const
    {
	double svalue = 0;
	if (lengths[0] > 1) {
	    double sum = 0;
	    for (unsigned long i = 0; i < lengths[0]; i++) {
		sum += args[0][i];
	    }
	    double mean = sum / lengths[0];

	    double var = 0;
	    for (unsigned long i = 0; i < lengths[0]; i++) {
		var += (args[0][i] - mean) * (args[0][i] - mean);
	    }
	    svalue = sqrt (var / (lengths[0] - 1));
	}
	return svalue;
    }

    bool SD::checkParameterLength (vector<unsigned long> const &lengths) const
    {
	return lengths[0] >= 2;
    }

}}
