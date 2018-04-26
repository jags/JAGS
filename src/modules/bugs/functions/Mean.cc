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

    bool Mean::isScale(vector<bool> const &mask, vector<bool> const &fix) const
    {
	return true;
    }

}}
