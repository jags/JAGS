#include <config.h>
#include <util/logical.h>
#include "Max.h"

#include <algorithm>

using std::vector;
using std::max_element;
using std::max;

namespace jags {
namespace bugs {

    Max::Max ()
	: ScalarVectorFunction("max", 0)
    {
    }

    double Max::scalarEval(vector<double const *> const &args,
			   vector<unsigned long> const &lengths) const
    {
	double ans = *max_element(args[0], args[0] + lengths[0]);
	for (unsigned long i = 1; i < args.size(); ++i) {
	    double maxi = *max_element(args[i], args[i] + lengths[i]);
	    ans = max(ans, maxi);
	}
	return ans;
    }

    bool Max::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

}}
