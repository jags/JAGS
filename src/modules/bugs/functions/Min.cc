#include <config.h>
#include <util/logical.h>
#include "Min.h"

#include <algorithm>

using std::min_element;
using std::min;
using std::vector;

namespace jags {
namespace bugs {

    Min::Min ()
	: ScalarVectorFunction ("min", 0)
    {
    }

    double Min::scalarEval(vector<double const *> const &args,
			   vector<unsigned long> const &lengths) const
    {
	double ans = *min_element(args[0], args[0] + lengths[0]);
	for (unsigned long i = 1; i < args.size(); ++i) {
	    double mini = *min_element(args[i], args[i] + lengths[i]);
	    ans = min(ans, mini);
	}
	return ans;
    }

    bool Min::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

}}
