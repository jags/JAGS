#include <config.h>
#include <util/dim.h>
#include "Order.h"

#include <algorithm>

using std::vector;
using std::stable_sort;

static bool lt_doubleptr (double const *arg1, double const *arg2) {
  return *arg1 < *arg2;
}

namespace jags {
namespace bugs {

    Order::Order ()
	: VectorFunction ("order", 1)
    {
    }

    void Order::evaluate(double *value, vector<double const *> const &args,
			vector<unsigned long> const &lengths) const
    {
	unsigned long N = lengths[0];

	//Create a vector of pointers to the elements of arg and sort it
	vector<double const *> argptrs(N);
	for (unsigned long i = 0; i < N; ++i) {
	    argptrs[i] = args[0] + i;
	}
	stable_sort(argptrs.begin(), argptrs.end(), lt_doubleptr);

	//Orders can be inferred from the sorted vector of pointers
	for (unsigned long i = 0; i < N; ++i) {
	    value[i] = argptrs[i] - args[0] + 1;
	}
    }

    unsigned long Order::length (vector<unsigned long> const &lengths,
				vector<double const *> const &parvalues) const
    {
	return lengths[0];
    }

    bool Order::isDiscreteValued(vector<bool> const &mask) const
    {
	return true;
    }

}}
