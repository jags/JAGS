#include <config.h>
#include <util/dim.h>
#include "Rank.h"

#include <algorithm>

using std::vector;
using std::stable_sort;

static inline bool lt_doubleptr (double const *arg1, double const *arg2) {
  return *arg1 < *arg2;
}

namespace jags {
namespace bugs {

    Rank::Rank ()
	: VectorFunction ("rank", 1)
    {
    }

    void Rank::evaluate(double *value, vector<double const *> const &args,
			vector<unsigned long> const &lengths) const
    {
	unsigned long N = lengths[0];

	//Create a vector of pointers to the elements of arg and sort it
	vector<double const*> argptrs(N);
	for (unsigned long i = 0; i < N; ++i) {
	    argptrs[i] = args[0] + i;
	}
	stable_sort(argptrs.begin(), argptrs.end(), lt_doubleptr);

	//Ranks can be inferred from the sorted vector of pointers
	for (unsigned long i = 0; i < N; ++i) {
	    value[argptrs[i] - args[0]] = i + 1;
	}
    }

    unsigned long Rank::length (vector<unsigned long> const &lengths,
				vector<double const *> const &) const
    {
	return lengths[0];
    }

    bool Rank::isDiscreteValued(vector<bool> const &) const
    {
	return true;
    }

}}
