#include <config.h>
#include "Sort.h"

#include <algorithm>

#include <util/dim.h>
#include <util/logical.h>

using std::vector;
using std::sort;

namespace jags {
namespace bugs {

    Sort::Sort ()
	: VectorFunction ("sort", 1)
    {
    }

    void Sort::evaluate (double *value, vector <double const *> const &args,
			 vector<unsigned long> const &lengths) const
    {
	for (unsigned long i = 0; i < lengths[0]; ++i) {
	    value[i] = args[0][i];
	}
	sort(value, value + lengths[0]);
    }

    unsigned long Sort::length (vector<unsigned long> const &parlengths,
				vector<double const *> const &) const
    {
	return parlengths[0];
    }

    bool Sort::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

}}
