#include <config.h>
#include "Combine.h"

#include <algorithm>
#include <numeric>

#include <util/dim.h>
#include <util/logical.h>

using std::vector;
using std::copy;
using std::accumulate;

namespace jags {
    namespace bugs {

	Combine::Combine ()
	    : VectorFunction ("c", 0)
	{
	}

	void Combine::evaluate(double *value, 
			       vector <double const *> const &args,
			       vector<unsigned long> const &lengths) const
	{
	    for (unsigned long i = 0; i < args.size(); ++i) {
		value = copy(args[i], args[i] + lengths[i], value);
	    }
	}

	bool Combine::isDifferentiable(unsigned long i) const
	{
	    return true;
	}
	
	void Combine::gradient(double *grad,
			       vector <double const *> const &args,
			       vector<unsigned long> const &lengths,
			       unsigned long i) const
	{
	    if (i >= args.size()) return;
	    
	    unsigned long P = accumulate(lengths.begin(), lengths.end(), 0U);
	    unsigned long Q = lengths[i];
	    
	    unsigned long p = 0;
	    for (unsigned long j = 0; j < i; ++j) {
		p += lengths[j];
	    }
	    for (unsigned long q = 0; q < Q; ++q, ++p) {
		grad[p + P * q] += 1;
	    }
	}
	
	unsigned long Combine::length (vector<unsigned long> const &lens,
				       vector<double const *> const &) const
	{
	    return accumulate(lens.begin(), lens.end(), 0U);
	}

	bool Combine::isDiscreteValued(vector<bool> const &mask) const
	{
	    return allTrue(mask);
	}

	bool Combine::isAdditive(vector<bool> const &mask, 
				 vector<bool> const &fixed) const
	{
	    //The Combine function behaves like an aggregate node
	    //and so has the same rules for preserving additive functions
	    //i.e. only one argument may be additive. 
	    bool found = false;
	    for (unsigned long i = 0; i < mask.size(); ++i) {
		if (mask[i]) {
		    if (found) return false;
		    else found = true;
		}
		if (!fixed.empty() && !fixed[i]) {
		    return false;
		}
	    }
	    return found;
	}
	
	bool Combine::isScale(vector<bool> const &mask, 
			      vector<bool> const &) const
	{
	    //The Combine function behaves like an aggregate node
	    //and so has the same rules for preserving scale functions
	    return allTrue(mask);
	}

	bool Combine::isLinear(vector<bool> const &, 
			       vector<bool> const &) const
	{
	    //The Combine function behaves like an aggregate node
	    //and so has the same rules for preserving linear functions
	    return true;
	}

	bool Combine::checkParameterLength(vector<unsigned long> const &)
	    const
	{
	    return true;
	}
    }

}
