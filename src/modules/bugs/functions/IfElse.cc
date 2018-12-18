#include <config.h>
#include "IfElse.h"
#include <util/nainf.h>

using std::vector;

namespace jags {
namespace bugs {

    IfElse::IfElse () : ScalarFunction("ifelse", 3)
    {
    }
    
    double IfElse::evaluate(vector<double const *> const &args) const
    {
	return *args[0] ? *args[1] : *args[2];
    }

    bool IfElse::isDifferentiable(unsigned long i) const
    {
	return i == 1 || i == 2;
    }
    
    double IfElse::gradient(vector<double const *> const &args,
			    unsigned long i) const
    {
	if (i == 1) {
	    return *args[0] ? 1 : 0;
	}
	else if (i == 2) {
	    return *args[0] ? 0 : 1;
	}
	else {
	    return 0;
	}
    }

    bool IfElse::isDiscreteValued(vector<bool> const &mask) const
    {
	return mask[1] && mask[2];
    }

    bool IfElse::isScale(vector<bool> const &mask, 
			 vector<bool> const &fixed) const
    {
	//Don't even bother with the case where the predicate is a
	//scale function. False negatives are allowed and this is a
	//silly case
	if (mask[0]) return false;
	
	if (fixed.empty()) {
	    //result is a scale function if both arguments are
	    return mask[1] && mask[2];
	}
	else {
	    //..but is never a fixed scale function
	    return false;
	}
    }

    bool IfElse::isLinear(vector<bool> const &mask, 
			  vector<bool> const &fixed) const
    {
	if (mask[0]) return false;
	
	if (fixed.empty()) {
	    //result is a linear function if either argument is
	    return mask[1] || mask[2];
	}
	else {
	    //..but is never a fixed scale function
	    return false;
	}
    }

    bool IfElse::isAdditive(vector<bool> const &mask, 
			    vector<bool> const &fixed) const
    {
	if (mask[0]) return false;

	if (fixed.empty()) {
	    //result is an addtive function if both arguments are
	    return mask[1] && mask[2];
	}
	else {
	    //..but it is never fixed
	    return false;
	}
	
    }
    


}}
