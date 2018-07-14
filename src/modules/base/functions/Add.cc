#include <config.h>
#include <util/logical.h>
#include "Add.h"

using std::vector;

namespace jags {
namespace base {

    Add::Add () : Infix ("+", 0)
    {
    }
    
    double Add::evaluate(vector <double const *> const &args) const
    {
	double out = *args[0];
	for (unsigned int i = 1; i < args.size(); ++i) {
	    out += *args[i];
	}
	return out;
    }
    
    bool Add::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

    bool
    Add::isAdditive(vector<bool> const &mask, vector<bool> const &fix) const
    {
	//There should be exactly one argument that is an additive function
	bool found = false; //Have we found it yet?
	for (unsigned int i = 0; i < mask.size(); ++i) {
	    if (mask[i]) {
		if (found) return false;
		else found = true;
	    }
	    else if (!fix.empty() && !fix[i]) {
		return false;
	    }
	}
	return found;
    }
    
    bool Add::isLinear(vector<bool> const &, vector<bool> const &) const
    {
	return true;
    }

    bool Add::isScale(vector<bool> const &mask, vector<bool> const &) const
    {
	return allTrue(mask);
    }

}}
