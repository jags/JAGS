#include <config.h>
#include <util/logical.h>
#include "Subtract.h"

using std::vector;

namespace jags {
namespace base {

    Subtract::Subtract () : Infix ("-")
    {
    }
    
    double
    Subtract::evaluate(vector<double const *> const &args) const
    {
	return *args[0] - *args[1];
    }
    
    bool Subtract::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

    bool Subtract::isAdditive(vector<bool> const &mask,
			      vector<bool> const &fix) const
    {
	if (!mask[0] || mask[1]) return false;
	if (fix.empty()) {
	    return true;
	}
	else {
	    return (!fix[0] && fix[1]);
	}
    }
    
    bool 
    Subtract::isScale(vector<bool> const &mask, vector<bool> const &) const
    {
	return allTrue(mask);
    }
    
    bool 
    Subtract::isLinear(vector<bool> const &, vector<bool> const &) const
    {
	return true;
    }

}}
