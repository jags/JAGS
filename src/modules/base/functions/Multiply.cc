#include <config.h>
#include <util/logical.h>
#include "Multiply.h"

using std::vector;

namespace jags {
namespace base {

    Multiply::Multiply () : Infix ("*", 0)
    {
    }

    double
    Multiply::evaluate(vector<double const *>const &args) const
    {
	double val = *args[0];
	/* We adopt the convention that any argument multiplied by zero
	   is zero. This is required to do calculations like x*log(x)
	   when x is zero.
	*/
	if (val == 0) {
	    return 0;
	}
	for (unsigned int i = 1; i < args.size(); ++i) {
	    double vali = *args[i];
	    if (vali == 0) {
		return 0;
	    }
	    val *= vali;
	}
	return val;
    }

    bool Multiply::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

    bool Multiply::isScale(vector<bool> const &mask,
			   vector<bool> const &fix) const
    {
	//Test for quadratic or higher terms
	unsigned long nfactor = 0;
	for (unsigned long i = 0; i < mask.size(); ++i) {
	    nfactor += mask[i];
	}
	if (nfactor > 1) {
	    return false;
	}

	if (!fix.empty()) {
	    for (unsigned int i = 0; i < fix.size(); ++i) {
		if (!mask[i] && !fix[i]) {
		    return false;
		}
	    }
	}

	return true;
    }

    bool Multiply::isPower(vector<bool> const &, 
			   vector<bool> const &) const
    {
	return true;
    }

}}
