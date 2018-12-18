#include <config.h>
#include "Divide.h"

using std::vector;

namespace jags {
namespace base {

    Divide::Divide () : Infix ("/")
    {
    }

    double Divide::evaluate(vector <double const *> const &args) const
    {
	return *args[0] / *args[1];
    }

    bool Divide::isDifferentiable(unsigned long i) const
    {
	return i < 2;
    }
    
    double Divide::gradient(vector<double const *> const &args,
			    unsigned long i) const
    {
	if (i == 1) {
	    double x = *args[1];
	    return -1/(x*x);
	}
	else {
	    //i == 0
	    return *args[0];
	}
    }
    
    bool Divide::checkParameterValue(vector<double const*> const &args) const
    {
	return *args[1] != 0;
    }

    bool 
    Divide::isScale(vector<bool> const &mask, vector<bool> const &fix) const
    {
	if (mask[1])
	    return false; //No reciprocal terms

	return fix.empty() || fix[1];
    }

    bool 
    Divide::isPower(vector<bool> const &, vector<bool> const &) const
    {
	return true;
    }
	
}}
