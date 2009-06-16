#include <config.h>
#include "Divide.h"

using std::vector;

namespace base {

    Divide::Divide () : Infix ("/")
    {
    }

    double Divide::evaluateScalar(vector <double const *> const &args) const
    {
	return *args[0] / *args[1];
    }

    bool Divide::checkScalarValue(vector<double const*> const &args) const
    {
	return *args[1] != 0;
    }

    bool 
    Divide::isScale(vector<bool> const &mask, vector<bool> const &fix) const
    {
	if (mask[1])
	    return false; //No reciprocal terms

	if (fix.empty())
	    return true;
	else
	    return fix[1];
    }

    bool 
    Divide::isPower(vector<bool> const &mask, vector<bool> const &fix) const
    {
	return true;
    }
	
}
