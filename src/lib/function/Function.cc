#include <config.h>
#include <function/Function.h>

using std::vector;
using std::string;

namespace jags {

Function::Function (string const &name, unsigned long npar)
  : _name (name), _npar (npar)
{
}

Function::~Function ()
{
  //virtual destructor
}

string const & Function::name () const
{
  return _name;
}

string Function::alias() const
{
    return "";
}

    unsigned long Function::npar() const
    {
	return _npar;
    }


bool Function::isDiscreteValued(vector<bool> const &) const
{
  return false;
}

bool Function::isScale(vector<bool> const &, vector<bool> const &) const
{
    return false;
}

bool Function::isLinear(vector<bool> const &mask,
			vector<bool> const &isfixed) const
{
    return isScale(mask, isfixed);
}

bool Function::isAdditive(vector<bool> const &, vector<bool> const &) const
{
    return false;
}
    
bool 
Function::isPower(vector<bool> const &, vector<bool> const &) const
{
    return false;
}

string Function::deparse(vector<string> const &par) const
{
    string name = this->name();
    name.append("(");
    for (unsigned long i = 0; i < par.size(); ++i) {
	if (i > 0) {
	    name.append(",");
	}
	name.append(par[i]);
    }
    name.append(")");

    return name;
}

    bool Function::checkParameterDiscrete(vector<bool> const &) const
    {
	return true;
    }

    bool Function::checkParameterFixed(std::vector<bool> const &) const
    {
	return true;
    }

} //namespace jags
