#include <config.h>
#include <function/Function.h>
#include <sarray/SArray.h>

using std::vector;
using std::string;

Function::Function (string const &name, unsigned int npar)
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

bool Function::checkParameterValue (vector <SArray const *> const &args) const
{
  return true;
}

#include <iostream>
bool Function::checkParameterLength (unsigned int npar) const
{
  return _npar == 0 ||  npar == _npar;
}

bool Function::isDiscreteValued(vector<bool> const &mask) const
{
  return false;
}

bool Function::isLinear(vector<bool> const &mask,
			vector<bool> const &isfixed) const
{
  return false;
}

bool Function::isScale(unsigned int index,
		       vector<bool> const &isfixed) const
{
    return false;
}
