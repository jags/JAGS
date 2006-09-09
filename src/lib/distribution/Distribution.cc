#include <config.h>
#include <distribution/Distribution.h>

using std::vector;
using std::string;

Distribution::Distribution(string const &name, unsigned int npar, 
			   bool canbound, bool discrete)
    : _name(name), _npar(npar), _canbound(canbound), _discrete(discrete)
{
}

Distribution::~Distribution()
{
}

string const &Distribution::name() const
{
  return _name;
}

bool 
Distribution::isDiscreteValued() const
{
  return _discrete;
}

unsigned int Distribution::npar() const
{
  return _npar;
}

SArray const *
Distribution::lowerBound(std::vector<SArray const *> const &args) const
{
  return _canbound ? args[_npar] : 0;
}

SArray const *
Distribution::upperBound(std::vector<SArray const *> const &args) const
{
  return _canbound ? args[_npar + 1] : 0;
}

bool Distribution::canBound() const
{
  return _canbound;
}

bool Distribution::isDeterministic() const
{
  return false;
}
