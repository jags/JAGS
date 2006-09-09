#include <config.h>
#include <rng/RNG.h>

using std::string;

RNG::RNG(string const &name)
  : _name(name)
{}

RNG::~RNG()
{}

#define i2_32m1 2.328306437080797e-10/* = 1/(2^32 - 1) */

double RNG::fixup(double x)
{
    /* ensure 0 and 1 are never returned */
    if(x <= 0.0) return 0.5*i2_32m1;
    if((1.0 - x) <= 0.0) return 1.0 - 0.5*i2_32m1;
    return x;
}

string const &RNG::name() const
{
   return _name;
}
