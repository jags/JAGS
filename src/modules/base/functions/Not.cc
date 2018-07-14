#include <config.h>
#include "Not.h"

using std::vector;
using std::string;

namespace jags {
namespace base {

Not::Not () : ScalarFunction ("!", 1)
{
}

double Not::evaluate(vector<double const*> const &args) const
{
    return *args[0] == 0;
}

bool Not::isDiscreteValued(vector<bool> const &) const
{
    return true;
}

string Not::deparse(vector<string> const &par) const
{
    return string("!") + par[0];
}

}}
