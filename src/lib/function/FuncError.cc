#include <config.h>
#include <function/FuncError.h>
#include <function/Function.h>

using std::string;

namespace jags {

FuncError::FuncError(Function const *efunc, string const &emsg)
    : runtime_error(emsg + " in function '" + efunc->name() + "'")
{
}

} //namespace jags
