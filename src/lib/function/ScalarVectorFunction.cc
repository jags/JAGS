#include <config.h>
#include <function/ScalarVectorFunction.h>

using std::vector;
using std::string;

namespace jags {

    ScalarVectorFunction::ScalarVectorFunction (string const &name, 
						unsigned long npar)
	: VectorFunction(name, npar)
    {
    }

    void 
    ScalarVectorFunction::evaluate(double *value, 
				   vector <double const *> const &args,
				   vector <unsigned long> const &lengths) const
    {
	*value = scalarEval(args, lengths);
    }

    unsigned long 
    ScalarVectorFunction::length(vector <unsigned long> const &lengths,
				 vector <double const *> const &values) const
    {
	return 1;
    }

}
