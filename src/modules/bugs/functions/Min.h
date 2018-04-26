#ifndef FUNC_MIN_H_
#define FUNC_MIN_H_

#include <function/ScalarVectorFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Minimum element of a vector
     * @see Max
     */
    class Min: public ScalarVectorFunction
    {
    public:
	Min ();
	double scalarEval (std::vector <double const *> const &args,
			   std::vector<unsigned long> const &lengths) const;
	/** Returns true if all arguments are true */
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };

}}

#endif /* FUNC_MIN_H_ */
