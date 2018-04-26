#ifndef FUNC_MAX_H_
#define FUNC_MAX_H_

#include <function/ScalarVectorFunction.h>

namespace jags {
namespace bugs {
    
    /**
     * @short Maximum element of a vector
     * @see Min
     */
    class Max: public ScalarVectorFunction
    {
    public:
	Max ();
	double scalarEval(std::vector <double const *> const &args,
			  std::vector<unsigned long> const &lengths) const;
	/** Returns true if all arguments are true */
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };

}}

#endif /* FUNC_MAX_H_ */
