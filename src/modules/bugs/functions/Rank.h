#ifndef FUNC_RANK_H_
#define FUNC_RANK_H_

#include <function/VectorFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Rank function
     *
     * Calculates the sample ranks of the values in a vector. The ranks
     * define a permutation that is the inverse of the one returned
     * by the Order function.
     *
     * <pre>
     * y <- rank(x[])
     * </pre>
     *
     * @see Order
     */
    class Rank : public VectorFunction
    {
    public:
	Rank ();
	void evaluate(double *value, 
		      std::vector<double const *> const &args,
		      std::vector<unsigned long> const &lengths) const;
	unsigned long length(std::vector<unsigned long> const &parlengths,
			    std::vector<double const *> const &parvalues) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };
    
}}

#endif /* FUNC_RANK_H_ */
