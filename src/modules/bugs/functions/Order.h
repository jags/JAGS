#ifndef FUNC_ORDER_H_
#define FUNC_ORDER_H_

#include <function/VectorFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Order function 
     *
     * Returns a permutation that sorts the elements of a vector in
     * ascending order. This is the inverse of the permutation
     * returned by the Rank function.
     *
     * <pre> 
     * y <- sort(x[])
     * </pre>
     *
     * @see Rank
     */
    class Order : public VectorFunction
    {
    public:
	Order ();
	void evaluate(double *value, 
		      std::vector<double const *> const &args,
		      std::vector<unsigned long> const &lengths) const;
	unsigned long length(std::vector<unsigned long> const &parlengths,
			    std::vector<double const *> const &parvalues) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };
    
}}

#endif /* FUNC_ORDER_H_ */
