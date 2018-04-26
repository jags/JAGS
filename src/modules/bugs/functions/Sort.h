#ifndef FUNC_SORT_H_
#define FUNC_SORT_H_

#include <function/VectorFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Sort function
     * Sorts the elements of a vector in ascending order
     * <pre>
     * y <- sort(x[])
     * </pre>
     */
    class Sort : public VectorFunction
    {
    public:
	Sort ();
	void evaluate(double *value, std::vector <double const *> const &args,
		      std::vector <unsigned long> const &lengths) const;
	unsigned long length(std::vector<unsigned long> const &parlengths,
			    std::vector<double const *> const &parvalues) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };

}}

#endif /* FUNC_SORT_H_ */
