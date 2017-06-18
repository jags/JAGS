#ifndef FUNC_SUM_H_
#define FUNC_SUM_H_

#include <function/ScalarVectorFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Sum of an array
     * Sum calculates the sum of the elements of an array
     * @see Mean
     * <pre>
     * y <- sum(x1[], x2[], ...)
     * </pre>
     */
    class Sum : public ScalarVectorFunction
    {
    public:
	Sum ();
	double scalarEval(std::vector <double const *> const &args,
			  std::vector<unsigned int> const &lengths) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
	bool isAdditive(std::vector<bool> const &mask,
			std::vector<bool> const &isfixed) const;
	bool isScale(std::vector<bool> const &mask,
		     std::vector<bool> const &isfixed) const;
	bool isLinear(std::vector<bool> const &mask,
		      std::vector<bool> const &isfixed) const;
    };
    
}}

#endif /* FUNC_SUM_H_ */
