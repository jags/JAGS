#ifndef FUNC_INPROD_H_
#define FUNC_INPROD_H_

#include <function/ScalarVectorFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Inner product of two vectors
     * The InProd function returns the sum of the product of corresponding
     * elements from two arrays with the same dimension
     * <pre>
     * y <- inprod(a[],b[])
     * </pre>
     */
    class InProd : public ScalarVectorFunction
    {
    public:
	InProd ();
	double scalarEval(std::vector <double const *> const &args,
			  std::vector<unsigned long> const &lengths) const;
	bool checkParameterLength (std::vector<unsigned long> const &len) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
	bool isScale(std::vector<bool> const &mask,
		     std::vector<bool> const &fix) const;
	bool isScaleMix(std::vector<bool> const &mask,
			std::vector<bool> const &fix) const;
    };

}}

#endif /* FUNC_INPROD_H_ */
