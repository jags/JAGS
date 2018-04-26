#ifndef FUNC_DINTERVAL_H_
#define FUNC_DINTERVAL_H_

#include <function/ScalarVectorFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Function version of dinterval distribution
     */
    class DIntervalFunc: public ScalarVectorFunction
    {
    public:
	DIntervalFunc ();
	double scalarEval (std::vector <double const *> const &args,
			   std::vector<unsigned long> const &lengths) const;
	bool checkParameterLength (std::vector<unsigned long> const &args) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };

}}

#endif /* FUNC_DINTERVAL_H_ */
