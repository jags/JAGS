#ifndef FUNC_SD_H_
#define FUNC_SD_H_

#include <function/ScalarVectorFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Standard deviation function
     * SD calculates the standard deviation of the elements of an array.
     * The array must have at least 2 elements.
     * @see mean
     * <pre>
     * y <- sd(x[])
     * </pre>
     */
    class SD : public ScalarVectorFunction
    {
    public:
	SD ();
	double scalarEval (std::vector<double const *> const &args,
			   std::vector<unsigned long> const &lengths) const;
	bool checkParameterLength (std::vector<unsigned long> const &len) const;
    };

}}

#endif /* FUNC_SD_H_ */
