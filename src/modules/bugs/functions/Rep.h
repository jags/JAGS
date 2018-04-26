#ifndef REP_H_
#define REP_H_

#include <function/VectorFunction.h>

namespace jags {
    namespace bugs {

    /**
     * @short Replicate elements of a vector
     *
     * <pre>
     * y <- rep(x, times)
     * </pre>
     *
     */
    class Rep : public VectorFunction
    {
    public:
	Rep();
	void evaluate(double *value, 
		      std::vector <double const *> const &values,
		      std::vector <unsigned long> const &lengths) const;
	unsigned long length(std::vector <unsigned long> const &lengths,
			    std::vector <double const *> const &values) const;
	bool checkParameterLength(std::vector<unsigned long> const &args) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
	bool checkParameterDiscrete(std::vector<bool> const &mask) const;
	bool checkParameterFixed(std::vector<bool> const &mask) const;
	bool checkParameterValue(std::vector<double const *> const &args,
				 std::vector<unsigned long> const &lens) const;
    };

    } /* namespace bugs */
} /* namespace jags */

#endif /* REP_H_ */
