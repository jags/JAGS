#ifndef FUNC_COMBINE_H_
#define FUNC_COMBINE_H_

#include <function/VectorFunction.h>

namespace jags {
    namespace bugs {

	/**
	 * @short Combine values into a vector
	 *
	 * Concatenates all arguments into a single vector
	 * <pre>
	 * y <- c(x[],y, z[,])
	 * </pre>
	 */
	class Combine : public VectorFunction
	{
	  public:
	    Combine ();
	    void evaluate(double *value, 
			  std::vector <double const *> const &args,
			  std::vector <unsigned long> const &lengths) const;
	    unsigned long 
		length(std::vector<unsigned long> const &parlengths,
		       std::vector<double const *> const &parvalues) const;
	    bool isDiscreteValued(std::vector<bool> const &mask) const;
	    bool isAdditive(std::vector<bool> const &mask,
			    std::vector<bool> const &isfixed) const;
	    bool isScale(std::vector<bool> const &mask,
			 std::vector<bool> const &isfixed) const;
	    bool isLinear(std::vector<bool> const &mask,
			  std::vector<bool> const &isfixed) const;
	    bool checkParameterLength(std::vector<unsigned long> const &args)
		const;
	};

    }
}

#endif /* FUNC_COMBINE_H_ */
