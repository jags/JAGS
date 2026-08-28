#ifndef INTERP_LIN_2D_H_
#define INTERP_LIN_2D_H_

#include <function/ArrayFunction.h>

namespace jags {
    namespace bugs {
	
	class InterpLin2D : public ArrayFunction
	{
	public:
	    InterpLin2D();
	    void evaluate(double *value, std::vector<double const *> const &args,
			  std::vector<std::vector<unsigned long>> const &dims) 
		const override;
	    std::vector<unsigned long> 
	    dim(std::vector<std::vector<unsigned long>> const &dims,
		std::vector<double const *> const &values) const override;
	    bool checkParameterDim(std::vector<std::vector<unsigned long>>
				   const &dims) const override;
	    bool checkParameterValue(std::vector<double const *> const &args,
				     std::vector<std::vector<unsigned long>> const &dims)
		const override;
	};

}}

#endif /* INTERP_LIN_2D_H_ */
