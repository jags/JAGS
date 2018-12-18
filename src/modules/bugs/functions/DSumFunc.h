#ifndef DSUM_FUNC_H_
#define DSUM_FUNC_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace bugs {
    
    /**
     * @short Sum of two discrete random variables
     */
    class DSumFunc : public ArrayFunction {
    public:
	DSumFunc();
	void evaluate(double *x,
		      std::vector <double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims) 
	    const;
	bool isDifferentiable(unsigned long i) const;
	void gradient(double *x, std::vector <double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims,
		      unsigned long i) const;
	bool checkParameterDim(std::vector<std::vector<unsigned long> > const 
			       &dims) const;
	std::vector<unsigned long> 
	    dim(std::vector <std::vector<unsigned long> > const &dims,
		std::vector<double const *> const &values) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
	bool isLinear(std::vector<bool> const &mask, 
		      std::vector<bool> const &fixed) const;
	bool isScale(std::vector<bool> const &mask, 
		     std::vector<bool> const &fixed) const;
    };
    
}}

#endif /* DSUM_FUNC_H_ */
