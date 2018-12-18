#ifndef TRANSPOSE_H_
#define TRANSPOSE_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace bugs {

    class Transpose : public ArrayFunction
    {
    public:
	Transpose();
	void evaluate(double *value, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims) 
	    const;
	std::vector<unsigned long> 
	    dim(std::vector<std::vector<unsigned long> > const &dims,
		std::vector<double const *> const &values) const;
	bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) const;
	bool isAdditive(std::vector<bool> const &mask,
			std::vector<bool> const &isfixed) const;
	bool isScale(std::vector<bool> const &mask,
		     std::vector<bool> const &isfixed) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
	bool isDifferentiable(unsigned long i) const;
	void gradient(double *value, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims,
		      unsigned long i) const;
	
    };

}}

#endif /* TRANSPOSE_H_ */
