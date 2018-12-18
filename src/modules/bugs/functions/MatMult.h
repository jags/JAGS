#ifndef MAT_MULT_H_
#define MAT_MULT_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace bugs {

    class MatMult : public ArrayFunction
    {
    public:
	MatMult();
	void evaluate(double *value, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims) 
	    const;
	bool isDifferentiable(unsigned long i) const;
	void gradient(double *grad, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims,
		      unsigned long i) const;
	std::vector<unsigned long> 
	    dim(std::vector<std::vector<unsigned long> > const &dims,
		std::vector<double const *> const &values) const;
	bool checkParameterDim(std::vector <std::vector<unsigned long> > const &dims) const;
	bool isScale(std::vector<bool> const &mask, 
		     std::vector<bool> const &fix) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };

}}

#endif /* MAT_MULT_H_ */
