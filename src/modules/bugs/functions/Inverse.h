#ifndef INVERSE_H_
#define INVERSE_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Inverts a symmetric positive-definite matrix using the Cholesky
     * decomposition
     * <pre>
     * y <- inverse.chol(x[,])
     * y <- inverse(x)
     * </pre>
     */
    class Inverse: public ArrayFunction
    {
    public:
	Inverse ();
	void evaluate(double *value, std::vector <double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims) 
	    const override;
	std::vector<unsigned long> 
	    dim(std::vector<std::vector<unsigned long> > const &args,
		std::vector<double const *> const &values) const override;
	bool checkParameterDim(std::vector<std::vector<unsigned long>> const &dims) const override;
	std::string alias() const override;
	bool hasGradient(unsigned long i) const override;
	void gradient(double *grad, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims,
		      unsigned long i) const override;
    };

}}

#endif /* INVERSE_H_ */
