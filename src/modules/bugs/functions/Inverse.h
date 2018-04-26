#ifndef INVERSE_H_
#define INVERSE_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Inverts a square symmetric matrix
     * <pre>
     * y[,] <- inverse(x[,])
     * </pre>
     */
    class Inverse: public ArrayFunction
    {
    public:
	Inverse ();
	void evaluate (double *value, std::vector <double const *> const &args,
		       std::vector<std::vector<unsigned long> > const &dims) 
	    const;
	std::vector<unsigned long> 
	    dim(std::vector<std::vector<unsigned long> > const &args,
		std::vector<double const *> const &values) const;
	bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) const;
	bool checkParameterValue(std::vector<double const *> const &args,
				 std::vector<std::vector<unsigned long> > const &dims) const;
    };

}}

#endif /* INVERSE_H_ */
