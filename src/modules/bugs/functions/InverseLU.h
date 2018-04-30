#ifndef INVERSE_LU_H_
#define INVERSE_LU_H_

#include <function/ArrayFunction.h>

namespace jags {
    namespace bugs {

	/**
	 * @short Inverts a square matrix using the LU decomposition
	 * <pre>
	 * y[,] <- inverse.lu(x[,])
	 * </pre>
	 */
	class InverseLU: public ArrayFunction
	{
	public:
	    InverseLU ();
	    void evaluate (double *value, std::vector <double const *> const &args,
			   std::vector<std::vector<unsigned long> > const &dims) const;
	    std::vector<unsigned long> 	dim(std::vector<std::vector<unsigned long> > const &args,
					    std::vector<double const *> const &values) const;
	    bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) const;
	};
	
    }
}

#endif /* INVERSE_LU_H_ */
