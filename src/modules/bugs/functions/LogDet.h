#ifndef LOG_DET_H_
#define LOG_DET_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Log determinant 
     * LogDet calculates the log determinant of a square matrix.  The
     * function assumes that the matrix is symmetric positive definite.
     * but currently does not test this.
     * <pre>
     * y <- logdet(x)
     * y <- log|x| for x an n x n symmetric positive definite matrix
     * </pre>
     */
    class LogDet : public ArrayFunction
    {
    public:
	LogDet ();
	void evaluate(double *x, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims) 
	    const;
	bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) const;
	std::vector<unsigned long>
	    dim(std::vector<std::vector<unsigned long> > const &dims,
		std::vector<double const *> const &values) const;
    };

}}

#endif /* LOG_DET_H_ */
