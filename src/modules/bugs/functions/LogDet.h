#ifndef LOG_DET_H_
#define LOG_DET_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Log determinant 
     * LogDet calculates the log determinant of a square matrix.  The
     * function assumes that the matrix is symmetric positive definite.
     * but currently does not test this. The matrix is assumed to be
     * symmetric as the function only reads the lower triangle.
     * <pre>
     * y <- logdet(x)
     * y = log|x| for x an n x n symmetric positive definite matrix
     * </pre>
     */
    class LogDet : public ArrayFunction
    {
    public:
	LogDet ();
	void evaluate(double *x, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned long> > const &dims) 
	    const override;
	bool checkParameterDim(std::vector<std::vector<unsigned long>>
			       const &dims) const override;
	std::vector<unsigned long>
	    dim(std::vector<std::vector<unsigned long> > const &dims,
		std::vector<double const *> const &values) const override;
	bool hasGradient(unsigned long i) const override;
	void gradient(double *grad, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned long>> const &dims,
		      unsigned long i) const override;
    };

}}

#endif /* LOG_DET_H_ */
