#ifndef VECTOR_LOG_DENSITY_H_
#define VECTOR_LOG_DENSITY_H_

#include <function/VectorFunction.h>

namespace jags {
    
    class VectorDist;

    /**
     * @short Log density function for a vector-valued Distribution
     */
    class VectorLogDensity : public VectorFunction
    {
	VectorDist const *_dist;
    public:
	VectorLogDensity(VectorDist const *dist);
	unsigned long length(std::vector<unsigned long> const &lengths,
			    std::vector<double const *> const &values) const;
	bool checkParameterLength(std::vector<unsigned long> const &lens) const;
	bool checkParameterValue(std::vector<double const *> const &args,
				 std::vector<unsigned long> const &lens) const;
	void evaluate(double *value,
		      std::vector <double const *> const &args,
		      std::vector<unsigned long> const &lens) const;
	bool isDifferentiable(unsigned long i) const;
    };

}

#endif /* VECTOR_LOG_DENSITY_H_ */
