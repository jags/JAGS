#ifndef ARRAY_LOG_DENSITY_H_
#define ARRAY_LOG_DENSITY_H_

#include <function/ArrayFunction.h>

namespace jags {
    
    class ArrayDist;

    /**
     * @short Log density function for an array-valued Distribution
     */
    class ArrayLogDensity : public ArrayFunction
    {
	ArrayDist const *_dist;
    public:
	ArrayLogDensity(ArrayDist const *dist);
	std::vector<unsigned long> dim(
	    std::vector<std::vector<unsigned long> > const &dims,
	    std::vector<double const *> const &values) const;
	bool checkParameterDim(
	    std::vector<std::vector<unsigned long> > const &dims) const;
	bool checkParameterValue(
	    std::vector<double const *> const &args,
	    std::vector<std::vector<unsigned long> > const &dims) const;
	void evaluate(
	    double *value,
	    std::vector <double const *> const &args,
	    std::vector<std::vector<unsigned long> > const &dims) const;
    };

}

#endif /* ARRAY_LOG_DENSITY_H_ */
