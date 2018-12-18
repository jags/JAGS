#ifndef SCALAR_LOG_DENSITY_H_
#define SCALAR_LOG_DENSITY_H_

#include <function/ScalarFunction.h>

namespace jags {

    class ScalarDist;

    /**
     * @short Log density function for a scalar-valued Distribution
     */
    class ScalarLogDensity : public ScalarFunction
    {
	ScalarDist const *_dist;
    public:
	ScalarLogDensity(ScalarDist const *dist);
	bool checkParameterValue(std::vector<double const *> const &args) const;
	double evaluate(std::vector <double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
    };

}

#endif /* SCALAR_LOG_DENSITY_H_ */
