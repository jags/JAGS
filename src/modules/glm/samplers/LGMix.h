#ifndef LG_MIX_H_
#define LG_MIX_H_

#include <vector>

namespace jags {

struct RNG;

namespace glm {

    /**
     * @short Finite mixture representation of negative log-gamma distribution
     *
     * Represents a negative log gamma distribution with integer shape
     * n > 0 and scale 1 as a finite mixture of normals. This is used
     * by the AuxMixPoisson and AuxMixBinomial Outcome classes.
     */
    class LGMix {
	double _n;
	long _r;
	unsigned int _ncomp;	
	double _weights[10];
	double _means[10];
	double _variances[10];
	void updateShapeExact(int n);
	void updateShapeApprox(double n);
	void updateShape(double n);
    public:
	LGMix(double n);
	/**
	 * Updates the mixture representation
	 *
	 * @param z Value sampled from negative log gamma distribution
	 *          with shape parameter n.
	 * @param n Value of n, which may be different from value used in
	 *          constructor. 
	 * @param rng Random number generator
	 */
	void update(double z, double n, RNG *rng);
	/**
	 * Returns the mean of the current normal component
	 */
	double mean() const;
	/**
	 * Returns the precision of the current normal component
	 */
	double precision() const;
	/**
	 * Gets the parameters of the mixture distribution. 
	 * 
	 * This is not required by the GLM machinery, but it is useful
	 * to expose the parameters of the mixture distribution for
	 * testing purposes.
	 */
	void getParameters(std::vector<double> &weights,
			   std::vector<double> &means,
			   std::vector<double> &variances);
    };
    
}}

#endif /* LG_MIX_H */
