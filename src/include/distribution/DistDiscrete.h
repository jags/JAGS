#ifndef DIST_DISCRETE_H_
#define DIST_DISCRETE_H_

#include <distribution/DistScalar.h>

/**
 * Abstract class for integer valued distributions
 *
 * @see DistReal
 * @short Integer valued distributions
 */
class DistDiscrete : public DistScalar
{
public:
  /**
   * Constructor for discrete distribution
   *
   * @param name BUGS language name of the distribution
   * @param npar Number of parameters, excluding upper and lower bounds
   * @param canbound Logical flag indicating whether distribution is boundable.
   */
  DistDiscrete(std::string const &name, unsigned int npar, 
	       Support support, bool canbound);
  double logLikelihood(SArray const &x,
		       std::vector<SArray const *> const &parameters) const;
  void randomSample(SArray &x,
		    std::vector<SArray const *> const &parameters,
                    RNG *rng) const;
};

#endif /* DIST_DISCRETE_H_ */
