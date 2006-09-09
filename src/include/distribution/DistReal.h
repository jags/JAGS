#ifndef DIST_REAL_H_
#define DIST_REAL_H_

#include <distribution/DistScalar.h>

/**
 * Abstract class for real valued distributions.
 *
 * @short Real valued distributions
 */
class DistReal : public DistScalar
{
 public:
  /**
   * Constructor
   *
   * @param name BUGS language name of distribution
   *
   * @param npar Number of parameters, excluding upper and lower bound
   *
   * @param support Support of distribution
   *
   * @param boundable Logical flag indicating whether the distribution
   * can be bounded using the T(,) construct.
   */
  DistReal(std::string const &name, unsigned int npar,
	   Support support, bool boundable);
  double logLikelihood(SArray const &x,
		       std::vector<SArray const *> const &parameters) const;
  void randomSample(SArray &x, std::vector<SArray const *> const &parameters,
		    RNG *rng) const;
};

#endif /* DIST_REAL_H_ */
