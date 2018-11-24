#ifndef DSCALED_WISHART_H_
#define DSCALED_WISHART_H_

#include <distribution/ArrayDist.h>

namespace jags {
namespace glm {

/**
 * <pre>
 * x[] ~ dwish(S[], k)
 * </pre>
 * @short Huang and Wand scaled Wishart distribution
 */
class DScaledWishart : public ArrayDist {
public:
  DScaledWishart();

  double logDensity(double const *x, PDFType type,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims) const;
  void randomSample(double *x,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims,
		    RNG *rng) const;
  static void sampleWishart(double *x,
			    double const *scale, unsigned long nrow,
			    double k, RNG *rng);
  /**
   * Checks that S is a vector and k is a scalar
   */
  bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) 
      const;
  /**
   * Checks that S and k are both positive
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
			   std::vector<std::vector<unsigned long> > const &dims)
      const;
  std::vector<unsigned long> 
      dim(std::vector<std::vector<unsigned long> > const &dims) const;
  void support(double *lower, double *upper,
	       std::vector<double const *> const &parameters,
	       std::vector<std::vector<unsigned long> > const &dims) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
  bool fullRank() const;
};

}}

#endif /* DSCALED_WISHART_H_ */
