#ifndef DWISH_H_
#define DWISH_H_

#include <distribution/ArrayDist.h>

namespace jags {
namespace bugs {

/**
 * <pre>
 * x[] ~ dwish(R[,], k)
 * </pre>
 * @short Wishart distribution
 */
class DWish : public ArrayDist {
public:
  DWish();

  double logDensity(double const *x, PDFType type,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims) const;
  void randomSample(double *x,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims,
		    RNG *rng) const;
  //FIXME: Can we retire this?
  static void randomSample(double *x,
                           double const *R, double k, unsigned long nrow,
                           RNG *rng);
  /**
   * Checks that R is a square matrix and k is a scalar
   */
  bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) 
      const;
  /**
   * Checks that R is symmetric and k >= nrow(R). There is
   * currently no check that R is positive definite
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

#endif /* DWISH_H_ */
