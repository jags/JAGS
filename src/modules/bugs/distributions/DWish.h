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

  double logDensity(double const *x, unsigned long length, PDFType type,
		       std::vector<double const *> const &parameters,
		       std::vector<std::vector<unsigned long> > const &dims,
		       double const *lower, double const *upper) const;
  void randomSample(double *x, unsigned long length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims,
		    double const *lower, double const *upper, RNG *rng) const;
  //FIXME: Can we retire this?
  static void randomSample(double *x, unsigned long length,
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
  void support(double *lower, double *upper, unsigned long length,
	       std::vector<double const *> const &parameters,
	       std::vector<std::vector<unsigned long> > const &dims) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
  unsigned long df(std::vector<std::vector<unsigned long> > const &dims) const;
};

}}

#endif /* DWISH_H_ */
