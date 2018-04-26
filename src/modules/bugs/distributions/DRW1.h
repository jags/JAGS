#ifndef DRW1_H_
#define DRW1_H_

#include <distribution/VectorDist.h>

namespace jags {
namespace bugs {

/**
 * @short First order random walk
 *
 * <pre>
 * y[] ~ drw1(tau, x[])
 * </pre>
 */
class DRW1 : public VectorDist {
public:
    DRW1();
    double logDensity(double const *x, unsigned long length, PDFType type,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned long> const &lengths,
		      double const *lower, double const *upper) const;
    void randomSample(double *x, unsigned long length,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned long> const &lengths,
		      double const *lower, double const *upper, RNG *rng) const;
    unsigned long length(std::vector<unsigned long> const &lengths) const;
    bool checkParameterLength(std::vector<unsigned long> const &lengths) const;
    bool checkParameterValue(std::vector<double const *> const &parameters,
			     std::vector<unsigned long> const &lengths)
	const;
    void support(double *lower, double *upper, unsigned long length,
		 std::vector<double const *> const &parameters,
		 std::vector<unsigned long> const &lengths) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    unsigned long df(std::vector<unsigned long> const &lengths) const;
};

}}

#endif /* DRW1_H_ */
