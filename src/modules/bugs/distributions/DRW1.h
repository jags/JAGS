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
    double logDensity(double const *x, unsigned int length, PDFType type,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned int> const &lengths,
		      double const *lower, double const *upper) const;
    void randomSample(double *x, unsigned int length,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned int> const &lengths,
		      double const *lower, double const *upper, RNG *rng) const;
    void typicalValue(double *x, unsigned int length,
		      std::vector<double const *> const &parameter,
		      std::vector<unsigned int> const &lengths,
		      double const *lower, double const *upper) const;
    unsigned int length(std::vector<unsigned int> const &lengths) const;
    bool checkParameterLength(std::vector<unsigned int> const &lengths) const;
    bool checkParameterValue(std::vector<double const *> const &parameters,
			     std::vector<unsigned int> const &lengths)
	const;
    void support(double *lower, double *upper, unsigned int length,
		 std::vector<double const *> const &parameters,
		 std::vector<unsigned int> const &lengths) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    unsigned int df(std::vector<unsigned int> const &lengths) const;
};

}}

#endif /* DRW1_H_ */
