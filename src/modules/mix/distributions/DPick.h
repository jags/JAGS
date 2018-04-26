#ifndef DPICK_H_
#define DPICK_H_

#include <distribution/ScalarDist.h>

namespace jags {
    namespace mix {

	class DPick : public ScalarDist {
	  public:
	    DPick();

	    double logDensity(double const x, PDFType type,
			      std::vector<double const *> const &parameters,
			      double const *lower, double const *upper) const;
	    double randomSample(std::vector<double const *> const &parameters,
				double const *lower, double const *upper,
				RNG *rng) const;
	    double typicalValue(std::vector<double const *> const &parameters,
				double const *lower, double const *upper) const;
	    void support(double *lower, double *upper, 
			 std::vector<double const *> const &parameters) const;
	    bool isSupportFixed(std::vector<bool> const &fixmask) const;
	    bool checkParameterValue(std::vector<double const *> const &par)
		const;
	    double typicalValue(std::vector<double const *> const &par) const;
	    unsigned long length(std::vector<unsigned long> const &lengths) const;
	    bool isDiscreteValued(std::vector<bool> const &mask) const;
	};
	
    }
}


#endif /* DPICK_H_ */
