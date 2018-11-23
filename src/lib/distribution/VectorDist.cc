#include <config.h>
#include <distribution/VectorDist.h>
#include <util/nainf.h>

using std::string;
using std::vector;

namespace jags {

VectorDist::VectorDist(string const &name, unsigned int npar)
  : Distribution(name, npar)
{
}

    double VectorDist::KL(vector<double const *> const &par1,
			  vector<double const *> const &par2,
			  vector<unsigned long> const &lengths,
			  RNG *rng, unsigned int nrep) const
    {
	double div = 0;

	unsigned long N = length(lengths);
	vector<double> v(N);
	for (unsigned int r = 0; r < nrep; ++r) {
	    randomSample(&v[0], par1, lengths, rng);
	    div += logDensity(&v[0], PDF_FULL, par1, lengths);
	    div -= logDensity(&v[0], PDF_FULL, par2, lengths);
	}
	return div / nrep;
    }

    double VectorDist::KL(std::vector<double const *> const &,
			  std::vector<double const *> const &,
			  std::vector<unsigned long> const &) const
    {
	return JAGS_NA;
    }
    
} //namespace jags
