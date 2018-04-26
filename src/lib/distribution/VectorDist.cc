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

unsigned long VectorDist::df(vector<unsigned long> const &par) const
{
    return length(par);
}

    double VectorDist::KL(vector<double const *> const &par1,
			  vector<double const *> const &par2,
			  vector<unsigned long> const &lengths,
			  double const *lower, double const *upper,
			  RNG *rng, unsigned int nrep) const
    {
	double div = 0;

	unsigned long N = length(lengths);
	vector<double> v(N);
	for (unsigned int r = 0; r < nrep; ++r) {
	    randomSample(&v[0], N, par1, lengths, lower, upper, rng);
	    div += logDensity(&v[0], N, PDF_FULL, par1, lengths, lower, upper);
	    div -= logDensity(&v[0], N, PDF_FULL, par2, lengths, lower, upper);
	}
	return div / nrep;
    }

    double VectorDist::KL(std::vector<double const *> const &par1,
			  std::vector<double const *> const &par2,
			  std::vector<unsigned long> const &lengths) const
    {
	return JAGS_NA;
    }
    
} //namespace jags
