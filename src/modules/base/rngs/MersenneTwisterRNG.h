#ifndef MERSENNE_TWISTER_RNG_H_
#define MERSENNE_TWISTER_RNG_H_

#include <rng/RmathRNG.h>

namespace jags {
namespace base {

    class MersenneTwisterRNG : public RmathRNG
    {
	unsigned int dummy[625];
	unsigned int *mt; /* the array for the state vector */
	int mti;
	void fixupSeeds(bool init);
	void MT_sgenrand(unsigned int seed);
    public:
	MersenneTwisterRNG(unsigned int seed, NormKind norm_kind);
	void init(unsigned int seed);
	bool setState(std::vector<int> const &state);
	void getState(std::vector<int> &state) const;
	double uniform();
    };

}}

#endif /* MERSENNE_TWISTER_RNG_H_ */
