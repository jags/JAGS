#ifndef RNG_STREAM_H_
#define RNG_STREAM_H_

#include <rng/RmathRNG.h>

#define m1    4294967087.0
#define m2    4294944443.0

namespace jags {
namespace lecuyer {
    
    /**
     * @short MRG32k3a combined multiple recursive generator
     */
    class RngStream : public RmathRNG {
	double Cg[6];
    public:
	/**
	 * Constructor for RngStream random number generator
	 *
	 * @param state Initial state for the RNG. The first three
	 * elements must be less than m1, and the second three must be
	 * less than m2. The first three elements may not be all zero
	 * and the second elements may not be all zero.
	 */
	RngStream(unsigned int state[6]);
	void init(unsigned int seed);
	bool setState(std::vector<int> const &state);
	void getState(std::vector<int> &state) const;
	double uniform();
	/**
	 * Generates a state vector from a random seed
	 */
	static void seedToState(unsigned int seed, unsigned int state[6]);
    };
 
}}

#endif /* RNGSTREAM_H_ */
 

