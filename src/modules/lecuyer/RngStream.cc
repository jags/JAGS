/*
 * Combined Multiple Recursive Generator MRG32k3a
 *
 * L'Ecuyer P (1999) Good parameters and implementations for combined
 * multiple recursive random number generators. Operations Research 47(1),
 * 159-164.
 * 
 * L'Ecuyer P, Simard R, Chen EJ, Kelton WD (2002) An object-oriented
 * random-number package with many long streams and substreams, Operations
 * Resarch 50(6), 1073-1075.

 * Adapted from the original C code by Pierre L'Ecuyer in the R package
 * rlecuyer 0.2.  The core of the code is copyright Pierre L'Ecuyer,
 * University of Montreal (lecuyer@iro.UMontreal.ca) 2001, and
 * released under the GPL license.
 *
 * The interface has been simplified and adapted to JAGS. The options
 * to produce antithetical and/or high- precision random variables
 * have been removed.  The code to advance the state of the RNG to the
 * start of the next stream or substream is in RngStreamFactory.cc
 */

#include "RngStream.h"
#include <module/ModuleError.h>

using std::vector;

#define norm  2.328306549295727688e-10 /* 2^-32 */
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0


static bool checkState (unsigned int state[6])
{
    for (int i = 0; i < 3; ++i) {
	if (state[i] >= m1) 
	    return false;
    }
    for (int i = 3; i < 6; ++i) {
	if (state[i] >= m2) 
	    return false;
    }
    if (state[0] == 0 && state[1] == 0 && state[2] == 0) {
	return false;
    }
    if (state[3] == 0 && state[4] == 0 && state[5] == 0) {
	return false;
    }

    return true;
}

static void fixState (unsigned int state[6])
{
    for (int i = 0; i < 3; ++i) {
	if (state[i] >= m1) 
	    state[i] = 0;
    }
    for (int i = 3; i < 6; ++i) {
	if (state[i] >= m2) 
	    state[i] = 0;
    }
    if (state[0] == 0 && state[1] == 0 && state[2] == 0) {
	state[0] = 1;
    }
    if (state[3] == 0 && state[4] == 0 && state[5] == 0) {
	state[3] = 1;
    }
}

namespace jags {
namespace lecuyer {

    RngStream::RngStream(unsigned int state[6])
	: RmathRNG("lecuyer::RngStream", KINDERMAN_RAMAGE)
    {
	if (!checkState(state)) {
	    throwLogicError("Invalid initial state in RngStream");
	}
	for (int i = 0; i < 6; ++i) {
	    Cg[i] = state[i];
	}
    }

    bool RngStream::setState(vector<int> const &state)
    {
	/* Pass from int to double via unsigned int */

	if (state.size() != 6)
	    return false;

	unsigned int cg[6];
	for (unsigned int i = 0; i < 6; ++i) {
	    cg[i] = static_cast<unsigned int>(state[i]);
	}
    
	if (!checkState(cg))
	    return false;
	
	for (int i = 0; i < 6; ++i) {
	    Cg[i] = cg[i];
	}
	return true;
    }

    void RngStream::getState(vector<int> &state) const
    {
	/* We pass from double to integer via unsigned integer
	   so as not to lose any information */
	   
	state.clear();

	unsigned int cg[6];
	for (int i = 0; i < 6; ++i) {
	    cg[i] = static_cast<unsigned int>(Cg[i]);
	}

	for (int i = 0; i < 6; ++i) {
	    state.push_back(static_cast<int>(cg[i]));
	}
    }

    double RngStream::uniform ()
    {
	/* Component 1 */
	double p1 = a12 * Cg[1] - a13n * Cg[0];
	long k = static_cast<long> (p1 / m1);
	p1 -= k * m1;
	if (p1 < 0.0)
	    p1 += m1;
	Cg[0] = Cg[1];
	Cg[1] = Cg[2];
	Cg[2] = p1;

	/* Component 2 */
	double p2 = a21 * Cg[5] - a23n * Cg[3];
	k = static_cast<long> (p2 / m2);
	p2 -= k * m2;
	if (p2 < 0.0)
	    p2 += m2;
	Cg[3] = Cg[4];
	Cg[4] = Cg[5];
	Cg[5] = p2;

	/* Combination */
	return ((p1 > p2) ? (p1 - p2) * norm : (p1 - p2 + m1) * norm);
    }

    void RngStream::init(unsigned int seed)
    {
	unsigned int state[6];
	seedToState(seed, state);
	for (int i = 0; i < 6; ++i) {
	    Cg[i] = state[i];
	}

    }

    void RngStream::seedToState(unsigned int seed, unsigned int state[6])
    {
        for (unsigned int i = 0; i < 6; i++) {
	    for(unsigned int j = 0; j < 50; j++) {
		seed = (69069 * seed + 1);
	    }
	    state[i] = seed;
        }
	fixState(state);
    }

}}
