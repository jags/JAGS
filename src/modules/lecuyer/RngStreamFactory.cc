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

 * Adapted from the original code by Pierre L'Ecuyer in the R package
 * rlecuyer 0.2.  The core of the code is copyright Pierre L'Ecuyer,
 * University of Montreal (lecuyer@iro.UMontreal.ca) 2001, and
 * released under the GPL license.
 *
 * The core algorithm is untouched, but the interface has been
 * simplified and adapted to JAGS. The options to produce antithetical
 * and/or high- precision random variables have been removed.  The
 * code to advance the state of the RNG to the start of the next stream
 * or substream is in RngStreamFactory.cc 
 */

#include "RngStreamFactory.h"
#include "RngStream.h"

#include <ctime>

using std::vector;
using std::time;
using std::string;

#define two17   131072.0               /* 2^17 */
#define two53   9007199254740992.0     /* 2^53 */

static const double A1p76[3][3] = {
          {      82758667.0, 1871391091.0, 4127413238.0 }, 
          {    3672831523.0,   69195019.0, 1871391091.0 }, 
          {    3672091415.0, 3528743235.0,   69195019.0 }
          };

static const double A2p76[3][3] = {
          {    1511326704.0, 3759209742.0, 1610795712.0 }, 
          {    4292754251.0, 1511326704.0, 3889917532.0 }, 
          {    3859662829.0, 4292754251.0, 3708466080.0 }
          };

static const double A1p127[3][3] = {
          {    2427906178.0, 3580155704.0,  949770784.0 }, 
          {     226153695.0, 1230515664.0, 3580155704.0 },
          {    1988835001.0,  986791581.0, 1230515664.0 }
          };

static const double A2p127[3][3] = {
          {    1464411153.0,  277697599.0, 1610723613.0 },
          {      32183930.0, 1464411153.0, 1022607788.0 },
          {    2824425944.0,   32183930.0, 2093834863.0 }
          };


static double MultModM (double a, double s, double c, double m)
   /* Compute (a*s + c) % m. m must be < 2^35.  Works also for s, c < 0 */
{
   double v;
   long a1;
   v = a * s + c;
   if ((v >= two53) || (v <= -two53)) {
      a1 = static_cast<long>(a / two17);
      a -= a1 * two17;
      v = a1 * s;
      a1 = static_cast<long>(v / m);
      v -= a1 * m;
      v = v * two17 + a * s + c;
   }
   a1 = static_cast<long>(v / m);
   if ((v -= a1 * m) < 0.0)
      return v += m;
   else
      return v;
}


static void MatVecModM (const double A[3][3], double s[3], double m)
   /* Returns s = A*s % m.  Assumes that -m < s[i] < m. */
{
   double x[3] = {0,0,0};
   for (int i = 0; i < 3; ++i) {
       for (int j = 0; j < 3; ++j) {
	   x[i] = MultModM (A[i][j], s[j], x[i], m);
       }
   }
   for (int i = 0; i < 3; ++i)
       s[i] = x[i];
}

namespace jags {
namespace lecuyer {

    RngStreamFactory::RngStreamFactory()
    {
        unsigned int seed = static_cast<unsigned int>(time(NULL));
	setSeed(seed);
    }
    
    RngStreamFactory::~RngStreamFactory()
    {
	for (unsigned int i = 0; i < _rngvec.size(); ++i) {
	    delete _rngvec[i];
	}
    }

    void RngStreamFactory::setSeed(unsigned int seed)
    {
	unsigned int state[6];
	RngStream::seedToState(seed, state);

	for (int i = 0; i < 6; ++i) {
	    Bg[i] = Ig[i] = state[i];
	}
    }

    void RngStreamFactory::nextStream ()
    {
	MatVecModM (A1p127, Ig, m1);
	MatVecModM (A2p127, &Ig[3], m2);
	for (int i = 0; i < 6; ++i) {
	    Bg[i] = Ig[i];
	}
    }

    void RngStreamFactory::nextSubstream ()
    {
	MatVecModM (A1p76, Bg, m1);
	MatVecModM (A2p76, &Bg[3], m2);
    }

    vector<RNG *> RngStreamFactory::makeRNGs(unsigned int n)
    {
	vector<RNG *> ans;
	for (unsigned int i = 0; i < n; ++i) {

	    unsigned int state[6];
	    for (int j = 0; j < 6; ++j) {
		state[j] = static_cast<unsigned int>(Bg[j]);
	    }
	    RNG *rng = new RngStream(state);
	    nextSubstream();
	    
	    _rngvec.push_back(rng);
	    ans.push_back(rng);
	}
	nextStream();
	return ans;
    }

    RNG * RngStreamFactory::makeRNG(string const &name)
    {
	if (name == "lecuyer::RngStream") {

	    unsigned int state[6];
	    for (int j = 0; j < 6; ++j) {
		state[j] = static_cast<unsigned int>(Bg[j]);
	    }
	    RNG *rng = new RngStream(state);
	    nextStream();

	    _rngvec.push_back(rng);
	    return rng;
	}
	else {
	    return 0;
	}
	
    }

    string RngStreamFactory::name() const
    {
	return "lecuyer::RngStream";
    }

}}
