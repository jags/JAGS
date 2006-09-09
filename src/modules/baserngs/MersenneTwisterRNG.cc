#include <ctime>
#include "MersenneTwisterRNG.h"

/* ===================  Mersenne Twister ========================== */
/* From http://www.math.keio.ac.jp/~matumoto/emt.html */

/* A C-program for MT19937: Real number version([0,1)-interval)
   (1999/10/28)
     genrand() generates one pseudorandom real number (double)
   which is uniformly distributed on [0,1)-interval, for each
   call. sgenrand(seed) sets initial values to the working area
   of 624 words. Before genrand(), sgenrand(seed) must be
   called once. (seed is any 32-bit integer.)
   Integer generator is obtained by modifying two lines.
     Coded by Takuji Nishimura, considering the suggestions by
   Topher Cooper and Marc Rieffel in July-Aug. 1997.

   Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.
   When you use this, send an email to: matumoto@math.keio.ac.jp
   with an appropriate reference to your work.

   REFERENCE
   M. Matsumoto and T. Nishimura,
   "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
   Pseudo-Random Number Generator",
   ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1, January 1998, pp 3--30.
*/

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

MersenneTwisterRNG::MersenneTwisterRNG(unsigned int seed, NormKind norm_kind)
  : RmathRNG("base::Mersenne-Twister", norm_kind), mt(dummy+1), mti(N+1)
{
  /* mti==N+1 means mt[N] is not initialized */
  init(seed);
}

void MersenneTwisterRNG::MT_sgenrand(unsigned int seed)
{
  int i;
  
  for (i = 0; i < N; i++) {
    mt[i] = seed & 0xffff0000;
    seed = 69069 * seed + 1;
    mt[i] |= (seed & 0xffff0000) >> 16;
    seed = 69069 * seed + 1;
  }
  mti = N;
}

/* Theoretically, there are 2^19937-1 possible states as an intial state.
   Essential bits in "seed_array[]" is following 19937 bits:
    (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1].
   (seed_array[0]&LOWER_MASK) is discarded.
   Theoretically,
    (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1]
   can take any values except all zeros.                             */

double MersenneTwisterRNG::uniform()
{
  unsigned int y;
  static const unsigned int mag01[2]={0x0, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
  mti = dummy[0];

  if (mti >= N) { /* generate N words at one time */
    int kk;

    if (mti == N+1)   /* if init() has not been called, */
      MT_sgenrand(4357); /* a default initial seed is used   */

    for (kk = 0; kk < N - M; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
      mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (; kk < N - 1; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
      mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

    mti = 0;
  }

  y = mt[mti++];
  y ^= TEMPERING_SHIFT_U(y);
  y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
  y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
  y ^= TEMPERING_SHIFT_L(y);
  dummy[0] = mti;

  /* reals: [0,1)-interval */
  return fixup( (double)y * 2.3283064365386963e-10 ); 
}

void MersenneTwisterRNG::init(unsigned int seed)
{
    /* Initial scrambling */
    for(unsigned int j = 0; j < 50; j++)
        seed = (69069 * seed + 1);

    for(unsigned int j = 0; j < 625; j++) {
       seed = (69069 * seed + 1);
       dummy[j] = seed;
    }
    fixupSeeds(true);
}

void MersenneTwisterRNG::fixupSeeds(bool initial)
{	
  bool notallzero = false;
  
  if(initial) dummy[0] = 624;
  /* No action unless user has corrupted .Random.seed */
  if(dummy[0] <= 0) dummy[0] = 624;
}

static bool notAllZero(unsigned int const *I)
{
  /* Check that not all seeds are zero. This is an invalid seed.
     Whereas R randomizes the seed based on the current time stamp,
     we want to return an error when this happens.
  */

  for (unsigned int j = 1; j <= 624; j++) {
    if(I[j] != 0) {
    }
  }
  return false;
}

bool MersenneTwisterRNG::setState(std::vector<int> const &state)
{
  if (state.size() != 625)
    return false;
  
  for (unsigned int j = 0; j < 625; ++j) {
    dummy[j] = static_cast<unsigned int>(state[j]);
  }
  fixupSeeds(false);

  return notAllZero(dummy);
}

void MersenneTwisterRNG::getState(std::vector<int> &state) const
{
  state.clear();
  state.reserve(625);
  for (unsigned int j = 0; j < 625; ++j) {
    state[j] = static_cast<int>(dummy[j]);
  }
}
