#include "SuperDuperRNG.h"

#define I1 I[0]
#define I2 I[1]

#define i2_32m1 2.328306437080797e-10/* = 1/(2^32 - 1) */

SuperDuperRNG::SuperDuperRNG(unsigned int seed, NormKind norm_kind)
  : RmathRNG("base::Super-Duper", norm_kind)
{
   init(seed);    
}

void SuperDuperRNG::fixupSeeds()
{
  if(I1 == 0) I1 = 1;
  /* I2 = Congruential: must be ODD */
  I2 |= 1;
}

double SuperDuperRNG::uniform()
{
  /* This is Reeds et al (1984) implementation;
   * modified using __unsigned__	seeds instead of signed ones
   */
  I1 ^= ((I1 >> 15) & 0377777); /* Tausworthe */
  I1 ^= I1 << 17;
  I2 *= 69069;		/* Congruential */
  return fixup((I1^I2) * i2_32m1); /* in [0,1) */
}


void SuperDuperRNG::init(unsigned int seed)
{
  /* Initial scrambling */
  for(unsigned int j = 0; j < 50; j++)
    seed = (69069 * seed + 1);
  
  for (unsigned int j = 0; j < 2; j++) {
    seed = (69069 * seed + 1);
    I[j] = seed;
  }
  fixupSeeds();
}

bool SuperDuperRNG::setState(std::vector<int> const &state)
{
  if (state.size() != 2) 
    return false;
  
  for (unsigned int j = 0; j < 2; j++) {
    I[j] = static_cast<unsigned int>(state[j]);
  }
  fixupSeeds();
  return true;
}


void SuperDuperRNG::getState(std::vector<int> &state) const
{
  state.clear();
  for (unsigned int j = 0; j < 2; j++) {
    state.push_back(static_cast<int>(I[j]));
  }
}
