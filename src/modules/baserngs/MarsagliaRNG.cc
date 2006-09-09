#include <config.h>
#include "MarsagliaRNG.h"

using std::vector;

#define i2_32m1 2.328306437080797e-10/* = 1/(2^32 - 1) */

MarsagliaRNG::MarsagliaRNG(unsigned int seed, NormKind normkind)
    : RmathRNG("base::Marsaglia-Multicarry", normkind)
{
   init(seed);
}

void MarsagliaRNG::fixupSeeds()
{
  if (I[0] == 0) I[0] = 1;
  if (I[1] == 0) I[1] = 1;
}

void MarsagliaRNG::init(unsigned int seed)
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

bool MarsagliaRNG::setState(vector<int> const &state)
{
  if (state.size() != 2)
    return false;

  I[0] = static_cast<unsigned int>(state[0]);
  I[1] = static_cast<unsigned int>(state[1]);
  fixupSeeds();
  return true;
}

void MarsagliaRNG::getState(vector<int> &state) const
{
  state.clear();
  state.push_back(static_cast<int>(I[0]));
  state.push_back(static_cast<int>(I[1]));
}

double MarsagliaRNG::uniform()
{
  I[0]= 36969*(I[0] & 0177777) + (I[0]>>16);
  I[1]= 18000*(I[1] & 0177777) + (I[1]>>16);
  /* result in in [0,1) */
  return fixup(((I[0] << 16)^(I[1] & 0177777)) * i2_32m1); 
}
