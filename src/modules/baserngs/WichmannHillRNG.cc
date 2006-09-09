#include "WichmannHillRNG.h"

void WichmannHillRNG::fixupSeeds()
{
  I[0] = I[0] % 30269; I[1] = I[1] % 30307; I[2] = I[2] % 30323;

  /* map values equal to 0 mod modulus to 1. */
  if(I[0] == 0) I[0] = 1;
  if(I[1] == 0) I[1] = 1;
  if(I[2] == 0) I[2] = 1;
}

WichmannHillRNG::WichmannHillRNG(unsigned int seed, NormKind normkind)
    : RmathRNG("base::Wichmann-Hill", normkind)
{
  init(seed);
}

double WichmannHillRNG::uniform()
{
  I[0] = I[0] * 171 % 30269;
  I[1] = I[1] * 172 % 30307;
  I[2] = I[2] * 170 % 30323;
  double value = I[0] / 30269.0 + I[1] / 30307.0 + I[2] / 30323.0;
  return fixup(value - (int) value); /* in [0,1) */
}

void WichmannHillRNG::init(unsigned int seed)
{
  /* Initial scrambling */
  for(unsigned int j = 0; j < 50; j++)
    seed = (69069 * seed + 1);
  
  
  for(unsigned int j = 0; j < 3; ++j) {
    seed = (69069 * seed + 1);
    I[j] = seed;
  }
  fixupSeeds();
}

bool WichmannHillRNG::setState(std::vector<int> const &state)
{
  if (state.size() != 3)
    return false;

  for (unsigned int i = 0; i < 3; ++i) {
    I[i] = static_cast<unsigned int>(state[i]);
  }
  fixupSeeds();
  return true;
}

void WichmannHillRNG::getState(std::vector<int> &state) const
{
  state.clear();
  for (unsigned int i = 0; i < 3; ++i) {
    state.push_back(static_cast<int>(I[i]));
  }
}
