#ifndef _WICHMANN_HILL_RNG_H_
#define _WICHMANN_HILL_RNG_H_

#include <rng/RmathRNG.h>

class WichmannHillRNG : public RmathRNG
{
    unsigned int I[3];
    void fixupSeeds();
public:
    WichmannHillRNG(unsigned int seed, NormKind norm_kind);
    double uniform();
    void init(unsigned int seed);
    bool setState(std::vector<int> const &state);
    void getState(std::vector<int> &state) const;
};

#endif
