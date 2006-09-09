#ifndef _MARSAGLIA_RNG_H_
#define _MARSAGLIA_RNG_H_

#include <rng/RmathRNG.h>

class MarsagliaRNG : public RmathRNG
{
    unsigned int I[2];
    void fixupSeeds();
public:
    MarsagliaRNG(unsigned int seed, NormKind norm_kind);
    void init(unsigned int seed);
    bool setState(std::vector<int> const &state);
    void getState(std::vector<int> &state) const;
    double uniform();
};

#endif
