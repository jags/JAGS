#ifndef RMATH_RNG_H_
#define RMATH_RNG_H_

#include <rng/RNG.h>

enum NormKind {AHRENS_DIETER, BOX_MULLER, KINDERMAN_RAMAGE};

/**
 *
 * 
 * 
 */
class RmathRNG : public RNG
{
    NormKind _N01_kind;
    double _BM_norm_keep;
public:
    RmathRNG(std::string const &name, NormKind norm_kind);
    double normal();
    double exponential();
};

#endif /* RMATH_RNG_H_ */
