#include <config.h>
#include <distribution/Distribution.h>
#include <graph/StochasticNode.h>
#include <graph/Graph.h>
#include <sampler/GraphView.h>
#include <rng/RNG.h>

#include "RealDSum.h"

//Initial step size
#define STEP 0.1

using std::vector;
using std::string;

namespace jags {
namespace bugs {

//Pick a random integer between 0 and n - 1
static unsigned int pick(unsigned long n, RNG *rng)
{
    double u = rng->uniform() * n;
    unsigned int i = 1;
    while (i < u) ++i;
    return i - 1;
}

RealDSum::RealDSum(GraphView const *gv, unsigned int chain)
    : RWDSum(gv, chain, STEP)
{
}

void RealDSum::step(vector<double> &value, unsigned long nrow,
		    unsigned long ncol, double s, RNG *rng) const
{
    //Randomly pick a row
    unsigned int r = pick(nrow, rng);

    //Randomly draw two columns 
    unsigned int c1 = pick(ncol, rng);
    unsigned int c2 = pick(ncol - 1, rng);
    if (c2 >= c1) ++c2;
    
    //Modify the chosen components while keeping the sum constant
    double eps = rng->normal() * s;
    value[c1 * nrow + r] += eps;
    value[c2 * nrow + r] -= eps;
}

}}
