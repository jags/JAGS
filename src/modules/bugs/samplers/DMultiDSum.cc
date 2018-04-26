#include <config.h>
#include <distribution/Distribution.h>
#include <graph/StochasticNode.h>
#include <graph/Graph.h>
#include <sampler/GraphView.h>
#include <rng/RNG.h>

#include "DMultiDSum.h"

#include <cmath>

using std::vector;
using std::fabs;
using std::string;

//Initial step size
#define STEP 1

namespace jags {
namespace bugs {

DMultiDSum::DMultiDSum(GraphView const *gv, unsigned int chain)
    : RWDSum(gv, chain, STEP)
{
}

//Pick a random integer between 0 and n - 1
static unsigned int pick(unsigned long n, RNG *rng)
{
    double u = rng->uniform() * n;
    unsigned int i = 1;
    while (i < u) ++i;
    return i - 1;
}

void DMultiDSum::step(vector<double> &value, 
		      unsigned long nrow, unsigned long ncol,
		      double s, RNG *rng) const
{
    //Randomly pick two rows
    unsigned int r1 = pick(nrow, rng);
    unsigned int r2 = pick(nrow - 1, rng);
    if (r2 >= r1) ++r2;
    
    //Randomly pick two columns
    unsigned int c1 = pick(ncol, rng);
    unsigned int c2 = pick(ncol - 1, rng);
    if (c2 >= c1) ++c2;
    
    double eps = rng->normal() * s;
    int inteps = static_cast<int>(fabs(eps)) + 1;
    
    //Modify the chosen elements while keeping row and column sums
    //constant
    value[nrow * c1 + r1] += inteps;
    value[nrow * c1 + r2] -= inteps;
    value[nrow * c2 + r1] -= inteps;
    value[nrow * c2 + r2] += inteps;
}

}}
