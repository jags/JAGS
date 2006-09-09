#ifndef REAL_SLICE_SAMPLER_H_
#define REAL_SLICE_SAMPLER_H_

#include <sampler/Slicer.h>

class StochasticNode;
struct RNG;

/**
 * Slice sampler for real-valued distributions
 */
class RealSliceSampler : public Slicer 
{
public:
  /**
   * Constructor for Slice Sampler
   * @param node Node to sample
   * @param width Initial width of slice
   * @param maxwidth Maximal width of slice as a multiple of the width
   * parameter
   * @param nburn Length of burnin
   */
  RealSliceSampler(StochasticNode *node, Graph const &graph,
		   unsigned int chain, double width = 1, long maxwidth = 10);
  static bool canSample(StochasticNode const *node, Graph const &graph);
  double value();
  void setValue(double value);
  double lowerLimit();
  double upperLimit();
  void update(RNG *rng);
};

#endif /* REAL_SLICE_SAMPLER_H_ */
