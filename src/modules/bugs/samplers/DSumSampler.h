#ifndef DSUM_SAMPLER_H_
#define DSUM_SAMPLER_H_

#include <sampler/Slicer.h>

#include <vector>
#include <set>

class StochasticNode;
class DeterministicNode;
class Node;

/**
 * @short Sample parents of dsum nodes
 */
class DSumSampler : public Slicer
{
  double _x;
  long _sum;
public:
  DSumSampler(std::vector<StochasticNode *> const &nodes,
	      Graph const &graph, unsigned int chain);
  ~DSumSampler();
  void update(RNG *);
  bool canSample(std::vector<StochasticNode *> const &nodes, 
		 Graph const &graph);
  void setValue(double x);
  double value();
  double lowerLimit();
  double upperLimit();
};

#endif /* DSUM_SAMPLER_H_ */
