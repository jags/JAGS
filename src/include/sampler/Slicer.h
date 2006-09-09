#ifndef SLICER_H_
#define SLICER_H_

#include <sampler/Sampler.h>

class StochasticNode;

/**
 * Base class for univariate slice sampler (Neal R, Slice Sampling,
 * Ann. Statist.  31 (2003), no. 3, 705-767)
 *
 * The Slicer class is written with a high level of abstraction to
 * allow it to be easily adapted: for example, a slice sampler can be
 * used to sample from a discrete distribution or from a univariate
 * "slice" of a multivariate distribution.
 */
class Slicer : public Sampler 
{
  double _width;
  bool _burnin;
  long _max;
  double _sumdiff;
  long _iter;
  bool accept(double xold, double xnew, double z, double L, double R);
 public:
  /**
   * Constructor for Slice Sampler
   * @param node Node to sample
   * @param Graph containing node within which sampling takes place.
   * @param width Initial width of slice
   * @param max Maximum number of times initial width of slice will increase.
   */
  Slicer(std::vector<StochasticNode*> const &nodes, Graph const &graph,
	 unsigned int chain, double width,  long max);
  /**
   * Update the current value using the "stepping" method
   */
  void updateStep(RNG *rng);
  /**
   * Update the current value using the "doubling" method
   */
  void updateDouble(RNG *rng);
  /**
   * Returns the current value of the sampler
   */
  virtual double value() = 0;
  /**
   * Sets the value of the sampler, and additionally updates all the
   * deterministic descendants
   */
  virtual void setValue(double x) = 0;
  /**
   * Returns the lowest possible value of the sampler
   */
  virtual double lowerLimit() = 0;
  /**
   * Returns the highest possible value of the sampler
   */
  virtual double upperLimit() = 0;
  void burninOff();

};

#endif /* SLICER_H_ */
