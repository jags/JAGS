#ifndef DMSTATE_H_
#define DMSTATE_H_

#include <distribution/DistDiscrete.h>

/**
 * Distribution of x in multi-state homogeneous Markov model at time t
 * when state at time 0 was xold and transition intensity matrix is
 * Lambda
 * <pre> 
 * x ~ dmstate(xold, dt, Lambda[,]) 
 * </pre>
 * @short Interval censored transitions in a multistate Markov model
 */
class DMState : public DistDiscrete {
public:
  DMState();
  ~DMState();
  double d(double x, std::vector<SArray const *> const &parameters, 
	   bool give_log) const;
  double p(double x, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	 bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  double l(std::vector<SArray const *> const &parameters, bool absolute) const;
  double u(std::vector<SArray const *> const &parameters, bool absolute) const;
  /**
   * Checks that x is scalar, t is scalar and Lambda is
   * a square matrix
   */
  bool checkParameterDim(std::vector<Index> const &dims) const;
  /** Checks that x is discrete */
  bool checkParameterDiscrete(std::vector<bool> const &mask) const;
  /**
   * Checks that t > 0, xold is in permitted range 1:nrow(Lambda)
   * and Lambda is a valid transition intensity matrix with
   * rows that sum to 0.
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;};

#endif /* DMSTATE_H_ */
