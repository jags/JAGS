#ifndef DSUM_H_
#define DSUM_H_

#include <distribution/Distribution.h>

/**
 * @short Sum of two discrete random variables
 */
class DSum : public Distribution {
 public:
  DSum();
  ~DSum();
  
  double logLikelihood (SArray const & x, 
		        std::vector<SArray const *> const &parameters) const;
  void randomSample (SArray &x, std::vector<SArray const *> const &parameters,
		     RNG *rng)
    const;
  bool checkParameterDim (std::vector<Index> const &dims) const;
  bool checkParameterDiscrete (std::vector<bool> const &mask) const;
  bool checkParameterValue(std::vector<SArray const *> const &parameters) 
    const;
  Index dim(std::vector<Index> const &dims) const;
  unsigned long df(std::vector<SArray const *> const &parameters) const;
  double lowerSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  double upperSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  void typicalValue(SArray &x, std::vector<SArray const *> const &parameters)
    const;
};

#endif /* DSUM_H_ */
