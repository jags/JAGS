#ifndef RNG_H_
#define RNG_H_

#include <string>
#include <vector>

/**
 * Abstract class for Random Number Generator
 *
 * @short Random number generator
 */
class RNG
{
  const std::string _name;
 public:
  RNG(std::string const &name);
  virtual ~RNG();
  virtual void init(unsigned int seed) = 0;
  virtual void getState(std::vector<int> &state) const = 0;
  virtual bool setState(std::vector<int> const &state) = 0;
  virtual double uniform() = 0;
  virtual double normal() = 0;
  virtual double exponential() = 0;
  static double fixup(double x);
  std::string const &name() const;
};

#endif /* RNG_H_ */
