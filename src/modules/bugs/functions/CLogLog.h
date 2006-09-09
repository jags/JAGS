#ifndef FUNC_CLOGLOG_H_
#define FUNC_CLOGLOG_H_

#include <function/ScalarFunc.h>

/**
 * @short the complementary log log link function
 * @see ICLogLog
 * <pre>
 * y <- cloglog(p)
 * y = log(-log(1-p)) for 0 < p < 1
 * </pre>
 */
class CLogLog : public ScalarFunc
{
public:
  CLogLog ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterValue (std::vector <SArray const *> const &args) const;
};

#endif /* FUNC_CLOGLOG_H_ */
