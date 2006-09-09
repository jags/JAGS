#ifndef FUNC_ICLOGLOG_H_
#define FUNC_ICLOGLOG_H_

#include <function/InverseLinkFunc.h>

/**
 * @short inverse complementary log log link
 * @see CLogLog
 * <pre>
 * cloglog(y) <- a + b*x
 * y <- icloglog(a + b*x)
 * </pre>
 */
class ICLogLog : public InverseLinkFunc
{
public:
  ICLogLog ();
  double eval (std::vector <SArray const *> const &args) const;
  double link (double y) const;
  double gradLink (double y) const;
};

#endif /* FUNC_ICLOGLOG_H_ */
