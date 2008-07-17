#ifndef FUNC_ICLOGLOG_H_
#define FUNC_ICLOGLOG_H_

#include <function/InverseLinkFunc.h>

namespace bugs {

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
	double evaluateScalar(std::vector<double const *> const &args) const;
	double link(double mu) const;
	double grad(double eta) const;
    };

}

#endif /* FUNC_ICLOGLOG_H_ */
