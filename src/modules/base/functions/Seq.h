#ifndef SEQ_H_
#define SEQ_H_

#include <function/VectorFunction.h>

namespace jags {
    namespace base {

    /**
     * @short Integer sequence
     *
     * The Seq function implements the infix operator ":" in the BUGS language.
     * The arguments must be integer.
     * <pre>
     * y <- a:b
     * </pre>
     *
     */
    class Seq : public VectorFunction
    {
    public:
	Seq();
	void evaluate(double *value, 
		      std::vector <double const *> const &values,
		      std::vector <unsigned long> const &lengths) const;
	unsigned long length(std::vector <unsigned long> const &lengths,
			     std::vector <double const *> const &values) const;
	bool checkParameterLength(std::vector<unsigned long> const &args) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
	bool checkParameterDiscrete(std::vector<bool> const &mask) const;
	bool checkParameterFixed(std::vector<bool> const &mask) const;
	std::string deparse(std::vector<std::string> const &par) const;
    };

    }
} /* namespace jags */

#endif /* SEQ_H_ */
