#ifndef FUNC_MULTIPLY_H_
#define FUNC_MULTIPLY_H_

#include "Infix.h"

namespace jags {
namespace base {

    /**
     * @short Product of two scalars 
     * 
     * The Multiply function implements the infix operator "*" in the
     * BUGS language.  The expression a * b * c is evaluated in a single 
     * function call.
     */
    class Multiply: public Infix
    {
    public:
	Multiply ();
	double evaluate(std::vector<double const *> const &args) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
	bool isScale(std::vector<bool> const &mask,
		     std::vector<bool> const &fixmask) const;
	bool isPower(std::vector<bool> const &mask, 
		     std::vector<bool> const &fix) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;

    };
    
}}

#endif /* FUNC_MULTIPLY_H_ */
