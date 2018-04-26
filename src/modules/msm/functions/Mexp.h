#ifndef MEXP_H_
#define MEXP_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace msm {

class Mexp : public ArrayFunction
{
public:
    Mexp();
    void evaluate (double *value, std::vector<double const *> const &args,
		   std::vector<std::vector<unsigned long> > const &dims) const;
    std::vector<unsigned long> 
	dim(std::vector<std::vector<unsigned long> > const &dims,
	    std::vector<double const *> const &values) const;
    bool checkParameterDim(std::vector <std::vector<unsigned long> > const &dims) const;
};

}}

#endif /* MEXP_H_ */
