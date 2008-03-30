#ifndef KL_NORM_H_
#define KL_NORM_H_

#include "KL.h"

namespace dic {

    class KLNorm : public KL
    {
	double divergence(std::vector<double const*> const &par0,
			  std::vector<double const*> const &par1) const;
	
    };

}

#endif /* KL_NORM_H_ */
