#ifndef KL_BERN_H_
#define KL_BERN_H_

#include "KL.h"

namespace dic {

    class KLBern : public KL
    {
	double divergence(std::vector<double const*> const &par0,
			  std::vector<double const*> const &par1) const;
	
    };

}

#endif /* KL_BERN_H_ */
