#ifndef KL_POIS_H_
#define KL_POIS_H_

#include "KL.h"

namespace dic {

    class KLPois : public KL
    {
	double divergence(std::vector<double const*> const &par0,
			  std::vector<double const*> const &par1) const;
	
    };

}

#endif /* KL_POIS_H_ */
