#ifndef KL_BIN_H_
#define KL_BIN_H_

#include "KL.h"

namespace dic {

    class KLBin : public KL
    {
	double divergence(std::vector<double const*> const &par0,
			  std::vector<double const*> const &par1) const;
	
    };

}

#endif /* KL_BIN_H_ */
