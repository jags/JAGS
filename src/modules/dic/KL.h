#ifndef KL_H_
#define KL_H_

#include <vector>

namespace dic {
    
    class KL
    {
    public:
	virtual ~KL();
	virtual double 
	    divergence(std::vector<double const*> const &param1,
		       std::vector<double const*> const &param2) const = 0;
	
    };

}

#endif /* KL_H_ */
