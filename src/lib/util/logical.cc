#include <util/logical.h>

#include <set>

using std::vector;
using std::set;

namespace jags {
    
    vector<bool> const * getUnique(vector<bool> const &mask)
    {
	static set<vector<bool> > _maskset;
	return &*_maskset.insert(mask).first;
    }

}
