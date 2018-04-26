#include <config.h>
#include <util/dim.h>

#include <set>

using std::vector;
using std::set;

namespace jags {

    vector<unsigned long> const &getUnique(vector<unsigned long> const &dim)
    {
	static set<vector<unsigned long> > _dimset;
	return *(_dimset.insert(dim).first);
    }

    vector<vector<unsigned long> > const &
    getUnique(vector<vector<unsigned long> > const &dimvec)
    {
	static set<vector<vector<unsigned long> > > _dimvecset;
	return *(_dimvecset.insert(dimvec).first);
    }

    unsigned long product(vector<unsigned long> const &x)
    {
	if (x.empty())
	    return 0;

	unsigned long y = x[0];
	for (unsigned long i = 1; i < x.size(); ++i) {
	    y *= x[i];
	}
	return y;
    }

    vector<unsigned long> drop(vector<unsigned long> const &dims)
    {
	vector<unsigned long> ans;
	bool empty = true;
	for (unsigned long i = 0; i < dims.size(); ++i) {
	    if (dims[i] != 0) {
		empty = false;
	    }
	    if (dims[i] != 1) {
		ans.push_back(dims[i]);
	    }
	}
	if (ans.empty() && !empty) {
	    ans.push_back(1U);
	}

	return ans;
    }

}
