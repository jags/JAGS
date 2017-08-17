#include <graph/MixTab.h>
#include <graph/Node.h>

#include <stdexcept>

using std::vector;
using std::map;
using std::logic_error;

namespace jags {

    static SimpleRange mkRange(map<vector<int>, Node const *> const &mixmap)
    {
	/* 
	   Calculates the smallest SimpleRange enclosing the index
	   values in the MixMap.  Also checks for consistency in the
	   lengths of the index values.
	*/
	map<vector<int>, Node const *>::const_iterator p = mixmap.begin();

	unsigned int N = p->first.size();
	vector <int> lower(p->first), upper(p->first);
	
	for (++p ; p != mixmap.end(); ++p) {
	    if (p->first.size() != N) {
		throw logic_error("index size mismatch in MixTab");
	    }
	    for (unsigned int j = 0; j < N; ++j) {
		int i = p->first[j];
		if (i < lower[j]) lower[j] = i;
		if (i > upper[j]) upper[j] = i;
	    }
	}
	
	return SimpleRange(lower, upper);
    }

    MixTab::MixTab(map<vector<int>, Node const *> const &mixmap)
	: _range(mkRange(mixmap)), _nodes(_range.length(), 0)
    {
	for (map<vector<int>, Node const *>::const_iterator p = mixmap.begin();
	     p != mixmap.end(); ++p)
	{
	    _nodes[_range.leftOffset(p->first)] = p->second;
	}
    }
    
    Node const * MixTab::getNode(vector<int> const &index) const
    {
	//FIXME: This needs an API change. See src/lib/SimpleRange.cc
	//if (!_range.contains(index)) return 0;
	unsigned int offset = _range.leftOffset(index);
	return _nodes[offset];
    }

    
    Range const &MixTab::range() const
    {
	return _range;
    }


/*
    double MixTab::density() 
    {
	MixTab::const_iterator p = _map.begin(); 
	vector <int> lower(p->first), upper(p->first);
	unsigned int N = p->first.size();
	
	// Find a range containing the indices
	for (++p ; p != _map.end(); ++p) {
	    for (unsigned int j = 0; j < N; ++j) {
		int i = p->first[j];
		if (i < lower[j]) lower[j] = i;
		if (i > upper[j]) upper[j] = i;
	    }
	}

	// Density is the number of indices actually used divided by
	// the number of possible indices in the enclosing range.
	double numerator = _map.size();
	double denominator = Range(lower, upper).length(); 
	return numerator / denominator;
    }
*/

}

