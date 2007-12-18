#include <config.h>
#include <distribution/DistTab.h>
#include <distribution/Distribution.h>
#include <functional>
#include <algorithm>

using std::string;
using std::binary_function;
using std::find_if;

typedef std::list<Distribution const*> DistList;

/* Adaptable binary predicate for find_if algorithm */
struct isName: 
    public binary_function<Distribution const *, string const *, bool> 
{
    bool operator()(Distribution const *dist, string const *name) const
    {
	return dist->name() == *name;
    }
};

void DistTab::insert(Distribution const *dist)
{
    DistList::iterator p = find_if(_dist_list.begin(), _dist_list.end(),
                                   bind2nd(isName(), &dist->name()));
    //Transfer distribution with same name (if any) to the masked list
    if (p != _dist_list.end()) {
	_masked_list.push_front(*p);
	_dist_list.erase(p);
    }
    _dist_list.push_front(dist);
}

Distribution const *DistTab::find(string const &distname) const
{
    DistList::const_iterator p = find_if(_dist_list.begin(), 
					      _dist_list.end(),
					      bind2nd(isName(), &distname));

    return (p == _dist_list.end()) ? 0 : *p;
}

void DistTab::erase(Distribution *dist)
{
    //Erase from the main list.
    //Note scoping operator to distinguish std::find from DistTab::find
    DistList::iterator p = std::find(_dist_list.begin(), _dist_list.end(),
				     dist);
    bool move_masked = false;
    if (p != _dist_list.end()) {
	_dist_list.erase(p);
	move_masked = true; //Move any masked distributions into the main list
    }

    //Erase from the masked list
    p = std::find(_masked_list.begin(), _masked_list.end(), dist);
    if (p != _masked_list.end()) {
        if (move_masked) {
	    _dist_list.push_front(*p);
	}
	_masked_list.erase(p);
    }
}

DistList const &DistTab::distributions() const
{
    return _dist_list;
}

DistList const &DistTab::masked() const
{
    return _masked_list;
}
