#include <config.h>
#include <sarray/SimpleRange.h>
#include <util/dim.h>

#include <algorithm>
#include <stdexcept>
#include <sstream>

using std::vector;
using std::min;
using std::max;
using std::out_of_range;
using std::invalid_argument;
using std::logic_error;
using std::string;
using std::ostringstream;
using std::sort;

static vector<int> asSigned(vector<unsigned int> const &orig)
{
    unsigned int n = orig.size();
    vector<int> ans(n);
    for (unsigned int i = 0; i < n; ++i) {
        ans[i] = static_cast<int>(orig[i]);
	if (ans[i] < 0) {
	    throw out_of_range("Dimension too large in Range constructor");
	}
    }
    return ans;
}

static vector<vector<int> > 
makeScope(vector<int> const &lower, vector<int> const &upper)
{
    unsigned int ndim = lower.size();
    if (upper.size() != ndim) {
	throw logic_error("Dimension mismatch in Range constructor");
    }
    
    vector<vector<int> > scope(ndim);
    for (unsigned int i = 0; i < ndim; ++i) {
	if (lower[i] <= upper[i]) {
	    for (int j = lower[i]; j <= upper[i]; ++j) {
		scope[i].push_back(j);
	    }
	}
    }
    return scope;
}

namespace jags {

    SimpleRange::SimpleRange(vector<int> const &lower, vector<int> const &upper)
	: Range(makeScope(lower, upper))
    {
    }
    
    SimpleRange::SimpleRange(vector<int> const &index)
	: Range(makeScope(index, index))
    {
    }
    
    SimpleRange::SimpleRange(vector<unsigned int> const &dim)
	: Range(makeScope(vector<int>(dim.size(), 1), asSigned(dim)))
    {
    }

    SimpleRange::SimpleRange()
	: Range()
    {
    }

    SimpleRange::~SimpleRange()
    {
    }

    bool SimpleRange::contains(vector<int> const &index) const
    {
	unsigned int N = ndim(false);
	if (N != index.size()) {
	    return false;
	}
	
	for (unsigned int i = 0; i < N; ++i) {
	    if (index[i] < _first[i] || index[i] > _last[i]) {
		return false;
	    }
	}
	return true;
    }
    
    bool SimpleRange::contains(SimpleRange const &other) const
    {
	unsigned int N = ndim(false);
	if (N != other.ndim(false)) {
	    return false;
	}
	
	for (unsigned int i = 0; i < N; ++i) {
	    if (other._first[i] < _first[i] || other._last[i] > _last[i]) {
		return false;
	    }
	}
	return true;
    }


    bool SimpleRange::contains(Range const &other) const
    {
	unsigned int ndim = scope().size();
	if (other.scope().size() != ndim) {
	    throw invalid_argument("SimpleRange::contains. Dimension mismatch");
	}
	
	for (unsigned int i = 0; i < ndim; ++i) {
	    vector<int> const &indices = other.scope()[i];
	    for (unsigned int j = 0; j < indices.size(); ++j) {
		if (indices[j] < _first[i] || indices[j] > _last[i]) {
		    return false;
		}
	    }
	}
	return true;
    }
    
    unsigned int SimpleRange::leftOffset(vector<int> const &index) const
    {
	unsigned int offset = 0;
	int step = 1;
	unsigned int ndim = _last.size();
	for (unsigned int i = 0; i < ndim; i++) {
	    if (index[i] < _first[i] || index[i] > _last[i]) {
		throw out_of_range("SimpleRange::leftOffset. Index outside of allowed range");
	    }
	    offset += step * (index[i] - _first[i]);
	    step *= _dim[i];
	}
	return offset;
    }

    unsigned int SimpleRange::rightOffset(vector<int> const &index) const
    {
	unsigned int offset = 0;
	int step = 1;
	for (int i = _last.size() - 1; i >= 0; --i) {
	    if (index[i] < _first[i] || index[i] > _last[i]) {
		throw out_of_range("SimpleRange::rightOffset. Index outside of allowed range");
	    }
	    offset += step * (index[i] - _first[i]);
	    step *= _dim[i];
	}
	return offset;
    }

    bool SimpleRange::operator<(SimpleRange const &rhs) const
    {
	//Sort by the first element
	if (_first < rhs._first) {
	    return true;
	}
	else if (rhs._first < _first) {
	    return false;
	}
	//Then last element
	return _last < rhs._last;
    }

    bool SimpleRange::operator==(SimpleRange const &rhs) const
    {
	return (_first == rhs._first) && (_last == rhs._last);
    }
    
    bool SimpleRange::operator!=(SimpleRange const &rhs) const
    {
	return (_first != rhs._first) || (_last != rhs._last);
    }
    
    string print(SimpleRange const &range)
    {
	if (isNULL(range)) {
	    return "";
	}

	vector<int> const & lower = range.lower();
	vector<int> const & upper = range.upper();
	ostringstream ostr;
	ostr << "[";
	for (unsigned int i = 0; i < range.ndim(false); ++i) {
	    if (i > 0)
		ostr << ",";
	    if (lower[i] == upper[i]) {
		ostr << lower[i];
	    }
	    else {
		ostr << lower[i] << ":" << upper[i];
	    }
	}
	ostr << "]";
	return ostr.str();
    }
    
} //namespace jags
