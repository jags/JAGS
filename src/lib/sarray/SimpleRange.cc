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

static vector<vector<unsigned long> > 
makeScope(vector<unsigned long> const &lower,
	  vector<unsigned long> const &upper)
{
    unsigned long ndim = lower.size();
    if (upper.size() != ndim) {
	throw logic_error("Dimension mismatch in Range constructor");
    }
    
    vector<vector<unsigned long> > scope(ndim);
    for (unsigned long i = 0; i < ndim; ++i) {
	if (lower[i] <= upper[i]) {
	    for (unsigned long j = lower[i]; j <= upper[i]; ++j) {
		scope[i].push_back(j);
	    }
	}
    }
    return scope;
}

namespace jags {

    SimpleRange::SimpleRange(vector<unsigned long> const &lower,
			     vector<unsigned long> const &upper)
	: Range(makeScope(lower, upper))
    {
    }

    /*
    SimpleRange::SimpleRange(vector<long> const &index)
	: Range(makeScope(index, index))
    {
    }
    */
    
    SimpleRange::SimpleRange(vector<unsigned long> const &dim)
	: Range(makeScope(vector<unsigned long>(dim.size(), 1), dim))
    {
    }

    SimpleRange::SimpleRange()
	: Range()
    {
    }

    SimpleRange::~SimpleRange()
    {
    }

    bool SimpleRange::contains(vector<unsigned long> const &index) const
    {
	unsigned long N = ndim(false);
	if (N != index.size()) {
	    return false;
	}
	
	for (unsigned long i = 0; i < N; ++i) {
	    if (index[i] < _first[i] || index[i] > _last[i]) {
		return false;
	    }
	}
	return true;
    }
    
    bool SimpleRange::contains(SimpleRange const &other) const
    {
	unsigned long N = ndim(false);
	if (N != other.ndim(false)) {
	    return false;
	}
	
	for (unsigned long i = 0; i < N; ++i) {
	    if (other._first[i] < _first[i] || other._last[i] > _last[i]) {
		return false;
	    }
	}
	return true;
    }


    bool SimpleRange::contains(Range const &other) const
    {
	unsigned long ndim = scope().size();
	if (other.scope().size() != ndim) {
	    throw invalid_argument("SimpleRange::contains. Dimension mismatch");
	}
	
	for (unsigned long i = 0; i < ndim; ++i) {
	    vector<unsigned long> const &indices = other.scope()[i];
	    for (unsigned long j = 0; j < indices.size(); ++j) {
		if (indices[j] < _first[i] || indices[j] > _last[i]) {
		    return false;
		}
	    }
	}
	return true;
    }
    
    unsigned long
    SimpleRange::leftOffset(vector<unsigned long> const &index) const
    {
	unsigned long offset = 0;
	unsigned long step = 1;
	unsigned long ndim = _last.size();
	for (unsigned long i = 0; i < ndim; i++) {
	    if (index[i] < _first[i] || index[i] > _last[i]) {
		throw out_of_range("SimpleRange::leftOffset. Index outside of allowed range");
	    }
	    offset += step * (index[i] - _first[i]);
	    step *= _dim[i];
	}
	return offset;
    }

    unsigned long
    SimpleRange::rightOffset(vector<unsigned long> const &index) const
    {
	unsigned long offset = 0;
	unsigned long step = 1;
	for (unsigned long j = _last.size(); j > 0; --j) {
	    unsigned long i = j - 1;
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

	vector<unsigned long> const & lower = range.lower();
	vector<unsigned long> const & upper = range.upper();
	ostringstream ostr;
	ostr << "[";
	for (unsigned long i = 0; i < range.ndim(false); ++i) {
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
