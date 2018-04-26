#include <config.h>
#include <sarray/Range.h>
#include <util/dim.h>

#include <stdexcept>
#include <sstream>

using std::vector;
using std::out_of_range;
using std::logic_error;
using std::string;
using std::ostringstream;

static vector<unsigned long> 
makeDim(vector<vector<unsigned long> > const &scope)
{
    vector<unsigned long> dims;
    for (unsigned long i = 0; i < scope.size(); i++) {
	dims.push_back(scope[i].size());
    }
    return dims;
}

namespace jags {

    Range::Range()
	: _scope(0), _dim(0), _dim_dropped(0), _first(0), _last(0), _length()
    {
    }

    Range::Range(vector<vector<unsigned long> > const &scope)
	: _scope(scope),
	  _dim(makeDim(_scope)), 
	  _dim_dropped(drop(_dim)),
	  _first(0),
	  _last(0),
	  _length(product(_dim))
    {
	//Check that we have no zero-length index vectors.  The
	//leftIndex and rightIndex member functions rely on this
	//assumption, as well as the RangeIterator class.

	for (unsigned long i = 0; i < scope.size(); ++i) {
	    if (scope[i].empty()) {
		throw logic_error("Zero dimension in Range constructor");
	    }
	    _first.push_back(scope[i].front());
	    _last.push_back(scope[i].back());
	}
    }

    Range::~Range()
    {
	//Virtual destructor
    }

    unsigned long Range::length() const
    {
	return _length;
    }

    unsigned long Range::ndim(bool drop) const
    {
	return drop ? _dim_dropped.size() : _dim.size();
    }

    vector<unsigned long> const &Range::dim(bool drop) const
    {
	return drop ? _dim_dropped : _dim;
    }
    
    vector<unsigned long> Range::leftIndex(unsigned long offset) const
    {
	if (offset >= _length) {
	    throw out_of_range("Range::leftIndex. Offset out of bounds");
	}
	
	unsigned long ndim = _scope.size();
	vector<unsigned long> index(ndim);
	for (unsigned long i = 0; i < ndim; ++i) {
	    index[i] = _scope[i][offset % _dim[i]];
	    offset /= _dim[i];
	}
	return index;
    }

    vector<unsigned long> Range::rightIndex(unsigned long offset) const
    {
	if (offset >= _length) {
	    throw out_of_range("Range::rightIndex. Offset out of bounds");
	}

	unsigned long ndim = _scope.size();
	vector<unsigned long> index(ndim);
	for (unsigned long j = ndim; j > 0; --j) {
	    unsigned long i = j - 1;
	    index[i] = _scope[i][offset % _dim[i]];
	    offset /= _dim[i];
	}
	return index;
    }

    vector<unsigned long> const &Range::first() const
    {
	return _first;
    }

    vector<unsigned long> const &Range::last() const
    {
	return _last;
    }

    bool Range::operator==(Range const &range) const
    {
	return _scope == range._scope;
    }

    bool Range::operator!=(Range const &range) const
    {
	return _scope != range._scope;
    }

    bool Range::operator<(Range const &rhs) const
    {
	//Sort by the first element
	if (_first < rhs._first) {
	    return true;
	}
	else if (rhs._first < _first) {
	    return false;
	}
	//Then last element
	else if (_last < rhs._last) {
	    return true;
	}
	else if (rhs._last < _last) {
	    return false;
	}
	//Then lexigraphic ordering of the scope
	return _scope < rhs._scope;
    }
	
    string printRange(Range const &range)
    {
	if (isNULL(range))
	    return "";
	
	ostringstream ostr;
	ostr << "[";
	for (unsigned long i = 0; i < range.ndim(false); ++i) {
	    if (i > 0) {
		ostr << ",";
	    }
	    vector<unsigned long> const &indices = range.scope()[i];
	    ostr << indices[0];
	    if (indices.size() > 1) {
		bool simple = true;
		unsigned long val = indices[0] + 1;
		for (unsigned long j = 1; j < indices.size(); ++j) {
		    if (indices[j] != val) {
			simple = false;
			break;
		    }
		    ++val;
		}
		if (simple) {
		    ostr << ":";
		}
		else {
		    ostr << "...";
		}
		ostr << indices.back();
	    }
	}
	ostr << "]";
	return ostr.str();
    }

    string printRange(vector<unsigned long> const &dim)
    {
	if (dim.empty()) return "";
	
	ostringstream ostr;
	ostr << "[";
	for (unsigned long i = 0; i < dim.size(); ++i) {
	    if (i > 0) {
		ostr << ",";
	    }
	    if (dim[i] == 1) {
		ostr << "1";
	    }
	    else {
		ostr << "1:" << dim[i];
	    }
	}
	ostr << "]";
	return ostr.str();
    }
    
    string printIndex(vector<unsigned long> const &index)
    {
	if (index.empty()) return "";
	
	ostringstream ostr;
	ostr << "[";
	for (unsigned long i = 0; i < index.size(); ++i) {
	    if (i > 0) {
		ostr << ",";
	    }
	    ostr << index[i];
	}
	ostr << "]";
	return ostr.str();
    }

    
    vector<vector<unsigned long> > const &Range::scope() const
    {
	return _scope;
    }

} //namespace jags
