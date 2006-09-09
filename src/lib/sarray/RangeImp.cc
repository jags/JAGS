#include <config.h>
#include "RangeImp.h"

#include <algorithm>
#include <stdexcept>
#include <sstream>

using std::min;
using std::max;
using std::out_of_range;
using std::invalid_argument;
using std::logic_error;
//using std::ostringstream;
using std::string;

static Index makeDim(Index const &lower, Index const &upper, bool drop)
{
  /* 
     Utility function for constructor: calculates the dimension of the
     range with, or without, dropping of scalar dimensions.
  */
  unsigned int ndim = upper.size();

  if (!drop) {
    //We only need to check this once
    if (lower.size() != ndim) {
      throw logic_error("RangeImp: Length mismatch in constructor");
    }
    else {
      for (unsigned int i = 0; i < ndim; i++) {
	if (upper[i] < lower[i]) {
	  throw out_of_range("RangeImp: upper < lower bound in constructor");
      }
      }
    }
  }

  long *dim = new long[ndim];

  int size = 0;
  for (unsigned int i = 0; i < ndim; i++) {
    long d = (upper[i] - lower[i] + 1);
    if (!drop || d > 1) {
      dim[size++] = d;
    }
  }

  if (drop && size == 0) {
    delete [] dim;
    return Index(1); //FIXME: Why?
  }
  else {
    Index index(size);
    for (int j = 0; j < size; j++) {
      index[j] = dim[j];
    }
    delete [] dim;
    return index;
  }
}

static unsigned long makeLength(Index const &lower, Index const &upper)
{
  /* 
     Utility function for constructor: calculates the length of the
     range.
  */
  unsigned int ndim = upper.size();
  if (ndim == 0) {
    return 0;
  }
  else {
    unsigned long length = 1;
    for (unsigned int i = 0; i < ndim; i++) {
      length *= (upper[i] - lower[i] + 1);
    }
    return length;
  }
}

RangeImp::RangeImp(Index const &lower, Index const &upper)
  : _lower(lower), _upper(upper), 
	   _dim(makeDim(lower, upper, false)),
	   _dim_dropped(makeDim(lower, upper, true)), 
	   _length(makeLength(lower, upper))
    
{
}

unsigned long RangeImp::length() const
{
  return _length;
}

unsigned int RangeImp::ndim(bool drop) const
{
  if (drop)
    return _dim_dropped.size();
  else
    return _dim.size();
}

Index const &RangeImp::dim(bool drop) const
{
  if (drop)
    return _dim_dropped;
  else
    return _dim;
}

bool RangeImp::contains(RangeImp const &other) const
{
  unsigned int ndim = _upper.size();
  if (other._lower.size() != ndim) {
    throw invalid_argument("Range::contains. Dimension mismatch");
  }
 
  for (unsigned int i = 0; i < ndim; ++i) {
    if (other._lower[i] < _lower[i] || other._upper[i] > _upper[i])
      return false;
  }
  return true;
}

bool RangeImp::overlaps(RangeImp const &other) const
{
  unsigned int dim = _upper.size();
  if (other.ndim(false) != dim) {
    throw invalid_argument("RangeImp::overlaps. Dimension mismatch");
  }

  for (unsigned int i = 0; i < dim; i++) {
    if (max(_lower[i], other._lower[i]) > min(_upper[i], other._upper[i])) {
      return false;
    }
  }
  return true;
}

unsigned long RangeImp::leftOffset(Index const &index) const
{
  if (!contains(RangeImp(index,index))) {
    throw out_of_range("RangeImp::leftOffset. Index outside of allowed range");
  }
  
  unsigned long offset = 0;
  long step = 1;
  unsigned int ndim = _upper.size();
  for (unsigned int i = 0; i < ndim; i++) {
      offset += step * (index[i] - _lower[i]);
      step *= (_upper[i] - _lower[i] + 1);
  }
  return offset;
}

Index RangeImp::leftIndex(unsigned long offset) const
{
  if (offset >= _length) {
    throw out_of_range("RangeImp::leftIndex. Offset exceeds length of range");
  }

  unsigned int ndim = _lower.size();
  Index index(_lower);
  for (unsigned int i = 0; i < ndim; ++i) {
    index[i] += offset % _dim[i];
    offset = offset / _dim[i];
  }
  return index;
}

unsigned long RangeImp::rightOffset(Index const &index) const
{
  if (!contains(RangeImp(index,index))) {
    throw out_of_range("RangeImp::rightOffset.Index outside of allowed range");
  }
 
  unsigned long offset = 0;
  long step = 1;
  for (int i = _upper.size() - 1; i >= 0; --i) {
    offset += step * (index[i] - _lower[i]);
    step *= (_upper[i] - _lower[i] + 1);
  }
  return offset;
}

Index RangeImp::rightIndex(unsigned long offset) const
{
  if (offset >= _length)
    throw out_of_range("RangeImp::rightIndex. Offset exceeds length of range");

  int ndim = _lower.size();
  Index index(_lower);
  for (int i = ndim - 1; i >= 0; --i) {
    index[i] += offset % _dim[i];
    offset = offset / _dim[i];
  }
  return index;
}

Index const & RangeImp::lower() const
{
  return _lower;
}

Index const & RangeImp::upper() const
{
  return _upper;
}

bool RangeImp::operator<(RangeImp const &rhs) const
{
  if (_lower < rhs._lower)
    return true;
  else if (rhs._lower < _lower)
    return false;
  else 
    return _upper < rhs._upper;
}
