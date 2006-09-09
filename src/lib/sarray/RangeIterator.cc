#include <config.h>
#include <sarray/RangeIterator.h>
#include <stdexcept>

using std::out_of_range;

RangeIterator::RangeIterator(Range const &range)
  : Index(range.lower()), _range(range), _atend(false)
{
}

RangeIterator::~RangeIterator()
{}

RangeIterator &RangeIterator::operator=(Index const &index)
{
    if (!_range.contains(index))
	throw out_of_range("Cannot set RangeIterator: index out of range");

    _atend = false;
    Index::operator=(index);
    return *this;
}

bool RangeIterator::atEnd() const
{
    return _atend;
}

RangeIterator &RangeIterator::nextLeft()
{
    int n = _range.ndim(false);
    int i = 0;
    for (; i < n; ++i) {
	long &ind = operator[](i);
	if (ind < _range.upper()[i]) {
	    ++ind;
	    break;
	}
	else {
	    ind = _range.lower()[i];
	}
    }
    if (i == n) {
	_atend = true;
    }
    return *this;
}

RangeIterator &RangeIterator::nextRight()
{
    int i = _range.ndim(false) - 1;
    for ( ; i >= 0; --i) {
	long &ind = operator[](i);
	if (ind < _range.upper()[i]) {
	    ++ind;
	    break;
	}
	else {
	    ind = _range.lower()[i];
	}
    }
    if (i < 0) {
	_atend = true;
    }
    return *this;
}
