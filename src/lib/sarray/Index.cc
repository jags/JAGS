#include <config.h>
#include <sarray/Index.h>

#include <stdexcept>
#include <algorithm>

using std::out_of_range;
using std::min;

Index::Index()
  : _size(0), _index(0)
{}

Index::Index(unsigned int size)
  : _size(size), _index(0)
{
  if (size > 0) {
    _index = new long[size];
    for (unsigned int i = 0; i < size; i++) {
      _index[i] = 1;
    }
  }
}

Index::Index(Index const &orig)
  : _size(orig._size), _index(0)
{
  if (_size > 0) {
    _index = new long[_size];
    for (unsigned int i = 0; i < _size; i++) {
      _index[i] = orig._index[i];
    }
  }
}

Index::~Index()
{
  delete [] _index;
}

Index &Index::operator=(Index const &rhs)
{
  if (rhs._size != _size) {
    _size = rhs._size;
    delete [] _index;
    _index = (_size > 0) ? new long[_size] : 0;
  }
  for (unsigned int i = 0; i < _size; i++) {
    _index[i] = rhs._index[i];
  }
  return *this;
}

unsigned int Index::size() const
{
  return _size;
}

bool Index::operator==(Index const &other) const
{
  if (_size != other._size) {
    return false;
  }
  else {
    for(unsigned int i = 0; i < _size; i++) {
      if (_index[i] != other._index[i])
	return false;
    }
  }
  return true;
}

bool Index::operator!=(Index const &other) const
{
  return !operator==(other);
}


long &Index::operator[](unsigned int i)
{
  if (i >= _size) {
    throw out_of_range("Attempt to access nonexistent element of Index");
  }

  return _index[i];
}

long Index::operator[](unsigned int i) const
{
  if (i >= _size) {
    throw out_of_range("Attempt to access nonexistent element of Index");
  }

  return _index[i];
}

bool Index::operator<(Index const &other) const
{
  unsigned int n = min(_size, other._size);
  for (unsigned int i = 0; i < n; ++i) {
    if (_index[i] < other._index[i]) {
      return true;
    }
    else if (_index[i] > other._index[i]) {
      return false;
    }
  }

  return _size < other._size;
}

