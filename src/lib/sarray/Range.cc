#include <config.h>
#include <sarray/Range.h>
#include "RangeImp.h"

#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <map>

using std::map;
using std::pair;
using std::string;
using std::ostringstream;

struct ltrangeimp
{
  bool operator()(RangeImp const *arg1, RangeImp const *arg2) const
  {
    return *arg1 < *arg2;
  }
};

map<RangeImp const*, long, ltrangeimp> &impMap()
{
  static map<RangeImp const*, long, ltrangeimp> _impmap;
  return _impmap;
}

inline long &ref(RangeImp const *p) { return impMap().find(p)->second; }

static RangeImp const *getPtr(Index const &lower, Index const &upper)
{
  RangeImp * newimp = new RangeImp(lower, upper);  
  map<RangeImp const*, long, ltrangeimp> &impmap = impMap();
  map<RangeImp const*, long, ltrangeimp>::iterator p = impmap.find(newimp);
  if (p == impmap.end()) {
    impmap.insert(pair<RangeImp*,long>(newimp,0));
    return newimp;
  }
  else {
    delete newimp;
    return p->first;
  }
}

Range::Range()
  : _p(getPtr(Index(0),Index(0)))
{
    ref(_p)++;
}

Range::Range(Index const &lower, Index const &upper)
  : _p(getPtr(lower,upper))
{
  ref(_p)++;
}

Range::Range(Range const &rhs)
  : _p(rhs._p)
{
  ref(_p)++;
}

Range::Range(Index const &index)
  : _p(getPtr(index,index))
{
  ref(_p)++;
}

Range &Range::operator=(Range const &rhs)
{
  if (_p == rhs._p)
    return *this; //Nothing to do

  if (--ref(_p) == 0) {
    impMap().erase(_p);
    delete _p;
  }
  _p = rhs._p;
  ref(_p)++;
  return *this;
}

Range::~Range()
{
  if (--ref(_p) == 0) {
    impMap().erase(_p);
    delete _p;
  }
}

bool Range::operator==(Range const &range) const
{ 
  return _p == range._p; 
}

bool Range::operator!=(Range const &range) const
{ 
  return _p != range._p; 
}

unsigned long Range::length() const 
{ 
  return _p->length(); 
}

bool Range::contains(Range const &range) const 
{ 
  return _p->contains(*range._p); 
}

long Range::leftOffset(Index const &index) const 
{ 
  return _p->leftOffset(index); 
}

Index Range::leftIndex(long n) const
{
  return _p->leftIndex(n);
}

long Range::rightOffset(Index const &index) const 
{ 
  return _p->rightOffset(index); 
}

Index Range::rightIndex(long n) const
{
  return _p->rightIndex(n);
}

Index const &Range::dim(bool drop) const 
{ 
  return _p->dim(drop); 
}

unsigned int Range::ndim(bool drop) const 
{ 
  return _p->ndim(drop); 
}

Index const & Range::lower() const 
{ 
  return _p->lower(); 
}

Index const & Range::upper() const 
{ 
  return _p->upper(); 
}

bool Range::operator<(Range const &rhs) const { 
  return *_p < *rhs._p; 
}

string print(Range const &range)
{
  ostringstream ostr;
  ostr << "[";
  for (unsigned int i = 0; i < range.ndim(false); ++i) {
    if (i > 0)
      ostr << ",";
    if (range.lower()[i] == range.upper()[i]) {
      ostr << range.lower()[i];
    }
    else {
      ostr << range.lower()[i] << ":" << range.upper()[i];
    }
  }
  ostr << "]";
  return ostr.str();
}
