#include <config.h>
#include <distribution/DistTab.h>
#include <distribution/Distribution.h>
#include <sarray/SArray.h>

#include <iostream>

using std::map;
using std::string;

DistTab::~DistTab()
{
    if (!_table.empty()) {
	std::cerr << "Warning: DistTab contained " << _table.size()
		  << "Distributions when destructor called" << std::endl;
    }
}

bool DistTab::insert(Distribution const *dist)
{
  string const &name = dist->name();
  if (_table.find(name) == _table.end()) {
    _table[name] = dist;
    return true;
  }
  else {
    return false;
  }
}

bool DistTab::erase(string const &name)
{
  map<const string, const Distribution*>::iterator p = _table.find(name);
  if (p == _table.end()) {
    return false;
  }
  else {
    _table.erase(p);
    return true;
  }
}

Distribution const *DistTab::find(string const &name) const
{
  map<const string, const Distribution*>::const_iterator p = _table.find(name);

  if (p == _table.end())
    return 0;
  else
    return p->second;
}

void DistTab::erase(Distribution *dist)
{
    string const &name = dist->name();
    map<const string, const Distribution*>::iterator p = _table.find(name);
    if (p != _table.end()) {
	if (p->second == dist) {
	    _table.erase(p);
	}
    }
}
