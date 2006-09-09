#include <config.h>
#include <function/FuncTab.h>
#include <function/Function.h>
#include <function/InverseLinkFunc.h>

#include <iostream>

using std::map;
using std::string;

FuncTab::~FuncTab ()
{
    if (!_table.empty()) {
	std::cerr << "Warning: FuncTab contains " << _table.size()
		  << " functions when destructor called" << std::endl;
    }
    if (!_link_table.empty()) {
	std::cerr << "Warning: FuncTab contains " << _link_table.size()
		  << " link functions when destructor called" << std::endl;
    }
}

bool
FuncTab::insert (Function const *func)
{
  string const &name = func->name ();
  if (_table.find (name) == _table.end ())
    {
      _table[name] = func;
      return true;
    }
  else
    {
      return false;
    }
}

bool
FuncTab::insert (InverseLinkFunc const *lfunc)
{
  string const &name = lfunc->name ();
  string const &linkname = lfunc->linkName ();
  if (_link_table.find (linkname) == _link_table.end () &&
      _table.find (name) == _table.end ())
    {
      _table[name] = lfunc;
      _link_table[linkname] = lfunc;

      return true;
    }
  else
    {
      return false;
    }
}

bool
FuncTab::erase (string const &name)
{
    /* Is it a link function?
       Manual search through the map of link functions required, since they
       are indexed by link name */
    map < const string, const InverseLinkFunc *>::iterator q;
    for (q = _link_table.begin(); q != _link_table.end(); ++q) {
	if (q->second->name() == name) {
	    _link_table.erase(q);
	    _table.erase(name);
	    return true;
	}
    }

    return( _table.erase(name) > 0);
}


Function const *
FuncTab::find (string const &name) const
{
  map < const string, const Function *>::const_iterator 
    p = _table.find (name);

  if (p == _table.end ())
    {
      return 0;
    }
  else
    {
      return p->second;
    }
}

Function const *
FuncTab::findInverse (string const &name) const
{
    map < const string, const InverseLinkFunc *>::const_iterator 
	p =_link_table.find (name);
  
    if (p == _link_table.end ())
	{
	    return 0;
	}
    else
	{
	    return p->second;
	}
}


void FuncTab::erase(Function *func)
{
    string const &name = func->name();
    map<const string, const Function*>::iterator p = _table.find(name);
    if (p != _table.end()) {
	if (p->second == func) {
	    erase(name);
	}
    }
}
