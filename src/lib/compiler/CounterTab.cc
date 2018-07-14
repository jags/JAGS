#include <config.h>
#include <compiler/CounterTab.h>

#include <string>

using std::vector;
using std::pair;
using std::string;

namespace jags {

CounterTab::CounterTab()
{
}

CounterTab::~CounterTab()
{
    unsigned long n = _table.size();
    for (unsigned int i = 0; i < n; i++) {
	popCounter();
    }
}

Counter * CounterTab::pushCounter(string const &name, Range const &range)
{
  Counter *counter = new Counter(range);
  pair<string, Counter*> cpair(name, counter);
  _table.push_back(cpair);
  return counter;
}

void CounterTab::popCounter()
{
  pair<string, Counter *> topcpair = _table.back();
  _table.pop_back();
  delete topcpair.second;
}

Counter *CounterTab::getCounter(string const &name) const
{
  vector<pair<string, Counter *> >::const_iterator p(_table.begin());

  for (; p != _table.end(); ++p) {
    if (name == p->first)
      return p->second;
  }
  return nullptr;
}

} //namespace jags
