#include <graph/DevianceNode.h>
#include <distribution/Distribution.h>
#include <sarray/Index.h>
#include <stdexcept>

using std::vector;
using std::set;
using std::logic_error;

static vector<Node*> mkParents(set<StochasticNode*> const &param)
{
   vector <Node*> par;
   unsigned long N = param.size();
   par.reserve(N);
   set<StochasticNode*>::const_iterator p;
   for (p = param.begin(); p != param.end(); ++p) {
      par.push_back(*p);
   }

   if (par.empty()) {
      throw logic_error("Cannot construct a deviance node with no parameters");
   }
   return par;
}
         
DevianceNode::DevianceNode(set<StochasticNode*> const &parameters)
  : DeterministicNode(Index(1), mkParents(parameters))
{
  _parameters.reserve(parameters.size());
  for (set<StochasticNode*>::const_iterator p = parameters.begin(); 
       p != parameters.end(); ++p) 
    {
      _parameters.push_back(*p);
    }
}

void DevianceNode::deterministicSample(unsigned int chain)
{
  double deviance = 0;
  for (vector<StochasticNode const*>::const_iterator p = _parameters.begin();
       p != _parameters.end(); ++p)
    {
       deviance -= (*p)->logDensity(chain);
    }
  deviance *= 2;
  setValue(&deviance, 1, chain);
}

bool DevianceNode::isLinear(std::set<Node*> const &parameters, bool fixed) const 
{
  return false;
}
bool DevianceNode::isScale(std::set<Node*> const &parameters, bool fixed) const 
{
  return false;
}
