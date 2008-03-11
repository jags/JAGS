#include <graph/DevianceNode.h>
#include <distribution/Distribution.h>
#include <stdexcept>

using std::vector;
using std::set;
using std::logic_error;
using std::string;

static vector<Node const *> mkParents(set<StochasticNode const *> const &param)
{
   vector <Node const *> par;
   unsigned long N = param.size();
   par.reserve(N);
   set<StochasticNode const *>::const_iterator p;
   for (p = param.begin(); p != param.end(); ++p) {
      par.push_back(*p);
      /* The parents of a deviance node include not only the given
	 parameters, but also their parents */
      vector <Node const *> const &parents = (*p)->parents();
      for (unsigned int i = 0; i < parents.size(); ++i) {
	  par.push_back(parents[i]);
      }
   }

   if (par.empty()) {
      throw logic_error("Cannot construct a deviance node with no parameters");
   }
   return par;
}
         
DevianceNode::DevianceNode(set<StochasticNode const *> const &parameters)
  : DeterministicNode(vector<unsigned int>(1,1), mkParents(parameters))
{
  _parameters.reserve(parameters.size());
  for (set<StochasticNode const *>::const_iterator p = parameters.begin(); 
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

bool DevianceNode::isLinear(GraphMarks const &linear_marks, bool fixed) const
{
    return false;
}

bool DevianceNode::isScale(GraphMarks const &scale_marks, bool fixed) const
{
    return false;
}

bool DevianceNode::checkParentValues(unsigned int nchain) const
{
    return true;
}

string DevianceNode::deparse(vector<string> const &parents) const
{
    return "deviance";
}
