#include <config.h>
#include <graph/AggNode.h>

#include <vector>
#include <stdexcept>

using std::vector;
using std::set;

AggNode::AggNode(Index const &dim, 
		 vector<Node *> const &nodes,
                 vector<unsigned long> const &offsets)
  : DeterministicNode(dim, nodes), _offsets(offsets)
{
  /* Check argument lengths */
  if (length() != nodes.size() || length() != offsets.size()) {
    throw std::length_error ("Length mismatch in Aggregate Node constructor");
  }

  /* Check that offsets are valid */
  for (unsigned long i = 0; i < length(); i++) {
    if (offsets[i] >= parents()[i]->length())
      throw std::out_of_range("Invalid offset in Aggregate Node constructor");
  }
  
  /* See if node is discrete-valued */
  bool isdiscrete = true;
  for (unsigned long i = 0; i < length(); i++) {
    if (!parents()[i]->isDiscreteValued()) {
      isdiscrete = false;
      break;
    }
  }
  if (isdiscrete) {
    setDiscreteValued();
  }
}

AggNode::~AggNode()
{
}

void AggNode::deterministicSample(unsigned int chain)
{
  unsigned long l = this->length();
  double *value = new double[l];
  vector<Node*> const &par = parents();
  for (unsigned long i = 0; i < l; ++i) {
    value[i] = par[i]->value(chain)[_offsets[i]];
  }
  setValue(value, l, chain);
  delete [] value;
}

AggNode const *asAggregate(Node *node)
{
  return dynamic_cast<AggNode const*>(node);
}

bool AggNode::isLinear(set<Node*> const &parameters, bool fixed) const
{
   return true;
}

bool AggNode::isScale(set<Node*> const &parameters, bool fixed) const
{
    return (parents().size() == 1);
}
