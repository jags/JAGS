#include <config.h>
#include <graph/DeterministicNode.h>
//#include <graph/NodeError.h>
#include <stdexcept>

using std::set;
using std::logic_error;

DeterministicNode::DeterministicNode(Index const &dim, std::vector<Node*> const &parents)
  : Node(dim, parents)
{}
