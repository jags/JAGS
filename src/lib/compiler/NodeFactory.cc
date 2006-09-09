#include <compiler/NodeFactory.h>
#include <graph/Node.h>
//FIXME: For JAGS_NA only
#include <sarray/SArray.h>

#include <cfloat>

static const double tol = 16 * DBL_EPSILON;

NodeFactory::~NodeFactory()
{}

Graph const &NodeFactory::graph() const
{
  return _graph;
}

bool equal(double arg1, double arg2)
{
  if (arg1 == JAGS_NA || arg2 == JAGS_NA)
    return false;
  else 
    return !(arg1 + tol < arg2 || arg2 + tol < arg1);
}

bool equal(Node const *node1, Node const *node2)
{
  if (node1 == node2) {
    return true;
  }
  else if (node1->isObserved() && node2->isObserved() &&
	   node1->length() == node2->length()) {
    unsigned long n = node1->length();
    double const *value1 = node1->value(0);
    double const *value2 = node2->value(0);
    for (unsigned long i = 0; i < n; ++i) {
      if (!equal(value1[i],value2[i])) {
	return false;
      }
    }
    return true;
  }
  else {
    return false;
  }
}
