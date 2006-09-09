#include <config.h>
#include <graph/NodeError.h>

NodeError::NodeError(Node const *enode, std::string const &emsg)
  : runtime_error(emsg), node(enode)
{
}
