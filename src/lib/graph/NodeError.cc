#include <config.h>
#include <graph/NodeError.h>
#include <model/SymTab.h>

using std::string;
using std::ostream;
using std::endl;

namespace jags {

    NodeError::NodeError(Node const *node, string const &msg)
	: runtime_error(msg), _node(node)
    {
    }

    void NodeError::printMessage(ostream &out, SymTab const &symtab) const
    {
	out << "Error in node " <<  symtab.getName(_node) << "\n" 
	    << what() << endl;
    }

} //namespace jags
