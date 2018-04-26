#ifndef MIX_COMPILER_H_
#define MIX_COMPILER_H_

class Node;
class ParseTree;
class Compiler;

#include <sarray/Range.h>

#include <string>
#include <utility>
#include <set>

namespace jags {

    typedef std::map<std::pair<std::string, Range>, std::set<unsigned long> > UMap;
    
    Node * getMixtureNode(ParseTree const * var, Compiler *compiler);
    void getMissingMixParams(ParseTree const * var, UMap &umap,
			     Compiler *compiler);

}

#endif /* MIX_COMPILER_H_ */
