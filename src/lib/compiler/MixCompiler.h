#ifndef MIX_COMPILER_H_
#define MIX_COMPILER_H_

class Node;
class ParseTree;
class Compiler;

Node * getMixtureNode(ParseTree const * var, Compiler *compiler);

#endif
