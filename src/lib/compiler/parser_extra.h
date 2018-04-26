#ifndef PARSER_EXTRA_H_
#define PARSER_EXTRA_H_

#include <compiler/ParseTree.h>
#include <cstdio>
#include <string>

/**
 * Parse a model file that describes the graphical model using the
 * BUGS language.  
 * 
 * @param file File containing the model description
 *
 * @param pvariables Pointer to a vector of ParseTree pointers. After
 * calling parse_bugs pvariables will point to a newly allocated
 * ParseTree representing the list of declared variables.  
 *
 * @param pdata Pointer to a ParseTree. After calling parse_bugs pdata
 * will point to a newly allocated ParseTree representing the
 * stochastic and logical relationships in the data block, if one exists.
 *
 * @param pmodel Pointer to a ParseTree. After calling parse_bugs pdata
 * will point to a newly allocated ParseTree representing the
 * stochastic and logical relationships in the data block, if one exists.
 *
 * @param message String that contains any error messages from the 
 * parser on exit.
 *
 * @return 0 on success and 1 on error.
 *
 */
int parse_bugs(std::FILE *file, std::vector<jags::ParseTree*> * &pvariables, 
	       jags::ParseTree * &pdata, jags::ParseTree * &pmodel,
	       std::string &message);

#endif /* PARSER_EXTRA_H_ */
