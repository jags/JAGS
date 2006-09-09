#ifndef READ_DATA_H_
#define READ_DATA_H_

#include <map>
#include <string>
#include <compiler/ParseTree.h>
#include <sarray/SArray.h>

bool readRData(std::vector<ParseTree*> const *r_assignment_list, 
	       std::map<std::string, SArray> &table,
               std::string &rngname);

#endif /* READ_DATA_H_ */
