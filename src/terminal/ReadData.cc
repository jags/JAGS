#include <model/NodeArray.h>
#include "ReadData.h"
#include <sarray/SArray.h>

#include <iostream>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

using std::map;
using std::string;
using std::vector;

bool readRData(std::vector<ParseTree*> const *array_list, 
	       std::map<std::string, SArray> &table,
	       std::string &rngname)
{
    /* Check validity of expressions */
    for (vector<ParseTree*>::const_iterator p = array_list->begin(); 
	 p != array_list->end(); ++p) 
      {
	if ((*p)->treeClass() == P_VAR) {
	  /*
	    Assignments of the form "foo" <- "bar" The only type
	    currently allowed is ".RNG.name" <- "bar"
	  */
	  if ((*p)->name() != ".RNG.name") {
	    cout << "Unrecognized string assignment. "
		 << "Expecting \".RNG.name\"";
	    return false;
	  }
		
	  if (((*p)->parameters().size() != 1) ||
	      ((*p)->parameters()[0]->treeClass() != P_VAR)) 
	    {
	      cout << "Invalid .RNG.name";
	      return false;
	    }
	      
	  rngname = (*p)->parameters()[0]->name();
	}
	else if ((*p)->treeClass() != P_ARRAY) {
	  cout << "Error reading S data.";
	  return false;
	}
      }  

    for (vector<ParseTree*>::const_iterator p = array_list->begin(); 
	 p != array_list->end(); ++p) 
	{
	  if ((*p)->treeClass() == P_VAR) {
	    /* Skip any string assignments and break if one comes
	       at the end */
	    ++p;
	    if (p == array_list->end())
	      break;
	  }

	  string const &name = (*p)->name();

	  /* Check to see if name is already in table */
	  if (table.find(name) != table.end()) {
	    cerr << "WARNING: Replacing " << name << endl;
	    table.erase(table.find(name));
	  }

	  /* Get the length of the data */
	  ParseTree const *vec = (*p)->parameters()[0]; 
	  unsigned long length = vec->parameters().size();
    
	  /* Get the number of dimensions of the array */
	  ParseTree const *pdim;
	  long ndim;
	  if ((*p)->parameters().size() == 2) {
	    // Array has dimension attribute
	    pdim = (*p)->parameters()[1];
	    ndim = pdim->parameters().size();
	  }
	  else {
	    // No dimension attribute
	    pdim = 0;
	    ndim = 1;
	  }

	  /* Get the dimensions of the array */
	  Index dim(ndim);
	  if (pdim) {
	    for (int i = 0; i < ndim; ++i) {
	      dim[i] = (long) (pdim->parameters()[i]->value() + 1.0E-6);
	    }
	    /* Check that dimension is consistent with length */
	    unsigned long dimprod = 1;
	    for (int i = 0; i < ndim; i++) {
	      dimprod *= dim[i];
	    }
	    if (dimprod != length) {
	      cerr << "Bad dimension for variable " << name << endl;
	      return false;
	    }
	  }
	  else {
	    dim[0] = length;
	  }

	  /* Get the data */
	  double *values = new double[length];
	  for (unsigned long i = 0; i < length; ++i) {
	    values[i] = vec->parameters()[i]->value();
	  }

	  /* Now assign it to an SArray */
	  SArray sarray(dim);
	  sarray.setValue(values, length);
	  delete [] values;
    
	  /* Since there is no default constructor for SArray, we can't
	     use the shorthand table[names[i]] = par;
	  */
	  table.insert(map<const string, SArray>::value_type(name, sarray));

	  /*
	    std::cout << "Reading " << name << "[";
	    for(unsigned int j = 0; j < dim.length(); j++) {
	    if (j > 0) {
	    std::cout << ",";
	    }
	    std::cout << dim[j];
	    }
	    std::cout << "]" << std::endl;
	  */
	}
    return true;
}










