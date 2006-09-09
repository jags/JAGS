#include <config.h>
#include <compiler/ConstantFactory.h>
#include <graph/ConstantNode.h>

#include <sstream>
#include <cmath>

using std::map;
using std::ostringstream;
using std::fabs;

bool lt(double arg1, double arg2)
{
  if (equal(arg1,arg2))
    return false;
  else
    return arg1 < arg2;
}

ConstantFactory::ConstantFactory(unsigned int nchain)
  : _nchain(nchain)
{
}

ConstantNode*
ConstantFactory::getConstantNode(double value)
{
  
  ConstantNode *cnode = 0;

  map<const double,ConstantNode*, ltdouble>::const_iterator i 
    = _constmap.find(value);
  if (i == _constmap.end()) {
    // Create a new constant node
    //ostringstream os;
    long ivalue;
    /* FIXME: This is, unfortunately, a big problem casting doubles to
       long that is used everywhere, and is currently broken for
       negative numbers
    */
    if (value >= 0) {
       ivalue = static_cast<long>(value + DBL_EPSILON);
    }
    else {
      ivalue = static_cast<long>(value - DBL_EPSILON);
    }
    if (fabs(value - ivalue) < DBL_EPSILON) {
      // Integer value 
      cnode = new ConstantNode(ivalue, _nchain);
      //os << ivalue;
    }
    else {
      // Floating point value
      cnode = new ConstantNode(value, _nchain);
      //os << value;
    }

    //cnode->setName(os.str());
    _graph.add(cnode);
    _constmap[value] = cnode;
    return cnode;
  }
  else {
    return i->second;
  }
}
