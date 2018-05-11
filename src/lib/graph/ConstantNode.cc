#include <config.h>
#include <graph/ConstantNode.h>

#include <sstream>
#include <cmath>
#include <stdexcept>

using std::vector;
using std::string;
using std::ostringstream;
using std::floor;
using std::logic_error;
using std::array;

namespace jags {

    ConstantNode::ConstantNode(double value, unsigned int nchain,
			       bool observed)
	: Node(vector<unsigned long>(1,1), nchain), _observed(observed)
    {
	for (unsigned int n = 0; n < nchain; ++n) {
	    setValue(&value, 1, n);
	}
    }

    ConstantNode::ConstantNode(vector<unsigned long> const &dim, 
			       vector<double> const &value,
			       unsigned int nchain, bool observed)
	: Node(dim, nchain), _observed(observed)
    {
	if (value.size() != _length) {
	    throw logic_error("Invalid value in ConstantNode");
	}

	for (unsigned int n = 0; n < nchain; ++n) {
	    setValue(&value[0], _length, n);
	}
    }

    array<int, 2> const &ConstantNode::depth() const
    {
	static const std::array<int, 2> _depth = {0, 0};
	return _depth;
    }
    
    bool ConstantNode::isDiscreteValued() const
    {
	double const *val = value(0);
	for (unsigned long i = 0; i < _length; ++i) {
	    if (val[i] != floor(val[i]))
		return false;
	}
	return true;
    }

    void ConstantNode::randomSample(RNG*, unsigned int) {}

    string ConstantNode::deparse(vector<string> const &parents) const
    {
	ostringstream os;
	if (length() == 1) {
	    os << *value(0);
	}
	else {
	    os << "c(" << value(0)[0] << "..."  << value(0)[length() - 1] << ")";
	}
	return os.str();
    }

    bool ConstantNode::checkParentValues(unsigned int) const
    {
	return true;
    }

    bool ConstantNode::isFixed() const
    {
	return true;
    }

    RVStatus ConstantNode::randomVariableStatus() const
    {
	return _observed ? RV_TRUE_OBSERVED :  RV_FALSE;
    }

    void ConstantNode::unlinkParents()
    {

    }
    
} //namespace jags
