#include <config.h>
#include <graph/GraphMarks.h>
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <distribution/DistError.h>
#include <distribution/Distribution.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <vector>
#include <stdexcept>
#include <string>
#include <stdexcept>
#include <cmath>
#include <set>

using std::vector;
using std::string;
using std::runtime_error;
using std::logic_error;
using std::set;

namespace jags {

static vector<Node const *> mkParents(vector<Node const *> const &parameters,
                                      Node const *lower, Node const *upper)
{
    //Add bounds to vector of parents, if they are non-zero
    vector<Node const *> parents = parameters;
    if (lower) {
        parents.push_back(lower);
    }
    if (upper) {
        parents.push_back(upper);
    }
    return parents;
}

    static bool mkDiscrete(Distribution const *dist,
			   vector<Node const *> const &parameters)
    {
	//Determine whether node is discrete-valued
	vector<bool> mask(parameters.size());
	for (unsigned long i = 0; i < parameters.size(); ++i) {
	    mask[i] = parameters[i]->isDiscreteValued();
	}
	if (!dist->checkParameterDiscrete(mask)) {
	    throw DistError(dist, 
			    "Failed check for discrete-valued parameters");
	}
	return(dist->isDiscreteValued(mask));
    }

StochasticNode::StochasticNode(vector<unsigned long> const &dim,
			       unsigned int nchain,
			       Distribution const *dist,
			       vector<Node const *> const &parameters,
			       Node const *lower, Node const *upper)
    : Node(dim, nchain, mkParents(parameters, lower, upper)), 
      _dist(dist), _lower(lower), _upper(upper), 
      _observed(false), 
      _discrete(mkDiscrete(dist, parameters)), _parameters(nchain)
{
    if (!checkNPar(dist, parameters.size())) {
	throw DistError(_dist, "Incorrect number of parameters");
    }

    //check boundaries
    if ((lower && lower->dim() != this->dim()) || 
	(upper && upper->dim() != this->dim()))
    {
	throw DistError(_dist, "Dimension mismatch when setting bounds");
    }
    if (!_dist->canBound() && (lower || upper)) {
	throw DistError(_dist, "Distribution cannot be bounded");
    }

    //Set up parameter vectors 
    for (unsigned int n = 0; n < nchain; ++n) {
	_parameters[n].reserve(parameters.size());
	for (unsigned long i = 0; i < parameters.size(); ++i) {
	    _parameters[n].push_back(parameters[i]->value(n));
	}
    }

    //Insert the current node as a child in all its parents
    //taking care to avoid repeats
    set<Node const*> pset;
    for (unsigned long i = 0; i < parents().size(); ++i) {
	Node const *p = parents()[i];
	if (pset.insert(p).second) {
	    p->addChild(this);
	}
    }

}

StochasticNode::~StochasticNode()
{
}

Distribution const *StochasticNode::distribution() const
{
    return _dist;
}

bool StochasticNode::isDiscreteValued() const
{
    return _discrete;
}

bool StochasticNode::isFixed() const
{
    return _observed;
}

RVStatus StochasticNode::randomVariableStatus() const
{
    return _observed ? RV_TRUE_OBSERVED : RV_TRUE_UNOBSERVED;    
}	

string StochasticNode::deparse(vector<string> const &parnames) const
{
    unsigned long npar = parnames.size();
    if (_upper) --npar;
    if (_lower) --npar;
    if (!checkNPar(_dist, npar)) {
	//We might be deparsing after throwing a NodeError, so we 
	//don't want to throw another exception.
	return _dist->name() + "(deparse error)";
    }

    string name = _dist->name();
    name.append("(");
    unsigned long i = 0; 
    for ( ; i < npar; ++i) {
	if (i != 0) {
	    name.append(",");
	}
	name.append(parnames[i]);
    }
    name.append(")");
    if (_lower || _upper) {
        name.append(" T(");
        if (_lower) {
            name.append(parnames[i++]);
        }
        name.append(",");
        if (_upper) {
            name.append(parnames[i++]);
        }
        name.append(")");
    }

    return name;
}

Node const *StochasticNode::lowerBound() const
{
    return _lower;
}

Node const *StochasticNode::upperBound() const
{
    return _upper;
}

/*
StochasticNode * 
StochasticNode::clone(vector<Node const *> const &parents) const
{
    vector<Node const *> param(parents);

    //Peel off upper and lower bounds
    Node const *lower = lowerBound();
    Node const *upper = upperBound();
    if (upper) {
	param.pop_back();
    }
    if (lower) {
	param.pop_back();
    }
    
    return clone(param, lower, upper);
}
*/
    
double const *StochasticNode::lowerLimit(unsigned int chain) const
{
    return _lower ? _lower->value(chain) : 0;
}

double const *StochasticNode::upperLimit(unsigned int chain) const
{
    return _upper ? _upper->value(chain) : 0;
}


bool isSupportFixed(StochasticNode const *node)
{
    if (node->lowerBound() && !node->lowerBound()->isFixed())
	return false;
    if (node->upperBound() && !node->upperBound()->isFixed())
	return false;
    
    vector<Node const *> parents = node->parents();
    if (node->upperBound())
	parents.pop_back();
    if (node->lowerBound())
	parents.pop_back();
    vector<bool> fixmask(parents.size());
    for (unsigned long i = 0; i < parents.size(); ++i) {
	fixmask[i] = parents[i]->isFixed();
    }
    return node->distribution()->isSupportFixed(fixmask);
}

void StochasticNode::support(double *lower, double *upper, unsigned long length,
			     unsigned int chain) const
{
    if (length != _length) {
        throw logic_error("Length mismatch in StochasticNode support");
    }

    // Get support without bounds
    sp(lower, upper, length, chain);

    if (_lower || _upper) {
        if (!distribution()->canBound()) {
            throw logic_error("Bounded node has non-boundable distribution");
        }
        if (_lower) {
            double const *lb = _lower->value(chain);
            for (unsigned long i = 0; i < length; ++i) {
                if (lower[i] < lb[i])
                    lower[i] = lb[i];
            }
        }
        if (_upper) {
            double const *ub = _upper->value(chain);
            for (unsigned long i = 0; i < length; ++i) {
                if (upper[i] > ub[i])
                    upper[i] = ub[i];
            }
        }
    }
}

bool isBounded(StochasticNode const *node)
{
    return node->lowerBound() || node->upperBound();
}

    void StochasticNode::setData(double const *value, unsigned long length)
    {
	for (unsigned int n = 0; n < _nchain; ++n) {
	    setValue(value, length, n);
	}
	_observed = true;
    }

    void StochasticNode::unlinkParents()
    {
	for (unsigned long i = 0; i < parents().size(); ++i) {
	    parents()[i]->removeChild(this);
	}
    }
    
} //namespace jags
