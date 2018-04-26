#include <config.h>
#include <model/Monitor.h>
#include <graph/StochasticNode.h>
#include <graph/Node.h>
#include <util/dim.h>

#include <stdexcept>

using std::string;
using std::vector;
using std::logic_error;
using std::copy;

namespace jags {

Monitor::Monitor(string const &type, vector<Node const *> const &nodes)
    : _type(type), _nodes(nodes)
{
}

Monitor::Monitor(string const &type, Node const *node)
    : _type(type), _nodes(vector<Node const*>(1,node))
{
}

Monitor::~Monitor()
{}

string const &Monitor::type() const
{
    return _type;
}

vector<Node const*> const &Monitor::nodes() const
{
    return _nodes;
}

string const &Monitor::name() const
{
    return _name;
}

void Monitor::setName(string const &name)
{
    _name = name;
}

vector<string> const &Monitor::elementNames() const
{
    return _elt_names;
}

void Monitor::setElementNames(vector<string> const &names)
{
    _elt_names = names;
}

SArray Monitor::dump(bool flat) const
{
    unsigned int nchain = poolChains() ? 1 : nodes()[0]->nchain();
    unsigned long nvalue = value(0).size();

    vector<double> v(nvalue * nchain);
    vector<double>::iterator p = v.begin();
    for (unsigned int ch = 0; ch < nchain; ++ch) {
	p = copy(value(ch).begin(), value(ch).end(), p);
    }

    vector<unsigned long> vdim = dim();
    unsigned long vlen = product(vdim);
    if (nvalue % vlen != 0) {
	throw logic_error("Inconsistent dimensions in Monitor");
    }
    unsigned long niter = nvalue / vlen;
    if (poolIterations() && niter != 1) {
	throw logic_error("Invalid number of iterations in Monitor");
    }

    if (flat) {
	vdim = vector<unsigned long>(1, vlen);
    }
	
    vector<string> names(vdim.size(), "");

    if (!poolIterations()) {
	vdim.push_back(niter);
	names.push_back("iteration");
    }
    if (!poolChains()) {
	vdim.push_back(nchain);
	names.push_back("chain");
    }
	
    SArray ans(vdim);
    ans.setValue(v);    
    ans.setDimNames(names);
    if (flat) {
	ans.setSDimNames(_elt_names, 0);
    }
    return(ans);
}

} //namespace jags
