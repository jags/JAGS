#include <config.h>
#include <compiler/MixtureFactory.h>
#include <compiler/NodeFactory.h>
#include <model/Model.h>
#include <graph/MixtureNode.h>
#include <model/NodeArray.h>

#include <utility>
#include <stdexcept>

using std::vector;
using std::map;
using std::logic_error;

namespace jags {

MixtureNode * 
MixtureFactory::getMixtureNode(vector<Node const *> const &index_nodes,
			       MixMap const &mixmap, Model &model)
{
    /* 
       Separate the index values from the parameters. Only the latter
       are used to index the mixture nodes
    */
    vector<Node const *> parameters;
    vector<vector<unsigned long> > index_values;
    for (map<vector<unsigned long>, Node const *>::const_iterator p = mixmap.begin();
         p != mixmap.end(); ++p)
    {
        index_values.push_back(p->first);
        parameters.push_back(p->second);
    }

    MixPair mpair(index_nodes, parameters);
    map<MixPair, MixtureNode*>::const_iterator p = _mix_node_map.find(mpair);

    if (p != _mix_node_map.end()) {
	return p->second;
    }
    else {
        if (parameters.size() != index_values.size()) {
            throw logic_error("Length mismatch in MixtureFactory::getMixtureNode");
        }
	
	MixtureNode *mix = new MixtureNode(index_nodes, model.nchain(), mixmap);
	_mix_node_map[mpair] = mix;
	model.addNode(mix);
	return mix;
    }
}

} //namespace jags
