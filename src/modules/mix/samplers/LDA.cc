#include <config.h>

#include "LDA.h"

#include <rng/RNG.h>
#include <graph/MixtureNode.h>
#include <graph/MixTab.h>
#include <graph/StochasticNode.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>
#include <sarray/RangeIterator.h>

#include <set>
#include <vector>
#include <algorithm>
#include <numeric>

using std::vector;
using std::set;
using std::string;
using std::upper_bound;
using std::accumulate;
using std::partial_sum;

namespace jags {

     bool isCat(StochasticNode const *snode) {
	 //FIXME when dcat is split into dcat.ordinal and dcat.nominal
	 return snode->length() == 1 && snode->distribution()->name() == "dcat";
    }

     bool isDirichlet(StochasticNode const *snode) {
	return snode->distribution()->name() == "ddirch";
    }

    namespace mix {

	LDA::LDA(vector<vector<StochasticNode*> > const &topics,
		 vector<vector<StochasticNode*> > const &words,
		 vector<StochasticNode*> const &topic_priors,
		 vector<StochasticNode*> const &word_priors,
		 GraphView const *gv, unsigned int ch)
	    : _nTopic(word_priors.size()), 
	      _nWord(word_priors[0]->length()),
	      _nDoc(topics.size()), 
	      _topicHyper(topic_priors[0]->parents()[0]->value(ch)), 
	      _wordHyper(word_priors[0]->parents()[0]->value(ch)),
	      _gv(gv), 
	      _chain(ch),
	    _topicTokens(_nDoc), 
	    _wordTokens(_nDoc), 
	    _topicsByWord(_nWord, vector<int>(_nTopic, 0)),
	    _topicsByDoc(_nDoc, vector<int>(_nTopic, 0)),
	    _docSums(_nDoc), _topicSums(_nTopic),
	    _wordsObserved(true)
	{
	    //Set up tables
	    for (unsigned int d = 0; d < _nDoc; ++d) {
		_docSums[d] = topics[d].size();
		for (unsigned int i = 0; i < _docSums[d]; ++i) {
		    int topic = static_cast<int>(*topics[d][i]->value(ch)) - 1;
		    _topicTokens[d].push_back(topic);
		    _topicsByDoc[d][topic]++;
		    _topicSums[topic]++;
		    int word = static_cast<int>(*words[d][i]->value(ch)) - 1;
		    _wordTokens[d].push_back(word);
		    _topicsByWord[word][topic]++;
		    if (isParameter(words[d][i])) _wordsObserved = false;
		}
	    }

	    //Sanity check on gv: 
	    vector<StochasticNode*> const &snodes = gv->nodes();
	    unsigned int offset = 0;
	    for (unsigned int doc = 0; doc < _nDoc; ++doc) {
		for (unsigned int i = 0; i < _docSums[doc]; ++i) {
		    if (topics[doc][i] != snodes[offset + i]) {
			throwLogicError("Bad GraphView in LD constructor");
		    }
		}
		offset += _docSums[doc];
	    }
	}

	void LDA::rebuildTable()
	{
	    for (unsigned int t = 0; t < _nTopic; ++t) {
		for (unsigned int w = 0; w < _nWord; ++w) {
		    _topicsByWord[w][t] = 0;
		}
	    }

	    vector<StochasticNode *> const &s = _gv->nodes();
	    unsigned int offset = 0;
	    for (unsigned int d = 0; d < _nDoc; ++d) {
		for (unsigned int i = 0; i < _docSums[d]; ++i) {
		    int topic = _topicTokens[d][i];
		    int word = static_cast<int>(*s[offset + i]->value(_chain)) 
			- 1;
		    _topicsByWord[word][topic]++;
		}
		offset += _docSums[d];
	    }
	}

	void LDA::update(RNG *rng)
	{
	    if (!_wordsObserved) {
		rebuildTable();
	    }

	    double wordHyperSum = 
		accumulate(_wordHyper, _wordHyper + _nWord, 0.0);
	    
	    vector<double> sump(_nTopic);
	    for (unsigned int doc = 0; doc < _nDoc; ++doc) {

		vector<int> &thisDocTopics = _topicsByDoc[doc];

		for (unsigned int i = 0; i < _docSums[doc]; ++i) {

		    //Find current topic and word values
		    int &topic = _topicTokens[doc][i];
		    int const &word = _wordTokens[doc][i];

		    vector<int> &thisWordTopics = _topicsByWord[word];

		    //Remove current value from tables
		    thisDocTopics[topic]--;
		    thisWordTopics[topic]--;
		    _topicSums[topic]--;

		    //Calculate probability vector
		    vector<double> prob(_nTopic);
		    for (unsigned int t = 0; t < _nTopic; ++t) {
			double prior = thisDocTopics[t] + _topicHyper[t];
			double likelihood = 
			    (thisWordTopics[t] + _wordHyper[word]) /
			    (_topicSums[t] + wordHyperSum);
			prob[t] = prior * likelihood;
		    }
		    
		    //Draw random sample from categorical distribution
		    partial_sum(prob.begin(), prob.end(), sump.begin());
		    double p = rng->uniform() * sump.back();
		    topic = upper_bound(sump.begin(), sump.end(), p) -
			sump.begin();
		    if (topic == _nTopic) --topic;

		    //Restore current value to tables
		    thisDocTopics[topic]++;
		    thisWordTopics[topic]++;
		    _topicSums[topic]++;
		}
	    }
	    
	    vector<double> value;
	    value.reserve(_gv->length());
	    for (unsigned int d = 0; d < _nDoc; ++d) {
		for (unsigned int i = 0; i < _docSums[d]; ++i) {
		    value.push_back(_topicTokens[d][i] + 1);
		}
	    }
	    _gv->setValue(value, _chain);
	}
	
	bool LDA::canSample(vector<vector<StochasticNode*> > const &topics,
			    vector<vector<StochasticNode*> > const &words,
			    vector<StochasticNode*> const &topic_priors,
			    vector<StochasticNode*> const &word_priors,
			    Graph const &graph)
	{
	    //Set up dimensions
	    unsigned int nDoc = topics.size();
	    if (nDoc == 0) return false;
	    if (words.size() != nDoc) return false;
	    if (topic_priors.size() != nDoc) return false;

	    unsigned int nTopic = word_priors.size();
	    if (nTopic == 0) return false;
	    
	    unsigned int nWord = word_priors[0]->length();
	    if (nWord == 0) return false;

	    //Topics and words are arranged hierarchically by document

	    //Topics and words have categorical distribution. Topics
	    //must be unobserved (because we are sampling them) and
	    //words must be observed (this may be relaxed in a future version).
	    for (unsigned int d = 0; d < nDoc; ++d) {
		if (topics[d].size() != words[d].size()) return false;
		for (unsigned int i = 0; i < topics[d].size(); ++i) {
		    if (!isCat(topics[d][i])) return false;
		    if (!isCat(words[d][i])) return false;
		    if (isObserved(topics[d][i])) return false;
		    if (isParameter(words[d][i])) return false;
		}
	    }
	    
	    //Word priors have a Dirichlet distribution.  They are fully
	    //unobserved, as we need to marginalize over their
	    //distribution
	    for (unsigned int t = 0; t < nTopic; ++t) {
		if (!isDirichlet(word_priors[t])) return false;
		if (isObserved(word_priors[t])) return false;
		if (word_priors[t]->length() != nWord) return false;
	    }
	    //Word priors have a common hyper-prior
	    Node const *word_hyper = word_priors[0]->parents()[0];
	    for (unsigned int t = 1; t < nTopic; ++t) {
		if (word_priors[t]->parents()[0] != word_hyper) return false;
	    }

	    //Topic priors have a Dirichlet distribution. They are fully
	    //unobserved, as we need to marginalize over their
	    //distribution.
	    for (unsigned int d = 0; d < nDoc; ++d) {
		if (!isDirichlet(topic_priors[d])) return false;
		if (isObserved(topic_priors[d])) return false;
		if (topic_priors[d]->length() != nTopic) return false;
	    }
	    //Topic priors have a common hyper-prior
	    Node const *topic_hyper = topic_priors[0]->parents()[0];
	    for (unsigned int d = 1; d < nDoc; ++d) {
		if (topic_priors[d]->parents()[0] != topic_hyper) return false;
	    }

	    //topics[d] is the stochastic children of topic_priors[d].
	    //There are no intermediate deterministic nodes.
	    for (unsigned int d = 0; d < nDoc; ++d) {
		
		SingletonGraphView tp_gv(topic_priors[d], graph);

		if (!tp_gv.deterministicChildren().empty()) return false;

		set<StochasticNode *> set1, set2;
		set1.insert(tp_gv.stochasticChildren().begin(),
			    tp_gv.stochasticChildren().end());
		set2.insert(topics[d].begin(), topics[d].end());
		if (set1 != set2) return false;
	    }

	    //The stochastic children of word_priors are in words[0]
	    //... words[nDoc] 
	    
	    GraphView wp_gv(word_priors, graph);

	    set<StochasticNode*> wp_schildren;
	    wp_schildren.insert(wp_gv.stochasticChildren().begin(),
				wp_gv.stochasticChildren().end());
	    set<StochasticNode*> all_words;
	    for (unsigned int d = 0; d < nDoc; ++d) {
		all_words.insert(words[d].begin(), words[d].end());
	    }
	    if (wp_schildren != all_words) return false;

	    //Create set to contain all mixture node children of word_priors
	    set<MixtureNode const*> wp_mixnodes;

	    //The intermediate deterministic nodes are all mixture
	    //nodes with a common mixTab.
	    vector<DeterministicNode*> const &wp_dchild = 
		wp_gv.deterministicChildren();
	    MixTab const *mtab = 0;
	    for (unsigned int i = 0; i < wp_dchild.size(); ++i) {
		MixtureNode const *m = asMixture(wp_dchild[i]);
		if (m == 0) return false;
		if (i == 0) {
		    mtab = m->mixTab();
		}
		else if (m->mixTab() != mtab) {
		    return false;
		}
		wp_mixnodes.insert(m);
	    }

	    //The only possible values of the mixTab are in word_priors
	    set<Node const*> wp_set;
#ifndef _RWSTD_NO_MEMBER_TEMPLATES
	    wp_set.insert(word_priors.begin(), word_priors.end());
#else
	    for (vector<StochasticNode*>::const_iterator p = 
		   word_priors.begin(); p != word_priors.end(); ++p)
	      {
		wp_set.insert(*p);
	      }
#endif
	    for (RangeIterator p(mtab->range()); !p.atEnd(); p.nextLeft()) {
		if (wp_set.count(mtab->getNode(p)) == 0) return false;
	    }

	    //Create set to contain all mixture node children of topics
	    set<MixtureNode const*> topic_mixnodes;

	    //Each topic is the index of a single mixture node and its
	    //single stochastic child is the corresponding word. All
	    //mixture nodes have the same mixTab
	    for (unsigned int d = 0; d < nDoc; ++d) {
		for (unsigned int i = 0; i < topics[d].size(); ++i) {
		    
		    SingletonGraphView topic_gv(topics[d][i], graph);

		    vector<DeterministicNode*> topic_dchild = 
			topic_gv.deterministicChildren();
		    if (topic_dchild.size() != 1) return false;
		    MixtureNode const *m = asMixture(topic_dchild[0]);
		    if (m->index_size() != 1) return false;
		    if (m->parents()[0] != topics[d][i]) return false;
		    if (m->mixTab() != mtab) return false;

		    topic_mixnodes.insert(m);

		    vector<StochasticNode*> topic_schild = 
			topic_gv.stochasticChildren();
		    if (topic_schild.size() != 1) return false;
		    if (topic_schild[0] != words[d][i]) return false;
		}
	    }

	    //Finally, check that the deterministic children of topics
	    //and of word_priors are the same
	    if (wp_mixnodes != topic_mixnodes) return false;

	    return true;
	}
    }
}

