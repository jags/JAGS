#ifndef NORM_MIX_METHOD_H_
#define NORM_MIX_METHOD_H_

#include <sampler/TemperedMetropolis.h>
#include <graph/StochasticNode.h>

#include <vector>

namespace jags {

class Graph;
class GraphView;

namespace mix {

    struct DirichletInfo;

/**
 * @short Sample method for normal mixture models
 *
 * NormMix implements a Metropolis-Hastings method for normal
 * mixture models with tempered transitions as recommended for mixture
 * models by Celeux et al (Celeux G, Hurn, M, Robert CP, Computational
 * and inferential difficulties with mixture posterior distributions,
 * J Am Statist Ass, 95, 957-970, 2000).
 */
    class NormMix : public TemperedMetropolis
    {
	GraphView const *_gv;
	unsigned int _chain;
	double *_lower;
	double *_upper;
	std::vector<DirichletInfo*> _di;

    public:
	/**
	 * Constructor for Mix Method with tempered transitions
	 *
	 * @param max_level Number of levels between minimum and maximum
	 * temperature
	 *
	 * @param max_temp Maximum temperature
	 *
	 * @param nrep Number of Metropolis-Hastings updates to do at
	 * each level
	 */
	NormMix(GraphView const *gv, unsigned int chain,
		unsigned int max_level, double max_temp, unsigned int nrep);
	~NormMix();
	void getValue(std::vector<double> &value) const;
	void setValue(std::vector<double> const &value);
	double logPrior() const;
	double logLikelihood() const;
	double logJacobian(std::vector<double> const &value) const;
	void step(std::vector<double> &value, double step, RNG *rng) const;
	static bool canSample(std::vector<StochasticNode *> const &snodes);
    };

}}

#endif /* NORM_MIX_METHOD_H_ */
