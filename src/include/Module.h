#ifndef _MODULE_H_
#define _MODULE_H_

#include <vector>

class InverseLinkFunc;
class Function;
class Distribution;
class SamplerFactory;
class RNGFactory;

class Module {
    std::vector<Function*> _functions;
    std::vector<Distribution*> _distributions;
    std::vector<SamplerFactory*> _sampler_factories;
    std::vector<RNGFactory*> _rng_factories;
public:
    Module();
    virtual ~Module();
    void insert(InverseLinkFunc*);
    void insert(Function*);
    void insert(Distribution*);
    void insert(SamplerFactory*);
    void insert(RNGFactory*);
    std::vector<Function*> const &functions() const;
    std::vector<Distribution*> const &distributions() const;
    std::vector<SamplerFactory*> const &samplerFactories() const;
    std::vector<RNGFactory*> const &rngFactories() const;
};

#endif /* _MODULE_H_ */
