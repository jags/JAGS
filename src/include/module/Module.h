#ifndef MODULE_H_
#define MODULE_H_

#include <vector>
#include <string>
#include <list>
#include <utility>

#include <function/FunctionPtr.h>
#include <distribution/DistPtr.h>

namespace jags {

class SamplerFactory;
class RNGFactory;
class MonitorFactory;
class LinkFunctiontion;
class ScalarFunction;
class VectorFunction;
class ArrayFunction;
class ScalarDist;
class RScalarDist;
class VectorDist;
class ArrayDist;

/**
 * @short Memory management for dynamically loadable modules
 *
 * Dynamically loadable modules may extend the capabilities of the JAGS library
 * by defining new subclasses of Function, Distribution, SamplerFactory, 
 * and RNGFactory. 
 *
 * Each module must instantiate a subclass of Module. The constructor for
 * this subcluass must dynamically allocate instances of the sub-classes 
 * defined by the module, and store them with the appropriate insert 
 * member function.
 */
class Module {
    std::string _name;
    bool _loaded;
    std::vector<FunctionPtr> _fp_list;
    std::vector<Function*> _functions;
    std::vector<std::pair<DistPtr, FunctionPtr> > _obs_functions;
    std::vector<DistPtr> _dp_list;
    std::vector<Distribution*> _distributions;
    std::vector<SamplerFactory*> _sampler_factories;
    std::vector<RNGFactory*> _rng_factories;
    std::vector<MonitorFactory*> _monitor_factories;
public:
    Module(std::string const &name);
    virtual ~Module();

    void insert(ScalarFunction*);
    void insert(LinkFunction*);
    void insert(VectorFunction*);
    void insert(ArrayFunction*);
    void insert(Distribution*);

    void insert(RScalarDist*);
    void insert(ScalarDist*);
    void insert(VectorDist*);
    void insert(ArrayDist*);

    void insert(ScalarDist*, ScalarFunction*);
    void insert(ScalarDist*, LinkFunction*);
    void insert(ScalarDist*, VectorFunction*);
    void insert(ScalarDist*, ArrayFunction*);

    void insert(VectorDist*, ScalarFunction*);
    void insert(VectorDist*, LinkFunction*);
    void insert(VectorDist*, VectorFunction*);
    void insert(VectorDist*, ArrayFunction*);

    void insert(ArrayDist*, ScalarFunction*);
    void insert(ArrayDist*, LinkFunction*);
    void insert(ArrayDist*, VectorFunction*);
    void insert(ArrayDist*, ArrayFunction*);

    void insert(SamplerFactory*);
    void insert(RNGFactory*);
    void insert(MonitorFactory*);

    std::vector<Function*> const &functions() const;
    std::vector<Distribution*> const &distributions() const;
    std::vector<SamplerFactory*> const &samplerFactories() const;
    std::vector<RNGFactory*> const &rngFactories() const;
    std::vector<MonitorFactory*> const &monitorFactories() const;
    
    void load();
    void unload();
    std::string const &name() const;
    static std::list<Module *> &modules();
    static std::list<Module *> &loadedModules();
};

} /* namespace jags */

#endif /* MODULE_H_ */
