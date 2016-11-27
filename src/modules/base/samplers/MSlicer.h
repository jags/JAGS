#ifndef M_SLICER_H_
#define M_SLICER_H_

#include <sampler/MutableSampleMethod.h>
#include <vector>

namespace jags {

    class StochasticNode;
    class SingletonGraphView;

    namespace base {

	/**
	 * Slice sampler for real-valued distributions
	 */
	class MSlicer : public MutableSampleMethod 
	{
	    SingletonGraphView const *_gv;
	    unsigned int _chain;
	    unsigned int _length;
	    std::vector<double> _width;
	    double _max;
	    std::vector<double> _value;
	    bool _adapt;
	    unsigned int _iter;
	    std::vector<double> _sumdiff;

	    double update0(RNG *rng, unsigned int i,
			   std::vector<double> const &lower,
			   std::vector<double> const &upper);
	    void update1(RNG *rng, std::vector<double> const &lower,
			 std::vector<double> const &upper);
	    double const * value() const;
	    void setValue(double value, unsigned int i);
	    void setValue(std::vector<double> const &value);
	    double logDensity() const;
	  public:
	    MSlicer(SingletonGraphView const *gv, unsigned int chain,
		    double width = 1, long maxwidth = 10);
	    void update(RNG *rng);
	    static bool canSample(StochasticNode const *node);
	    bool isAdaptive() const;
	    void adaptOff();
	    bool checkAdaptation() const;
	};

    }
}

#endif /* M_SLICER_H_ */

