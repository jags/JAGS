#include <config.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <util/nainf.h>

#include "MSlicer.h"

#include <cmath>
#include <cfloat>

using std::vector;
using std::string;

#define MIN_ADAPT 50

#include <iostream>

namespace jags {
    namespace base {

	MSlicer::MSlicer(SingletonGraphView const *gv, unsigned int chain,
			 double width, long maxwidth)
	    : _gv(gv), _chain(chain), _length(gv->length()),
	      _width(_length, width), _max(maxwidth),
	      _value(_length), _adapt(true), _iter(0), _sumdiff(_length)
	{
	    if (!canSample(gv->node())) {
		throwLogicError("Invalid MSlicer");
	    }
	    gv->checkFinite(chain);
	    gv->getValue(_value, chain);
	}
	
	bool MSlicer::canSample(StochasticNode const *node)
	{
	    if (node->isDiscreteValued() || node->length() <= 1)
		return false;

	    if (node->df() < node->length())
		return false; 

	    return true;
	}

	double MSlicer::update0(RNG *rng, unsigned int i,
				vector<double> const &lower,
				vector<double> const &upper)
	{
	    // Generate auxiliary variable
	    double g0 = logDensity();
	    double z = g0 - rng->exponential();;

	    // Generate random interval of width "_width[i]" about current value
	    double xold = _value[i];
	    double L = xold - rng->uniform() * _width[i]; 
	    double R = L + _width[i];

	    // Stepping out 

	    // Randomly set number of steps in left and right directions,
	    // subject to the limit in the maximal size of the interval
	    int j = static_cast<int>(rng->uniform() * _max);
	    int k = _max - 1 - j;


	    if (L < lower[i]) {
		L = lower[i];
	    }
	    else {
		setValue(L, i);
		while (j-- > 0 && logDensity() > z) {
		    L -= _width[i];
		    if (L < lower[i]) {
			L = lower[i];
			break;
		    }
		    setValue(L, i);
		}
	    }

	    if (R > upper[i]) {
		R = upper[i];
	    }
	    else {
		setValue(R, i);
		while (k-- > 0 && logDensity() > z) {
		    R += _width[i];
		    if (R > upper[i]) {
			R = upper[i];
			break;
		    }
		    setValue(R, i);
		}
	    }

	    // Keep sampling from the interval until acceptance (the loop is
	    // guaranteed to terminate).
	    double xnew;
	    for(;;) {
		xnew =  L + rng->uniform() * (R - L);
		setValue(xnew, i);
		double g = logDensity();
		if (g >= z - DBL_EPSILON) {
		    // Accept point
		    break;
		}
		else {
		    // shrink the interval
		    if (xnew < xold) {
			L = xnew;
		    }
		    else {
			R = xnew;
		    }
		}
	    }

	    return xnew;
	}

	void MSlicer::update1(RNG *rng, vector<double> const &lower,
			      vector<double> const &upper)
	{
	    // Generate auxiliary variable
	    double g0 = logDensity();
	    double z = g0 - rng->exponential();;

	    // Generate random interval of width "_width" about current value
	    vector<double> L(_length), R(_length);
	    for (unsigned int i = 0; i < _length; ++i) {
		L[i] = _value[i] - rng->uniform() * 2 * _width[i]; 
		R[i] = L[i] + 2 * _width[i];
	    }

	    // Keep sampling from the interval until acceptance (the loop is
	    // guaranteed to terminate).
	    vector<double> xold(_value);
	    vector<double> xnew(_length);
	    for(;;) {
		for (unsigned int i = 0; i < _length; ++i) {
		    xnew[i] =  L[i] + rng->uniform() * (R[i] - L[i]);
		}
		setValue(xnew);
		double g = logDensity();
		if (g >= z - DBL_EPSILON) {
		    // Accept point
		    break;
		}
		else {
		    // shrink the interval
		    for (unsigned int i = 0; i < _length; ++i) {
			if (xnew[i] < xold[i]) {
			    L[i] = xnew[i];
			}
			else {
			    R[i] = xnew[i];
			}
		    }
		}
	    }
	}

	void MSlicer::update(RNG *rng)
	{
	    // Test current value
	    double g0 = logDensity();
	    if (!jags_finite(g0)) {
		if (g0 > 0) {
		    throwNodeError(_gv->node(),
				   "Slicer stuck at value with infinite density");
		}
		else {
		    throwNodeError(_gv->node(),
				   "Current value is inconsistent with data");
		}
	    }

	    vector<double> lower(_length), upper(_length);
	    _gv->node()->support(&lower[0], &upper[0], _length, _chain);

	    if (_adapt) {
		//During adaptation we do univariate slices

		++_iter;
		for (unsigned int i = 0; i < _length; ++i) {
		    double xold = _value[i];
		    double xnew = update0(rng, i, lower, upper);
		    double delta = fabs(xnew - xold) - _width[i];
		    _width[i] += 2 * delta / (_iter+1);
		}
	    }
	    /*
	    else {
				
		for (unsigned int i = 0; i < _length; ++i) {
		    update0(rng, i, lower, upper);
		}

	    }
	    */

	    update1(rng, lower, upper);
	}

	double const * MSlicer::value() const
	{
	    return _gv->node()->value(_chain);
	}
 
	void MSlicer::setValue(double value, unsigned int i)
	{
	    _value[i] = value; 
	    _gv->setValue(&_value[0], _value.size(), _chain);
	}

	void MSlicer::setValue(vector<double> const &value)
	{
	    _value = value; 
	    _gv->setValue(&_value[0], _value.size(), _chain);
	}


	double MSlicer::logDensity() const
	{
	    return _gv->logFullConditional(_chain);
	}

	bool MSlicer::isAdaptive() const
	{
	    return true;
	}

	void MSlicer::adaptOff() {
	    _adapt = false;
	    std::cout << "chain " << _chain << " node " << _gv->node() <<
		std::endl;
	    for (unsigned int i = 0; i < _length; ++i) {
		std::cout << _width[i] << " ";
	    }
	    std::cout << std::endl;
	}

	bool MSlicer::checkAdaptation() const
	{
	    return _iter > MIN_ADAPT;
	}
	
    }
}
