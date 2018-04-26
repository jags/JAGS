#include <config.h>
#include <sampler/Slicer.h>
#include <rng/RNG.h>
#include <util/nainf.h>

#include <cmath>
#include <cfloat>
#include <stdexcept>

//Minimum length of adaptive phase before we adjust width
#define MIN_ADAPT 50

using std::vector;
using std::fabs;
using std::runtime_error;

namespace jags {

Slicer::Slicer(double width, unsigned int max)
    : _width(width), _adapt(true), _max(max), _sumdiff(0), _iter(0), 
      _state(SLICER_OK)
{
}

bool Slicer::updateStep(RNG *rng)
{
    // Test current value
    double g0 = logDensity();
    if (!jags_finite(g0)) {
	if (g0 > 0) {
	    _state = SLICER_POSINF;
	    return false;
	    //throw runtime_error("Singularity in likelihood found by Slicer");
	}
	else {
	    _state = SLICER_NEGINF;
	    return false;
	    //FIXME: Not very informative
	    //throw runtime_error("Error in Slicer: Current value is inconsistent with data");
	}
    }

    // Generate auxiliary variable
    double z = g0 - rng->exponential();;

    // Generate random interval of width "_width" about current value
    double xold = value();
    double L = xold - rng->uniform() * _width; 
    double R = L + _width;

    double lower = JAGS_NEGINF, upper = JAGS_POSINF;
    getLimits(&lower, &upper);

    // Stepping out 

    // Randomly set number of steps in left and right directions,
    // subject to the limit in the maximal size of the interval
    int j = static_cast<int>(rng->uniform() * (_max + 1));
    int k = static_cast<int>(_max) - j;

    if (L < lower) {
	L = lower;
    }
    else {
	setValue(L);
	while (j-- > 0 && logDensity() > z) {
	    L -= _width;
	    if (L < lower) {
		L = lower;
		break;
	    }
	    setValue(L);
	}
    }

    if (R > upper) {
	R = upper;
    }
    else {
	setValue(R);
	while (k-- > 0 && logDensity() > z) {
	    R += _width;
	    if (R > upper) {
		R = upper;
		break;
	    }
	    setValue(R);
	}
    }

    // Keep sampling from the interval until acceptance (the loop is
    // guaranteed to terminate).
    double xnew;
    for(;;) {
	xnew =  L + rng->uniform() * (R - L);
	setValue(xnew);
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

    if (_adapt) {
	_sumdiff += _iter * fabs(xnew - xold);
	++_iter;
	if (_iter > MIN_ADAPT) {
	    _width = 2 * _sumdiff / _iter / (_iter - 1);  
	}
    }

    return true;
}

bool Slicer::updateDouble(RNG *rng)
{
  using namespace std;

  // Test current value
  double g0 = logDensity();
  if (!jags_finite(g0)) {
    if (g0 < 0) {
	_state = SLICER_NEGINF;
	return false;
	//throw runtime_error("Error in Slicer: Current value is inconsistent with data");
    }
    else {
	_state = SLICER_POSINF;
	return false;
	//throw runtime_error("Singularity in likelihood found by Slicer");
    }
  }

  // Generate auxiliary variable
  double z = g0 - rng->exponential();

  // Generate random interval of width "_width" about current value
  double xold = value();
  double L = xold - rng->uniform() * _width; 
  double R = L + _width;

  double lower = JAGS_NEGINF, upper = JAGS_POSINF;
  getLimits(&lower, &upper);

  // Doubling 
  bool left_ok = false, right_ok = false;
  for (unsigned int i = 0; i < _max; ++i) {
    if (rng->uniform() < 0.5) {
      if (L >= lower) {
	L = 2*L - R;
        if (L < lower) {
           left_ok = true;
        }
        else {
	    setValue(L);
	    left_ok = logDensity() < z;
        }
      }
      else {
	left_ok = true;
      }
    }
    else {
      if (R <= upper) {
	R = 2*R - L;
        if (R > upper) {
           right_ok = true;
        }
        else {
	    setValue(R);
	    right_ok = logDensity() < z;
        }
      }
      else {
	right_ok = true;
      }
    }
    if (left_ok && right_ok)
      break;
  }
  
  // Keep sampling from the interval until acceptance (the loop is
  // guaranteed to terminate).
  double Lbar = L, Rbar = R;
  double xnew;
  for(;;) {
    xnew =  Lbar + rng->uniform() * (Rbar - Lbar);
    if (xnew >= lower && xnew <= upper) {
	setValue(xnew);
	double g = logDensity();
	if (g >= z && accept(xold, xnew, z, L, R, lower, upper)) {
	    // The accept function will alter the current value. So we
	    // must reset it.
	    setValue(xnew);
	    break;
	}
    }
    // shrink the interval
    if (xnew <= xold) {
	Lbar = xnew;
    }
    else {
      Rbar = xnew;
    }
  }

  if (_adapt) {
    _sumdiff += _iter * fabs(xnew - xold);
    ++_iter;
    if (_iter > MIN_ADAPT) {
      _width = 2 * _sumdiff / _iter / (_iter - 1);  
    }
  }
  return true;
}

bool Slicer::accept(double xold, double xnew, double z, double L, double R,
                    double lower, double upper)
{
  //Acceptance step for doubling update method

  bool d = false;
  while ((R - L) > 1.1 * _width) {
    double M = (L + R)/2;
    if ((xold < M && xnew >= M) || (xold >= M && xnew < M))
      d = true;
    if (xnew < M) {
      R = M;
    }
    else {
      L = M;
    }
    if (d) {
      bool right_ok = true;
      if (R <= upper) {
	  setValue(R);
	  right_ok = logDensity() < z;
      }
      bool left_ok = true;
      if (L >= lower) {
	  setValue(L);
         left_ok = logDensity() < z;
      }
      if (left_ok && right_ok) {
	return false;
      }
    }
  }
  return true;
}  

void Slicer::adaptOff()
{
  _adapt = false;
}

bool Slicer::checkAdaptation() const
{
    //FIXME We could try a bit harder than this
    return (_iter > MIN_ADAPT);
}

bool Slicer::isAdaptive() const
{
    return true;
}

SlicerState Slicer::state() const
{
    return _state;
}

} //namespace jags
