#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

/* Helper functions for linear interpolation used by InterpLin, InterpLin2D,
   and InterpLin3D */

#include <cstddef>

namespace jags {
    namespace bugs {

	/*
	 * Structure to hold coordinates of a point within a grid
	 * l, u are the indices of the bracketing values in the grid.
	 * t is the fractional distance in [0,1] between lower and upper bounds
	 */
	struct GridIndex {
	    size_t l;
	    size_t u;
	    double t;
	};

	/*
	 * Reruns the coordinates of val within a grid of length n, which is
	 * assumed to be increasing.
	 */
	 GridIndex find_grid_index(double const *grid, unsigned long n, double val);

	 /*
	  * Checks whether array x of length n is increasing
	  */
	bool increasing(double const *x, unsigned long n);

	/*
	 * Convenience function to do linear interpolation
	 * l, u are lower and upper bounds
	 * t is the fractional distance in [0,1] between l and u
	 */
	inline double interplin(double t, double l, double u) {
	    return (1 - t) * l + t * u;
	}

    }
}

#endif /* INTERPOLATION_H_ */

