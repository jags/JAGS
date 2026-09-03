#include "interpolation.h"

#include <algorithm>
#include <stdexcept>

using std::logic_error;
using std::lower_bound;
using std::distance;

namespace jags {
    namespace bugs {
	
	GridIndex find_grid_index(double const *grid, unsigned long n,
				  double val)
	{
	    if (n < 2)
		throw logic_error("Grid must contain at least 2 points");

	    // locate first element >= val
	    auto it = lower_bound(grid, grid + n, val);
	    size_t u = distance(grid, it);
	    if (u == 0) {
		return {0, 1, 0.0};
	    }
	    else if (u == n) {
		return {n-2, n-1, 1.0};
	    }
	    else {
		size_t l = u - 1;
		double t = (val - grid[l]) / (grid[u] - grid[l]);
		return {l, u, t};
	    }
	}

	bool increasing(double const *x, unsigned long n)
	{
	    for (unsigned long i = 1; i < n; ++i) {
		if (x[i] <= x[i-1]) return false;
	    }
	    return true;
	}

    }
}

