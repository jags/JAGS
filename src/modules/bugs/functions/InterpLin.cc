#include <config.h>

#include "InterpLin.h"
#include "interpolation.h"

//#include <util/nainf.h>
//#include <util/dim.h>

using std::vector;

namespace jags {
    namespace bugs {

	static double linear_interp(const double *f_vector,
				    double const *x_grid, unsigned long Nx,
				    double x)
	{  
	    GridIndex ix = find_grid_index(x_grid, Nx, x);
	    
	    double f0 = f_vector[ix.l];
	    double f1 = f_vector[ix.u];
	    
	    // Interpolate over x
	    return interplin(ix.t, f0, f1);
	}
	

	InterpLin::InterpLin() : ScalarVectorFunction("interp.lin", 3)
	{}
    
	double InterpLin::scalarEval(vector<double const *> const &args,
				     vector<unsigned long> const &lengths) const
	{
	    double x = args[0][0];
	    double const *xgrid = args[1];
	    double const *f_vector = args[2];

	    unsigned long Nx = lengths[1];

	    return linear_interp(f_vector, xgrid, Nx, x);
	}
	
	bool InterpLin::checkParameterLength(vector<unsigned long> const &lengths)
	    const
	{
	    return lengths[0] == 1 && lengths[1] > 1 && lengths[2] == lengths[1];
	}
	
	bool 
	InterpLin::checkParameterValue(vector <double const *> const &args,
				       vector <unsigned long> const &lengths) const
	{
	    return increasing(args[1], lengths[1]);
	}
	
    }
}
