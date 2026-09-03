#include <config.h>

#include "InterpLin2D.h"
#include "interpolation.h"

#include <util/dim.h>

using std::vector;

namespace jags {
    namespace bugs {

	static double bilinear_interp(double const *p,
				      double const *x_grid, unsigned long Nx, 
				      double const *y_grid, unsigned long Ny,
				      double const *f_matrix)
	{  
	    GridIndex gx = find_grid_index(x_grid, Nx, p[0]);
	    GridIndex gy = find_grid_index(y_grid, Ny, p[1]);
	    
	    // extract values from the outcome matrix
	    auto F = [&f_matrix, &Nx](unsigned long ix,
				      unsigned long iy)
	    {
		return f_matrix[ix + Nx * iy];
	    };
	    
	    double f00 = F(gx.l, gy.l);
	    double f10 = F(gx.u, gy.l);
	    double f01 = F(gx.l, gy.u);
	    double f11 = F(gx.u, gy.u);
    
	    // Interpolate over x
	    double f0 = interplin(gx.t, f00, f10);
	    double f1 = interplin(gx.t, f01, f11);

	    // Interpolate over y
	    double f = interplin(gy.t, f0, f1);

	    return f;
	}


	InterpLin2D::InterpLin2D()
	    : ArrayFunction("interp.lin2d", 4)
	{
	}

	void 
	InterpLin2D::evaluate (double *value, vector<double const *> const &args,
			       vector<vector<unsigned long>> const &dims) const
	{
	    *value =  bilinear_interp(args[0],
				      args[1], dims[1][0],
				      args[2], dims[2][0],
				      args[3]);
	}

	vector<unsigned long> 
	InterpLin2D::dim(vector <vector<unsigned long>> const &dims,
			 vector<double const *> const &) const
	{
	    return  vector<unsigned long>(1,1UL);
	}
	
	bool 
	InterpLin2D::checkParameterDim(vector<vector<unsigned long>> const &dims) const
	{
	    // Check that the coordinate parameter is a vector of length 2
	    if (dims[0].size() != 1 || dims[0][0] != 2UL) {
		return false;
	    }
	    
	    // Check that grid parameters are vectors and the value
	    // parameter is a matrix. NB This enforces grid lengths >=
	    // 2 which is required by find_grid_index.
	    if (!isVector(dims[1]) || !isVector(dims[2]) || !isMatrix(dims[3])) {
		return false;
	    }

	    // Check that the lengths of the grid parameters conform
	    // to the value dimensions
	    if (dims[3][0] != dims[1][0] || dims[3][1] != dims[2][0]) {
		return false;
	    }

	    return true;
	}
	
	bool 
	InterpLin2D::checkParameterValue(vector<double const *> const &args,
					 vector<vector<unsigned long>> const &dims) const
	{
	    if (!increasing(args[1], dims[1][0])) return false;
	    if (!increasing(args[2], dims[2][0])) return false;
	    return true;
	}

    }
}
