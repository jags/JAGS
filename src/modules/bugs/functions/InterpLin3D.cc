#include <config.h>

#include "InterpLin3D.h"
#include "interpolation.h"

#include <util/dim.h>

//debuggin
#include <iostream>

using std::vector;

namespace jags {
    namespace bugs {

	static double trilinear_interp(
	    double const *p,
	    double const *x_grid, unsigned long Nx, 
	    double const *y_grid, unsigned long Ny,
	    double const *z_grid, unsigned long Nz,
	    double const *f_array)
	{    
	    GridIndex gx = find_grid_index(x_grid, Nx, p[0]);
	    GridIndex gy = find_grid_index(y_grid, Ny, p[1]);
	    GridIndex gz = find_grid_index(z_grid, Nz, p[2]);

	    /* extract corners of bounding cuboid */

	    auto F = [&f_array, &Nx, &Ny](unsigned long ix,
					  unsigned long iy,
					  unsigned long iz)
	    {
		return f_array[ix + Nx * (iy + Ny * iz)];
	    };

	    double f000 = F(gx.l, gy.l, gz.l);
	    double f100 = F(gx.u, gy.l, gz.l);
	    double f010 = F(gx.l, gy.u, gz.l);
	    double f110 = F(gx.u, gy.u, gz.l);
    
	    double f001 = F(gx.l, gy.l, gz.u);
	    double f101 = F(gx.u, gy.l, gz.u);
	    double f011 = F(gx.l, gy.u, gz.u);
	    double f111 = F(gx.u, gy.u, gz.u);

	    /* trilinear interpolation between corners */

	    // Interpolate over x
	    double f00 = interplin(gx.t, f000, f100);
	    double f10 = interplin(gx.t, f010, f110);
	    double f01 = interplin(gx.t, f001, f101);
	    double f11 = interplin(gx.t, f011, f111);

	    // Interpolate over y
	    double f0 = interplin(gy.t, f00, f10);
	    double f1 = interplin(gy.t, f01, f11);

	    // Interpolate over z
	    double f = interplin(gz.t, f0, f1);
	    std::cout << f << std::endl;
	    return f;
	}

	InterpLin3D::InterpLin3D()
	    : ArrayFunction("interp.lin3d", 5)
	{
	}

	void 
	InterpLin3D::evaluate (double *value, vector<double const *> const &args,
			       vector<vector<unsigned long>> const &dims) const
	{
	    *value =  trilinear_interp(args[0],
				       args[1], dims[1][0],
				       args[2], dims[2][0],
				       args[3], dims[3][0],
				       args[4]);
	}

	vector<unsigned long> 
	InterpLin3D::dim (vector<vector<unsigned long>> const &dims,
			  vector<double const *> const &) const
	{
	    return  vector<unsigned long>(1,1UL);
	}

	bool 
	InterpLin3D::checkParameterDim (vector<vector<unsigned long>> const &dims) const
	{
	    // Check that the coordinate parameter is a vector of length 3
	    if (dims[0].size() != 1 || dims[0][0] != 3UL) {
		return false;
	    }

	    // Check that grid parameters are vectors
	    // NB This enforces grid lengths >= 2 which is required by find_grid_index.
	    if (!isVector(dims[1]) || !isVector(dims[2]) || !isVector(dims[3]) || !isArray(dims[4])) {
		return false;
	    }

	    // Check that value array has 3 dimensions
	    if (dims[4].size() != 3) {
		return false;
	    }
	    
	    // Check that the lengths of the grid parameters conform with the value dimensions
	    if (dims[4][0] != dims[1][0] || dims[4][1] != dims[2][0] || dims[4][2] != dims[3][0]) {
		return false;
	    }

	    return true;
	}

	bool 
	InterpLin3D::checkParameterValue(vector<double const *> const &args,
					 vector<vector<unsigned long>> const &dims) const
	{
	    // Check that the grids are strictly increasing
	    if (!increasing(args[1], dims[1][0])) return false;
	    if (!increasing(args[2], dims[2][0])) return false;
	    if (!increasing(args[3], dims[3][0])) return false;
	    return true;
	}


    }
}
