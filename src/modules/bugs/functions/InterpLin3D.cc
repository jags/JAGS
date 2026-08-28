#include <config.h>

#include "InterpLin3D.h"
#include "interpolation.h"

#include <util/dim.h>

using std::vector;

namespace jags {
    namespace bugs {

	static double trilinear_interp(
	    double const *f_array,
	    double const *x_grid, unsigned long Nx, 
	    double const *y_grid, unsigned long Ny,
	    double const *z_grid, unsigned long Nz,
	    double x, double y, double z)
	{    
	    GridIndex ix = find_grid_index(x_grid, Nx, x);
	    GridIndex iy = find_grid_index(y_grid, Ny, y);
	    GridIndex iz = find_grid_index(z_grid, Nz, z);

	    // -----------------------------------------
	    // extract corners
	    // -----------------------------------------

	    auto F = [&](unsigned long ix,
			 unsigned long iy,
			 unsigned long iz)
	    {
		return f_array[ix + Nx * (iy + Ny * iz)];
	    };

	    double f000 = F(ix.l, iy.l, iz.l);
	    double f100 = F(ix.u, iy.l, iz.l);
	    double f010 = F(ix.l, iy.u, iz.l);
	    double f110 = F(ix.u, iy.u, iz.l);
    
	    double f001 = F(ix.l, iy.l, iz.u);
	    double f101 = F(ix.u, iy.l, iz.u);
	    double f011 = F(ix.l, iy.u, iz.u);
	    double f111 = F(ix.u, iy.u, iz.u);

	    // trilinear interpolation

	    // Interpolate over x
	    double f00 = interplin(ix.t, f000, f100);
	    double f10 = interplin(ix.t, f010, f110);
	    double f01 = interplin(ix.t, f001, f101);
	    double f11 = interplin(ix.t, f011, f111);

	    // Interpolate over y
	    double f0 = interplin(iy.t, f00, f10);
	    double f1 = interplin(iy.t, f01, f11);

	    // Interpolate over z
	    double f = interplin(iz.t, f0, f1);

	    return f;
	}

	InterpLin3D::InterpLin3D()
	    : ArrayFunction("interp.lin2d", 5)
	{
	}

	void 
	InterpLin3D::evaluate (double *value, vector<double const *> const &args,
			       vector<vector<unsigned long>> const &dims) const
	{
	    double x = args[0][0];
	    double y = args[0][1];
	    double z = args[0][2];

	    unsigned long Nx = dims[1][0];
	    unsigned long Ny = dims[2][0];
	    unsigned long Nz = dims[3][0];

	    double const *x_grid = args[1];
	    double const *y_grid = args[2];
	    double const *z_grid = args[3];
	    double const *f_matrix = args[4];
	    
	    
	    *value =  trilinear_interp(f_matrix,
				       x_grid, Nx, 
				       y_grid, Ny,
				       z_grid, Nz,
				       x, y, z);
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
	    if (!isVector(dims[0]) || dims[0].size() != 3) {
		return false;
	    }

	    // Check that grid parameters are vectors and value parameter is an array
	    // NB This enforces grid lengths >= 2 which is required by find_grid_index.
	    if (!isVector(dims[1]) || !isVector(dims[2]) || !isVector(dims[3]) || !isArray(dims[4])) {
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
