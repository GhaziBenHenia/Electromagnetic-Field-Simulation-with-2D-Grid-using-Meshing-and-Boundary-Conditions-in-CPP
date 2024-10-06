#include "Grid2D.h"
#include <iostream>

Grid2D::Grid2D(int nx_, int ny_, double dx_, double dy_)
	: nx(nx_), ny(ny_), dx(dx_), dy(dy_),
	Ex(nx_, std::vector<double>(ny_, 0.0)),
	Ey(nx_, std::vector<double>(ny_, 0.0)),
	Bz(nx_, std::vector<double>(ny_, 0.0)) {}

// Initialize all fields to zero
void Grid2D::initializeGrid() {
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			Ex[i][j] = 0.0;
			Ey[i][j] = 0.0;
			Bz[i][j] = 0.0;
		}
	}
}

// Display grid information
void Grid2D::displayGridInfo() const {
	std::cout << "Grid: " << nx << "x" << ny << " with spacing dx=" << dx << ", dy=" << dy << std::endl;
}
