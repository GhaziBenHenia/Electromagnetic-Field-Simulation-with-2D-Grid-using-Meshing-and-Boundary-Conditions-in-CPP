#pragma once
#include <vector>

class Grid2D {
public:
	int nx, ny;              // Number of grid points in x and y directions
	double dx, dy;           // Spacing between grid points
	std::vector<std::vector<double>> Ex;  // Electric field (x-component)
	std::vector<std::vector<double>> Ey;  // Electric field (y-component)
	std::vector<std::vector<double>> Bz;  // Magnetic field (z-component for 2D)

	Grid2D(int nx_, int ny_, double dx_, double dy_);

	// Function to initialize grid values
	void initializeGrid();

	// Function to display grid info
	void displayGridInfo() const;
};
