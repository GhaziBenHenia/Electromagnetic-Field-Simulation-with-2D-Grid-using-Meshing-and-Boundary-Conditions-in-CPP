#include "ElectromagneticWave.h"
#include <cmath>
#include <iostream>
#include "BoundaryConditions.h"

ElectromagneticWave::ElectromagneticWave(int numNodes_)
	: ElectromagneticField(numNodes_) {}

// Method to simulate wave propagation
void ElectromagneticWave::simulateWavePropagation(double timeStep, int steps, BoundaryConditions& bc) {
	for (int step = 0; step < steps; ++step) {
		std::cout << "Time Step " << step << std::endl;

		// Apply boundary conditions
		bc.applyDirichlet(*this);

		// Compute the curl of the electric field
		Eigen::VectorXd curlE = computeCurlE();

		// Update magnetic field using Faraday's Law: dB/dt = -curl(E)
		magneticField -= timeStep * curlE;

		// Apply absorbing boundary conditions to the magnetic field
		applyAbsorbingBoundaryConditions();

		// Compute the curl of the magnetic field
		Eigen::VectorXd curlB = computeCurlB();

		// Update electric field using Ampère's Law: dE/dt = curl(B)
		electricField += timeStep * curlB;

		// Apply absorbing boundary conditions to the electric field
		applyAbsorbingBoundaryConditions();
	}

	std::cout << "Wave propagation simulation completed." << std::endl;
}

// Apply absorbing boundary conditions to reduce reflections from the boundaries
void ElectromagneticWave::applyAbsorbingBoundaryConditions() {
	int nx = sqrt(numNodes);

	// Apply conditions at the boundary to prevent wave reflection
	for (int i = 0; i < nx; ++i) {
		// Top and bottom boundaries
		electricField(i) = 0.0;  // Set to zero to fully absorb
		magneticField(i) = 0.0;  // Set to zero to fully absorb
		electricField(numNodes - 1 - i) = 0.0;  // Set to zero to fully absorb
		magneticField(numNodes - 1 - i) = 0.0;  // Set to zero to fully absorb
	}

	// Apply absorbing boundary conditions for the sides
	for (int i = 0; i < nx; ++i) {
		// Left and right boundaries
		electricField(i * nx) = 0.0;  // Set to zero to fully absorb
		magneticField(i * nx) = 0.0;  // Set to zero to fully absorb
		electricField((i + 1) * nx - 1) = 0.0;  // Set to zero to fully absorb
		magneticField((i + 1) * nx - 1) = 0.0;  // Set to zero to fully absorb
	}
}
