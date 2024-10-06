#pragma once
#include "ElectromagneticField.h"

class BoundaryConditions;


class ElectromagneticWave : public ElectromagneticField {
public:

	ElectromagneticWave(int numNodes_);

	// Method for simulating wave propagation
	void simulateWavePropagation(double timeStep, int steps, BoundaryConditions& bc);

	void applyAbsorbingBoundaryConditions();
};
