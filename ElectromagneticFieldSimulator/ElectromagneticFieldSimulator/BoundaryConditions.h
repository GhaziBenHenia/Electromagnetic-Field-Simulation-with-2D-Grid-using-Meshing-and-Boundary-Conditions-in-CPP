#pragma once
#include "ElectromagneticField.h"

class BoundaryConditions {
public:
	// Constructor
	BoundaryConditions();

	// Apply Dirichlet boundary conditions
	void applyDirichlet(ElectromagneticField& emField);

	// Apply Neumann boundary conditions
	//void applyNeumann(ElectromagneticField& emField);
};
