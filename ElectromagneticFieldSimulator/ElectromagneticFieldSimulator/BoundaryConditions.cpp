#include "BoundaryConditions.h"
#include "ElectromagneticField.h"
#include <iostream>

BoundaryConditions::BoundaryConditions() {}

// Apply Dirichlet boundary conditions: set electric and magnetic fields to zero at the boundary nodes
void BoundaryConditions::applyDirichlet(ElectromagneticField& emField) {
	const std::vector<std::vector<double>>& nodes = emField.getNodes();
	Eigen::MatrixXd& stiffnessMatrix = emField.getGlobalStiffnessMatrix();
	Eigen::VectorXd& loadVector = emField.getGlobalLoadVector();

	int nx = sqrt(nodes.size());

	for (int i = 0; i < nodes.size(); ++i) {
		int x = i % nx;
		int y = i / nx;

		// Apply Dirichlet conditions at the boundary nodes (top, bottom, left, and right edges)
		if (x == 0 || x == nx - 1 || y == 0 || y == nx - 1) {
			stiffnessMatrix.row(i).setZero();
			stiffnessMatrix.col(i).setZero();
			stiffnessMatrix(i, i) = 1.0;
			loadVector(i) = 0.0;
		}
	}
}


// Apply Neumann boundary conditions: set field derivatives to zero
//void BoundaryConditions::applyNeumann(ElectromagneticField& emField) {
//    std::cout << "Neumann boundary condition applied (not yet implemented)." << std::endl;
//}
