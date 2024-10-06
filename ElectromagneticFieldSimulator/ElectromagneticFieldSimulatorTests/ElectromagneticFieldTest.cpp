#include <gtest/gtest.h>
#include "ElectromagneticField.h"
#include "BoundaryConditions.h"
#include <fstream>

// Test that the mesh is generated correctly
TEST(ElectromagneticFieldTest, MeshGeneration) {
	// ElectromagneticField object with a 5x5 grid (25 nodes)
	ElectromagneticField emField(25);

	// A 1x1 grid with 5x5 nodes (nx = 5, ny = 5)
	emField.generateMesh(1.0, 1.0, 5, 5);

	// Check that the mesh has the correct number of nodes and elements
	EXPECT_EQ(emField.getNodes().size(), 25);       // 5x5 grid has 25 nodes
	EXPECT_EQ(emField.getElements().size(), 32);    // 16 cells with 2 triangles each = 32 triangles

	// Check some node positions
	EXPECT_NEAR(emField.getNodes()[0][0], 0.0, 1e-6);   // x-coordinate of node 0
	EXPECT_NEAR(emField.getNodes()[0][1], 0.0, 1e-6);   // y-coordinate of node 0
	EXPECT_NEAR(emField.getNodes()[24][0], 1.0, 1e-6);  // x-coordinate of node 24
	EXPECT_NEAR(emField.getNodes()[24][1], 1.0, 1e-6);  // y-coordinate of node 24
}

// Test the assembly and solution of both the electric and magnetic fields with boundary conditions
TEST(ElectromagneticFieldTest, SolveFieldsWithBoundaryConditions) {
	// ElectromagneticField object with a 3x3 grid (9 nodes)
	ElectromagneticField emField(9);
	BoundaryConditions bc;

	emField.generateMesh(1.0, 1.0, 3, 3);
	emField.assembleSystemMatrix();

	// Apply Dirichlet boundary conditions
	bc.applyDirichlet(emField);

	// Time step for the simulation
	double dt = 1e-9;

	// Set a custom load vector
	Eigen::VectorXd customLoad(9);
	customLoad.setZero();
	customLoad(0) = 1.0;  // Apply a load to node 0

	// Solving the fields with the user-specified time step and external load vector
	emField.solveFields(dt, customLoad, bc);

	EXPECT_EQ(emField.getElectricField().size(), 9);  // Should match the number of nodes
	EXPECT_EQ(emField.getMagneticField().size(), 9);  // Should match the number of nodes

	EXPECT_NEAR(emField.getElectricField()(0), 1.0, 1e-6);  // Boundary node 0
	EXPECT_NEAR(emField.getElectricField()(8), 0.0, 1e-6);  // Boundary node 8
	EXPECT_NEAR(emField.getMagneticField()(0), 0.0, 1e-6);  // Boundary node 0
	EXPECT_NEAR(emField.getMagneticField()(8), 0.0, 1e-6);  // Boundary node 8
}

// Test to debug magnetic field instability
TEST(ElectromagneticFieldTest, MagneticFieldInstabilityDebug) {
	ElectromagneticField emField(9);

	emField.generateMesh(1.0, 1.0, 3, 3);
	emField.assembleSystemMatrix();
	BoundaryConditions bc;

	Eigen::VectorXd customLoadVector(9);
	customLoadVector.setZero();
	customLoadVector(0) = 1.0;

	double dt = 1e-12;

	std::cout << "Starting Magnetic Field: \n" << emField.getMagneticField() << std::endl;

	for (int step = 0; step < 10; ++step) {
		emField.solveFields(dt, customLoadVector, bc);
		std::cout << "Time Step " << step << " - Electric Field: \n" << emField.getElectricField() << std::endl;
		std::cout << "Time Step " << step << " - Magnetic Field: \n" << emField.getMagneticField() << std::endl;
	}

	// Check final field values
	EXPECT_NEAR(emField.getElectricField()(0), 1.0, 1e-6);
	EXPECT_NEAR(emField.getElectricField()(8), 0.0, 1e-6);
	EXPECT_NEAR(emField.getMagneticField()(0), 0.0, 1e-6);
	EXPECT_NEAR(emField.getMagneticField()(8), 0.0, 1e-6);
}

TEST(ElectromagneticFieldTest, SimulationWithDifferentValues) {
	// Scenario 1: 3x3 grid with load at node 0
	{
		ElectromagneticField emField(9);  // 3x3 grid (9 nodes)
		BoundaryConditions bc;

		emField.generateMesh(1.0, 1.0, 3, 3);
		emField.assembleSystemMatrix();

		bc.applyDirichlet(emField);

		Eigen::VectorXd loadVector(9);
		loadVector.setZero();
		loadVector(0) = 1.0;

		double dt = 1e-9;

		for (int step = 0; step < 10; ++step) {
			emField.solveFields(dt, loadVector, bc);

			std::cout << "Scenario 1 - Step " << step << ": Electric Field: \n" << emField.getElectricField() << std::endl;
			std::cout << "Scenario 1 - Step " << step << ": Magnetic Field: \n" << emField.getMagneticField() << std::endl;

			EXPECT_NEAR(emField.getElectricField()(0), 1.0, 1e-6);  // Boundary condition at node 0
			EXPECT_NEAR(emField.getMagneticField()(0), 0.0, 1e-6);  // Boundary condition for magnetic field
		}
	}

	// Scenario 2: 5x5 grid with load applied at the center node
	{
		ElectromagneticField emField(25);  // 5x5 grid (25 nodes, 40 elements)
		BoundaryConditions bc;

		emField.generateMesh(1.0, 1.0, 5, 5);
		emField.assembleSystemMatrix();

		bc.applyDirichlet(emField);

		Eigen::VectorXd loadVector(25);
		loadVector.setZero();
		loadVector(12) = 1.0;

		double dt = 1e-9;

		for (int step = 0; step < 10; ++step) {
			emField.solveFields(dt, loadVector, bc);

			std::cout << "Scenario 2 - Step " << step << ": Electric Field: \n" << emField.getElectricField() << std::endl;
			std::cout << "Scenario 2 - Step " << step << ": Magnetic Field: \n" << emField.getMagneticField() << std::endl;

			EXPECT_GT(emField.getElectricField()(12), 0.0);  // Check that it is greater than zero
			EXPECT_LT(emField.getElectricField()(12), 1.0);  // It should be less than the applied load

			EXPECT_NEAR(emField.getElectricField()(0), 0.0, 1e-6);   // Boundary condition at node 0
			EXPECT_NEAR(emField.getMagneticField()(12), 0.0, 1e-6);  // Magnetic field at the center should evolve
		}
	}

	// Scenario 3: 4x4 grid with loads applied at the corner nodes
	{
		ElectromagneticField emField(16);  // 4x4 grid (16 nodes)
		BoundaryConditions bc;

		emField.generateMesh(1.0, 1.0, 4, 4);  // Generate 4x4 mesh
		emField.assembleSystemMatrix();

		bc.applyDirichlet(emField);

		// Apply a load at the four corner nodes
		Eigen::VectorXd loadVector(16);
		loadVector.setZero();
		loadVector(0) = 1.0;   // Load at top-left corner
		loadVector(3) = 1.0;   // Load at top-right corner
		loadVector(12) = 1.0;  // Load at bottom-left corner
		loadVector(15) = 1.0;  // Load at bottom-right corner

		double dt = 1e-8;

		for (int step = 0; step < 10; ++step) {
			emField.solveFields(dt, loadVector, bc);

			std::cout << "Scenario 3 - Step " << step << ": Electric Field: \n" << emField.getElectricField() << std::endl;
			std::cout << "Scenario 3 - Step " << step << ": Magnetic Field: \n" << emField.getMagneticField() << std::endl;

			EXPECT_GT(emField.getElectricField()(0), 0.0);  // Load at node 0 (top-left)
			EXPECT_GT(emField.getElectricField()(3), 0.0);  // Load at node 3 (top-right)
			EXPECT_GT(emField.getElectricField()(12), 0.0); // Load at node 12 (bottom-left)
			EXPECT_GT(emField.getElectricField()(15), 0.0); // Load at node 15 (bottom-right)

			EXPECT_EQ(emField.getElectricField()(5), 0.0);  // Middle node should have no load
		}
	}

}
