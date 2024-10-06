#include <gtest/gtest.h>
#include "ElectromagneticField.h"
#include "BoundaryConditions.h"
#include <fstream>

TEST(ElectromagneticFieldTest, VisualizationTest) {
	// Set up electromagnetic field simulation
	ElectromagneticField emField(25);  // Example grid
	BoundaryConditions bc;
	emField.generateMesh(1.0, 1.0, 5, 5);
	emField.assembleSystemMatrix();

	// Apply boundary conditions and load
	bc.applyDirichlet(emField);
	Eigen::VectorXd loadVector(25);
	loadVector.setZero();
	loadVector(12) = 5;
	double dt = 1e-9;

	// Solve fields
	for (int step = 0; step < 10; ++step) {
		emField.solveFields(dt, loadVector, bc);
	}

	// Save results to CSV for visualization
	emField.saveResultsToFile("./output_field.csv");


}

TEST(ElectromagneticFieldTest, LargeGridComplexLoad) {
	ElectromagneticField emField(100);  // 10x10 grid has 100 nodes
	BoundaryConditions bc;
	emField.generateMesh(1.0, 1.0, 10, 10);  // Generate a 1x1 grid with 10x10 nodes
	emField.assembleSystemMatrix();

	bc.applyDirichlet(emField);

	// Complex load distribution across multiple nodes
	Eigen::VectorXd loadVector(100);
	loadVector.setZero();
	loadVector(45) = 1.0;   // Apply a load at a central node
	loadVector(55) = 0.5;   // Apply a smaller load at another node
	loadVector(66) = 0.3;   // Apply a small load at another node
	loadVector(33) = -0.5;  // Apply a negative load to a node

	// Define a time step for the simulation
	double dt = 1e-9;

	// Run the simulation for more steps (e.g., 20 steps)
	for (int step = 0; step < 20; ++step) {
		emField.solveFields(dt, loadVector, bc);
	}

	// Save results to CSV for visualization
	emField.saveResultsToFile("./large_grid_output.csv");

}