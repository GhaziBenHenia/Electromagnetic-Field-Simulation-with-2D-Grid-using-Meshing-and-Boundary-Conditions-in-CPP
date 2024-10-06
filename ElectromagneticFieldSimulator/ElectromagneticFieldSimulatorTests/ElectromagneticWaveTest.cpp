#include <gtest/gtest.h>
#include "ElectromagneticWave.h"
#include "BoundaryConditions.h"
#include <fstream>  
#include <Eigen/Dense>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Test that wave propagation initializes correctly
TEST(ElectromagneticWaveTest, Initialization) {
	// ElectromagneticWave object with a 3x3 grid (9 nodes, 8 elements)
	ElectromagneticWave wave(9);

	// Generate the mesh
	wave.generateMesh(1.0, 1.0, 3, 3);

	// Check that the mesh has the correct number of nodes and elements
	EXPECT_EQ(wave.getNodes().size(), 9);       // 3x3 grid has 9 nodes
	EXPECT_EQ(wave.getElements().size(), 8);    // 4 cells with 2 triangles each = 8 triangles
}

// Test boundary conditions during wave propagation with custom load
TEST(ElectromagneticWaveTest, AbsorbingBoundaryConditionsWithCustomLoad) {
	// Create an ElectromagneticWave object with a 3x3 grid (9 nodes, 8 elements)
	ElectromagneticWave wave(9);

	wave.generateMesh(1.0, 1.0, 3, 3);
	wave.assembleSystemMatrix();

	// Custom global load vector
	Eigen::VectorXd customLoad = Eigen::VectorXd::Constant(9, 1.0);
	wave.setGlobalLoadVector(customLoad);

	BoundaryConditions bc;

	wave.applyAbsorbingBoundaryConditions();

	// Check that the boundary nodes have small or zero field values
	Eigen::VectorXd electricField = wave.getElectricField();
	Eigen::VectorXd magneticField = wave.getMagneticField();

	int nx = sqrt(9);  // 3x3 grid

	for (int i = 0; i < nx; ++i) {
		EXPECT_NEAR(electricField(i), 0.0, 1e-5);  // Boundary electric field should be close to zero
		EXPECT_NEAR(magneticField(i), 0.0, 1e-5);  // Boundary magnetic field should be close to zero
		EXPECT_NEAR(electricField(9 - 1 - i), 0.0, 1e-5);  // Boundary electric field should be close to zero
		EXPECT_NEAR(magneticField(9 - 1 - i), 0.0, 1e-5);  // Boundary magnetic field should be close to zero
	}
}

TEST(ElectromagneticWaveTest, WavePropagationVisualizationWithCustomLoad) {
	// Parameters
	int numNodes = 100;
	double dt = 0.01;      // Time step for wave propagation
	int timeSteps = 10;

	// ElectromagneticWave object with a 10x10 grid (100 nodes, 180 elements)
	ElectromagneticWave wave(numNodes);

	wave.generateMesh(1.0, 1.0, 10, 10);

	wave.assembleSystemMatrix();

	// Custom global load vector with spatial variation
	Eigen::VectorXd customLoad(numNodes);
	int gridSize = 10;  // 10x10 grid
	for (int i = 0; i < numNodes; ++i) {
		int x = i % gridSize;
		int y = i / gridSize;

		// Example: a sinusoidal load based on the grid position
		customLoad(i) = std::sin(M_PI * x / (gridSize - 1)) * std::sin(M_PI * y / (gridSize - 1));
	}
	wave.setGlobalLoadVector(customLoad);

	// boundary conditions
	BoundaryConditions bc;

	std::ofstream resultFile("wave_propagation.csv");
	resultFile << "Time,Node,X,Y,ElectricField,MagneticField\n";

	// Simulate wave propagation using the dedicated method
	for (int step = 0; step < timeSteps; ++step) {
		double currentTime = step * dt;
		wave.simulateWavePropagation(dt, 1, bc);

		// Output results for this time step
		const Eigen::VectorXd& electricField = wave.getElectricField();
		const Eigen::VectorXd& magneticField = wave.getMagneticField();
		const auto& nodes = wave.getNodes();

		for (int node = 0; node < numNodes; ++node) {
			resultFile << currentTime << "," << node << "," << nodes[node][0] << "," << nodes[node][1] << ","
				<< electricField(node) << "," << magneticField(node) << "\n";
		}
	}
	resultFile.close();

	// Check that the boundary nodes have small or reasonable field values
	Eigen::VectorXd electricFieldCheck = wave.getElectricField();
	Eigen::VectorXd magneticFieldCheck = wave.getMagneticField();

	int nx = sqrt(numNodes);  // 10x10 grid

	for (int i = 0; i < nx; ++i) {
		EXPECT_NEAR(electricFieldCheck(i), 0.0, 1e-3);  // Boundary electric field should be reasonably small
		EXPECT_NEAR(magneticFieldCheck(i), 0.0, 1e-3);  // Boundary magnetic field should be reasonably small
		EXPECT_NEAR(electricFieldCheck(numNodes - 1 - i), 0.0, 1e-3);  // Boundary electric field should be reasonably small
		EXPECT_NEAR(magneticFieldCheck(numNodes - 1 - i), 0.0, 1e-3);  // Boundary magnetic field should be reasonably small
	}
}
