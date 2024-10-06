#pragma once
#include <vector>
#include <Eigen/Dense>

class BoundaryConditions;

class ElectromagneticField {
protected:
	int numNodes;
	int numElements;

	Eigen::MatrixXd globalStiffnessMatrix;
	Eigen::VectorXd electricField;
	Eigen::VectorXd magneticField;
	Eigen::VectorXd globalLoadVector;

	std::vector<std::vector<double>> nodes;
	std::vector<std::vector<int>> elements;

	// Helper methods for computing the curl of the electric and magnetic fields
	Eigen::VectorXd computeCurlE() const;
	Eigen::VectorXd computeCurlB() const;

public:

	ElectromagneticField(int numNodes_);

	// Method to generate a structured 2D triangular mesh
	void generateMesh(double width, double height, int nx, int ny);

	// Method to assemble the global system matrix
	void assembleSystemMatrix();

	// Method to solve the system for both electric and magnetic fields
	void solveFields(double dt, const Eigen::VectorXd& globalLoad, BoundaryConditions& bc);

	// Save results to a CSV file for visualization
	void saveResultsToFile(const std::string& filename) const;

	// Getter methods
	Eigen::MatrixXd& getGlobalStiffnessMatrix() { return globalStiffnessMatrix; }
	const Eigen::VectorXd& getElectricField() const { return electricField; }
	const Eigen::VectorXd& getMagneticField() const { return magneticField; }
	const std::vector<std::vector<double>>& getNodes() const { return nodes; }
	const std::vector<std::vector<int>>& getElements() const { return elements; }
	void setGlobalLoadVector(const Eigen::VectorXd& load) { globalLoadVector = load; }
	Eigen::VectorXd& getGlobalLoadVector() { return globalLoadVector; }
};
