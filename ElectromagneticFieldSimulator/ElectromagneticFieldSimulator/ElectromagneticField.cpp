#include "ElectromagneticField.h"
#include "BoundaryConditions.h"
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <Eigen/SVD>

ElectromagneticField::ElectromagneticField(int numNodes_)
	: numNodes(numNodes_),
	globalStiffnessMatrix(Eigen::MatrixXd::Zero(numNodes_, numNodes_)),
	electricField(Eigen::VectorXd::Zero(numNodes_)),
	magneticField(Eigen::VectorXd::Zero(numNodes_)),
	globalLoadVector(Eigen::VectorXd::Zero(numNodes_)),
	numElements(0) {}

// Generate a structured 2D triangular mesh
void ElectromagneticField::generateMesh(double width, double height, int nx, int ny) {
	nodes.clear();
	elements.clear();

	double dx = width / (nx - 1);
	double dy = height / (ny - 1);

	for (int j = 0; j < ny; ++j) {
		for (int i = 0; i < nx; ++i) {
			nodes.push_back({ i * dx, j * dy });
		}
	}

	for (int j = 0; j < ny - 1; ++j) {
		for (int i = 0; i < nx - 1; ++i) {
			int n0 = j * nx + i;
			int n1 = n0 + 1;
			int n2 = n0 + nx;
			int n3 = n2 + 1;

			elements.push_back({ n0, n1, n2 });
			elements.push_back({ n1, n3, n2 });
		}
	}

	numElements = elements.size();
}

// Assemble the global stiffness matrix
void ElectromagneticField::assembleSystemMatrix() {
	for (int elem = 0; elem < numElements; ++elem) {
		int n0 = elements[elem][0];
		int n1 = elements[elem][1];
		int n2 = elements[elem][2];

		Eigen::Vector2d p0(nodes[n0][0], nodes[n0][1]);
		Eigen::Vector2d p1(nodes[n1][0], nodes[n1][1]);
		Eigen::Vector2d p2(nodes[n2][0], nodes[n2][1]);

		double area = 0.5 * std::abs((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
		if (std::abs(area) < 1e-10) continue;

		Eigen::Vector2d gradN0, gradN1, gradN2;
		gradN0 << p1[1] - p2[1], p2[0] - p1[0];
		gradN1 << p2[1] - p0[1], p0[0] - p2[0];
		gradN2 << p0[1] - p1[1], p1[0] - p0[0];

		Eigen::Matrix3d localStiffnessMatrix;
		localStiffnessMatrix << gradN0.dot(gradN0), gradN0.dot(gradN1), gradN0.dot(gradN2),
			gradN1.dot(gradN0), gradN1.dot(gradN1), gradN1.dot(gradN2),
			gradN2.dot(gradN0), gradN2.dot(gradN1), gradN2.dot(gradN2);
		localStiffnessMatrix *= (1.0 / (4.0 * area));

		globalStiffnessMatrix(n0, n0) += localStiffnessMatrix(0, 0);
		globalStiffnessMatrix(n0, n1) += localStiffnessMatrix(0, 1);
		globalStiffnessMatrix(n0, n2) += localStiffnessMatrix(0, 2);
		globalStiffnessMatrix(n1, n0) += localStiffnessMatrix(1, 0);
		globalStiffnessMatrix(n1, n1) += localStiffnessMatrix(1, 1);
		globalStiffnessMatrix(n1, n2) += localStiffnessMatrix(1, 2);
		globalStiffnessMatrix(n2, n0) += localStiffnessMatrix(2, 0);
		globalStiffnessMatrix(n2, n1) += localStiffnessMatrix(2, 1);
		globalStiffnessMatrix(n2, n2) += localStiffnessMatrix(2, 2);
	}

	globalStiffnessMatrix += Eigen::MatrixXd::Identity(numNodes, numNodes) * 1e-6;
}

// Compute the curl of the electric field
Eigen::VectorXd ElectromagneticField::computeCurlE() const {
	Eigen::VectorXd curlE = Eigen::VectorXd::Zero(numNodes);
	int nx = sqrt(numNodes);
	int ny = nx;

	double dx = 1.0 / (nx - 1);
	double dy = 1.0 / (ny - 1);

	for (int j = 1; j < ny - 1; ++j) {
		for (int i = 1; i < nx - 1; ++i) {
			int idx = j * nx + i;

			double dEy_dx = (electricField[idx + 1] - electricField[idx - 1]) / (2 * dx);
			double dEx_dy = (electricField[idx + nx] - electricField[idx - nx]) / (2 * dy);

			curlE(idx) = dEy_dx - dEx_dy;
		}
	}

	return curlE;
}

// Compute the curl of the magnetic field
Eigen::VectorXd ElectromagneticField::computeCurlB() const {
	Eigen::VectorXd curlB = Eigen::VectorXd::Zero(numNodes);
	int nx = sqrt(numNodes);
	int ny = nx;

	double dx = 1.0 / (nx - 1);
	double dy = 1.0 / (ny - 1);

	for (int j = 1; j < ny - 1; ++j) {
		for (int i = 1; i < nx - 1; ++i) {
			int idx = j * nx + i;

			double dBz_dy = (magneticField[idx + nx] - magneticField[idx - nx]) / (2 * dy);
			double dBz_dx = (magneticField[idx + 1] - magneticField[idx - 1]) / (2 * dx);

			curlB(idx) = dBz_dy - dBz_dx;
		}
	}

	return curlB;
}

// Solve the system for both electric and magnetic fields
void ElectromagneticField::solveFields(double dt, const Eigen::VectorXd& globalLoad, BoundaryConditions& bc) {
	bc.applyDirichlet(*this);

	if (globalLoad.size() != numNodes) {
		std::cerr << "Error: Load vector size does not match the number of nodes." << std::endl;
		return;
	}

	electricField = globalStiffnessMatrix.ldlt().solve(globalLoad);

	bc.applyDirichlet(*this);

	Eigen::VectorXd curlE = computeCurlE();

	for (int i = 0; i < numNodes; ++i) {
		magneticField(i) = magneticField(i) - dt * curlE(i);
	}

	bc.applyDirichlet(*this);

	Eigen::VectorXd curlB = computeCurlB();

	for (int i = 0; i < numNodes; ++i) {
		electricField(i) = electricField(i) + dt * curlB(i);
	}

	bc.applyDirichlet(*this);
}

// Save results to a CSV file for visualization
void ElectromagneticField::saveResultsToFile(const std::string& filename) const {
	std::ofstream file(filename);
	if (file.is_open()) {
		file << "Time,Node,X,Y,Ex,Ey,Bz\n";
		for (int i = 0; i < numNodes; ++i) {
			file << i << "," << nodes[i][0] << "," << nodes[i][1] << ","
				<< electricField(i) << "," << electricField(i) << "," << magneticField(i) << "\n";
		}
		file.close();
	}
}
