#pragma once

#include "Derivatives.h"
#include "LinearSystem.h"
#include "ScalarField.h"
#include "SGridBasedData.h"

#include <numbers>

namespace PhysX::SurfaceTension {

template <int Dim>
double getCurvature(const Vector<Dim, double> &pos, const ScalarField<Dim> &sdf, double dx)
{
	DECLARE_DIM_TYPES(Dim)

	double kappa = 0;
	double invDx = 1. / dx;
	for (int axis = 0; axis < Dim; axis++)
		kappa += sdf.gradient(pos + VectorDd::Unit(axis) * dx).normalized()[axis] - sdf.gradient(pos - VectorDd::Unit(axis) * dx).normalized()[axis];
	kappa *= invDx * .5;
	return std::abs(kappa) > invDx ? (kappa < 0 ? -1 : 1) * invDx : kappa;
}

template <int Dim>
void solve(
	SGridBasedData<Dim, double> &velocity,
	const ScalarField<Dim> &sdf,
	const SGridBasedData<Dim, double> &boundaryFraction,
	const double coeff, const double maxSteps,
	const double dt)
{
	DECLARE_DIM_TYPES(Dim)

	const double dx = velocity.grids[0].spacing;
	const double invDx = velocity.grids[0].invSpacing;
	double eps = maxSteps * dx;
	int n = 0;

	SGridBasedData<Dim, int> grid2mat(velocity.grids, -VectorDi::Ones());
	std::vector<std::pair<int, int>> mat2grid;
	std::vector<Tripletd> elements;

	const auto curvature = [&](const VectorDd &pos)->double {
		return getCurvature(pos, sdf, dx);
	};
	const auto dirac = [&](const double phi)->double {
		if (phi < -eps) return 0;
		else if (phi > eps) return 0;
		else return 0.5 * (1.0 + std::cos(std::numbers::pi * phi / eps)) / eps;
	};

	// Find interface faces.
	forEach(velocity.grids, [&](const int axis, const VectorDi &face) {
		const VectorDd pos = velocity.grids[axis].position(face);
		const double phi = sdf(pos);
		if (-eps < phi && phi < eps && boundaryFraction[axis][face] < 1) {
			grid2mat[axis][face] = n++;
			mat2grid.push_back(std::pair(axis, int(velocity.grids[axis].index(face))));
		}
	});

	SparseMatrixd matLaplacian(n, n);
	VectorXd vel(n);
	VectorXd rhs(n);

	for (int r = 0; r < n; r++) {
		const int axis = mat2grid[r].first;
		const VectorDi face = velocity.grids[axis].coordinate(mat2grid[r].second);
		const VectorDd pos = velocity.grids[axis].position(face);
		const VectorDd normal = sdf.gradient(pos).normalized();
		const double kappa = curvature(pos);
		// eps = std::min(std::max(1. / kappa, 3. * dx), maxSteps * dx);
		const double delta = std::max(dirac(sdf(pos)), 1e-7 * invDx);
		const VectorDd jacobian = Derivatives::getFirst<2>(velocity[axis], face);
		const MatrixDd hessian = Derivatives::getSecond<2>(velocity[axis], face);
		double diaCoef = 1 / delta;
		// vel[r] = rhs[r] = velocity[axis][face] / delta - coeff * dt * (kappa * normal[axis]);
		vel[r] = rhs[r] = velocity[axis][face] / delta - coeff * dt * (kappa * normal[axis] + dt * ((normal.transpose() * hessian * normal).value() + kappa * jacobian.dot(normal)));
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbFace = Grid<Dim>::neighbor(face, i);
			if (!velocity.grids[axis].isValid(nbFace) || boundaryFraction[axis][nbFace] == 1) continue;
			const double weight = coeff * dt * dt / (dx * dx);
			diaCoef += weight;
			if (int c = grid2mat[axis][nbFace]; c != -1)
				elements.push_back(Tripletd(r, c, -weight));
			else rhs[r] += weight * velocity[axis][nbFace];
		}
		elements.push_back(Tripletd(r, r, diaCoef));
	}
	matLaplacian.setFromTriplets(elements.begin(), elements.end());

	LinearSystem::solve(matLaplacian, vel, rhs, 2000, 1e-5);
	std::cout << " (Capillary)";

	for (int r = 0; r < n; r++) {
		const int axis = mat2grid[r].first;
		const VectorDi face = velocity.grids[axis].coordinate(mat2grid[r].second);
		velocity[axis][face] = vel[r];
	}
}

}
