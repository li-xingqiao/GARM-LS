#pragma once

#include "Derivatives.h"
#include "FastMarching.h"

namespace PhysX::Reinitialization {

template <typename Scheme, int Dim>
inline void solve(GridBasedData<Dim, double> &phi, const int maxSteps)
{
	Scheme(phi.grid).perform(phi, maxSteps);
}

template <typename Scheme, int Dim>
inline void solve(GridBasedData<Dim, double> &phi, GridBasedData<Dim, Vector<Dim, double>> &psi, const int maxSteps, const double cfl = 0.9)
{
	Scheme(phi.grid).perform(phi, psi, maxSteps, cfl);
}

template <int Dim, int Len>
void advance(const GridBasedData<Dim, double> &field, GridBasedData<Dim, double> &newField, const double dt, std::pair<double, double>(*getDrv)(const std::array<double, Len> &, const double))
{
	DECLARE_DIM_TYPES(Dim)

	const auto square = [](const double x) { return x * x; };
	const auto sign = [](const double val, const double sp) { return val / std::sqrt(val * val + sp * sp); };

	parallelForEach(field.grid, [&](const VectorDi &coord) {
		std::pair<double, double> drv[Dim];
		for (int axis = 0; axis < Dim; axis++)
			drv[axis] = getDrv(field.stencil<Len - 1>(coord - VectorDi::Unit(axis) * (Len >> 1), axis), field.grid.invSpacing);
		double sum = 0;
		if (field[coord] < 0) {
			for (int axis = 0; axis < Dim; axis++) {
				sum += square(std::min(drv[axis].first, 0.0));
				sum += square(std::max(drv[axis].second, 0.0));
			}
		}
		else {
			for (int axis = 0; axis < Dim; axis++) {
				sum += square(std::max(drv[axis].first, 0.0));
				sum += square(std::min(drv[axis].second, 0.0));
			}
		}
		sum = std::sqrt(sum);
		newField[coord] = field[coord] + dt * sign(field[coord], field.grid.spacing * sum) * (1 - sum);
	});
}

template <int Order, int Dim>
inline void solve(GridBasedData<Dim, double> &phi, const int maxSteps, const double cfl = .5)
{
	DECLARE_DIM_TYPES(Dim)
	static_assert(Order == 1 || Order == 3 || Order == 5, "Order must be 1, 3 or 5.");

	const double dt = phi.grid.spacing * cfl;
	for (int iter = 0; iter < maxSteps; iter++) {
		if constexpr (Order == 1) {
			GridBasedData<Dim, double> newPhi(phi.grid);
			advance(phi, newPhi, dt, Derivatives::upwind1);
			phi = newPhi;
		}
		else if constexpr (Order == 3) {
			GridBasedData<Dim, double> newPhi1(phi.grid);
			GridBasedData<Dim, double> newPhi2(phi.grid);
			advance(phi, newPhi1, dt, Derivatives::eno3);
			advance(newPhi1, newPhi2, dt, Derivatives::eno3);
			parallelForEach(phi.grid, [&](const VectorDi &coord) { phi[coord] = (phi[coord] + newPhi2[coord]) / 2; });
		}
		else {
			GridBasedData<Dim, double> newPhi1(phi.grid);
			GridBasedData<Dim, double> newPhi2(phi.grid);
			advance(phi, newPhi1, dt, Derivatives::weno5);
			advance(newPhi1, newPhi2, dt, Derivatives::weno5);
			parallelForEach(phi.grid, [&](const VectorDi &coord) { newPhi1[coord] = (phi[coord] * 3 + newPhi2[coord]) / 4; });
			advance(newPhi1, newPhi2, dt, Derivatives::weno5);
			parallelForEach(phi.grid, [&](const VectorDi &coord) { phi[coord] = (phi[coord] + newPhi2[coord] * 2) / 3; });
		}
	}
}

}
