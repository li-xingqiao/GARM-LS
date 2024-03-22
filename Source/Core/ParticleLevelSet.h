#pragma once

#include "Advection.h"
#include "GridBasedData.h"
#include "Interpolator.h"
#include "Reinitialization.h"
#include "ScalarField.h"
#include "VectorField.h"

#include <algorithm>

namespace PhysX::PLS {

template <int Dim>
void advectParticles(std::vector<Vector<Dim, double>> &pPos, const VectorField<Dim> &flow, const double dt) {
	DECLARE_DIM_TYPES(Dim)
	for (auto &pos : pPos) {
		const VectorDd vel0 = flow(pos);
		const VectorDd vel1 = flow(pos + vel0 * dt * .5);
		const VectorDd vel2 = flow(pos - vel0 * dt + 2 * vel1 * dt);
		pos += (vel0 + 4 * vel1 + vel2) * dt / 6;
	}
}

template <int Dim>
void correctLevelSet(GridBasedData<Dim, double> &phi, const std::vector<Vector<Dim, double>> &pPos, const std::vector<double> &pRad, const std::vector<int> &pSign) {
	DECLARE_DIM_TYPES(Dim)
	GridBasedData<Dim, double> phiPositive = phi;
	GridBasedData<Dim, double> phiNegative = phi;
	for (int i = 0; i < pPos.size(); i++) {
		const VectorDd pos = pPos[i];
		const double rad = pRad[i];
		const int sign = pSign[i];
		const double tPhi = ScalarFieldView<LinearIntrpl<Dim>>(phi)(pos);
		if (tPhi * sign < 0 && std::abs(tPhi) > rad) {
			for (const auto &coord : LinearIntrpl<Dim>::points(phi.grid, pos)) {
				if (sign > 0) {
					phiPositive[coord] = std::max(phiPositive[coord], rad - (pos - phi.grid.position(coord)).norm());
				} else {
					phiNegative[coord] = std::min(phiNegative[coord], (pos - phi.grid.position(coord)).norm() - rad);
				}
			}
		}
	}
	parallelForEach(phi.grid, [&](const VectorDi &coord) {
		if (std::abs(phiPositive[coord]) <= std::abs(phiNegative[coord])) {
			phi[coord] = phiPositive[coord];
		} else {
			phi[coord] = phiNegative[coord];
		}
	});
}

template <int Dim>
void resample(const GridBasedData<Dim, double> &phi, std::vector<Vector<Dim, double>> &pPos, std::vector<double> &pRad, std::vector<int> &pSign, const int nPPC) {
	DECLARE_DIM_TYPES(Dim)
	pPos.clear();
	pRad.clear();
	pSign.clear();
	const double dx = phi.grid.spacing;
	forEach(phi.grid, [&](const VectorDi &cell) {
		if (std::abs(phi[cell]) <= 3 * dx) {
			for (int i = 0; i < nPPC; i++) {
				const VectorDd pos = phi.grid.position(cell) + VectorDd::Random() * .5 * dx;
				const double phiVal = ScalarFieldView<LinearIntrpl<Dim>>(phi)(pos);
				const int sign = phiVal < 0 ? -1 : 1;
				const double rad = std::clamp(sign * phiVal, .1 * dx, .5 * dx);
				pPos.push_back(pos);
				pRad.push_back(rad);
				pSign.push_back(sign);
			}
		}
	});
}

template <int Dim>
void reradius(const GridBasedData<Dim, double> &phi, std::vector<Vector<Dim, double>> &pPos, std::vector<double> &pRad, std::vector<int> &pSign) {
	DECLARE_DIM_TYPES(Dim)
	const double dx = phi.grid.spacing;
	for (int i = 0; i < pPos.size(); i++) {
		const VectorDd pos = pPos[i];
		double &rad = pRad[i];
		const int sign = pSign[i];
		rad = std::clamp(sign * ScalarFieldView<LinearIntrpl<Dim>>(phi)(pos), .1 * dx, .5 * dx);
	}
}

}
