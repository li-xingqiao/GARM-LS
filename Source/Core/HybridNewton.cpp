#include "HybridNewton.h"

#include "Advection.h"
#include "FastMarching.h"

namespace PhysX {

template <int Dim>
void HybridNewton<Dim>::perform(GridBasedData<Dim, double> &phi, GridBasedData<Dim, Vector<Dim, double>> &psi, const int maxSteps, const double cfl)
{
	initializeInterface(phi, psi);
#ifdef _ENABLE_FMM_
	sweep(phi, psi, maxSteps);
#else
	advect(phi, psi, maxSteps, cfl);
#endif
}

template <int Dim>
void HybridNewton<Dim>::perform(GridBasedData<Dim, double> &phi, GridBasedData<Dim, Vector<Dim, double>> &psi, const GridBasedData<Dim, VectorDd> &points, const int maxSteps, const double cfl)
{
	initializeInterface(phi, psi, points);
#ifdef _ENABLE_FMM_
	sweep(phi, psi, maxSteps);
#else
	advect(phi, psi, maxSteps, cfl);
#endif
}

template <int Dim>
void HybridNewton<Dim>::advect(GridBasedData<Dim, double> &phi, GridBasedData<Dim, Vector<Dim, double>> &psi, const int maxSteps, const double cfl) const
{
	const double dt = phi.grid.spacing * cfl;
	for (int iter = 0; iter < maxSteps; iter++) {
		const auto phiView = ScalarFieldView<HermiteIntrpl<Dim, double>>(phi, psi);
		auto newPhi = phi;
		auto newPsi = psi;
		parallelForEach(phi.grid, [&](const VectorDi &coord) {
			if (!_mark[coord]) {
				const int sgn = phi[coord] < 0 ? -1 : 1;
				const VectorDd pos = phi.grid.position(coord);
				const VectorDd targetPos = Advection::trace<3>(pos, Velocity(phiView), -sgn * dt);
				newPhi[coord] = phiView(targetPos) + sgn * dt;
				newPsi[coord] = phiView.gradient(targetPos).normalized();
			}
		});
		phi = newPhi, psi = newPsi;
	}
}

template <int Dim>
void HybridNewton<Dim>::sweep(GridBasedData<Dim, double> &phi, GridBasedData<Dim, Vector<Dim, double>> &psi, const int maxSteps)
{
	if (maxSteps == 0) return;
	FastMarching<Dim> fmm(phi.grid);
	fmm.perform(phi, maxSteps, _mark);
	const auto phiView = ScalarFieldView<CubicIntrpl<Dim>>(phi);
	parallelForEach(phi.grid, [&](const VectorDi &coord) {
		if (!_mark[coord]) {
			// Calculate the gradient.
			const VectorDd pos = phi.grid.position(coord);
			psi[coord] = phiView.gradient(pos).normalized();
		}
	});
}

template <int Dim>
void HybridNewton<Dim>::prepareInterface(const GridBasedData<Dim, double> &phi)
{
	// Find interface cells. (It may be over-smoothed.)
	parallelForEach(phi.grid, [&](const VectorDi &coord) {
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
			if (phi.grid.isValid(nbCoord) && phi[coord] * phi[nbCoord] <= 0) {
				_mark[coord] = true;
				_across[coord] = nbCoord;
				return;
			}
		}
	});
	extrapolate(1);
}


template <int Dim>
void HybridNewton<Dim>::initializeInterface(GridBasedData<Dim, double> &phi, GridBasedData<Dim, Vector<Dim, double>> &psi)
{
	auto phiView = ScalarFieldView<HermiteIntrpl<Dim, double>>(phi, psi);
	// auto phiView = ScalarFieldView<CubicIntrpl<Dim>>(phi);
	initializeInterface(phi, psi, phiView);
}

template <int Dim>
void HybridNewton<Dim>::initializeInterface(GridBasedData<Dim, double> &phi, GridBasedData<Dim, Vector<Dim, double>> &psi, const GridBasedData<Dim, Vector<Dim, double>> &points)
{
	auto phiView = ScalarFieldView<HermiteIntrpl<Dim, double>>(phi, psi);
	// auto phiView = ScalarFieldView<CubicIntrpl<Dim>>(phi);
	// auto phiView = ScalarFieldView<LinearIntrpl<Dim>>(phi);
	initializeInterface(phi, psi, points, phiView);
}

template <int Dim>
void HybridNewton<Dim>::extrapolate(const int maxSteps)
{
	GridBasedData<Dim, uchar> temp(_mark.grid);
	for (int iter = 0; iter < maxSteps; iter++) {
		temp = _mark;
		parallelForEach(_mark.grid, [&](const VectorDi &coord) {
			if (!temp[coord]) {
				_mark[coord] = false;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCoord = Grid<Dim>::neighbor(coord, i);
					if (temp.grid.isValid(nbCoord) && temp[nbCoord]) {
						_mark[coord] = true;
						_across[coord] = _across[nbCoord];
						break;
					}
				}
			}
			else _mark[coord] = true;
		});
	}
}

// template <int Dim>
// template <typename Field>
// Vector<Dim, double> HybridNewton<Dim>::getClosestPosition(const Vector<Dim, int> &coord, const VectorDd &startPos, const VectorDd &initPos, const Field &phiView)


template class HybridNewton<2>;
template class HybridNewton<3>;

}
