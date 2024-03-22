#pragma once

#include <format>
#include "Reinitialization.h"

namespace PhysX {

template <int Dim, int Order>
class LambdaCorr
{
	DECLARE_DIM_TYPES(Dim)

protected:
	GridBasedData<Dim, uchar> _mark;
	// inline double getSign(double val, double sp) { return (std::abs(val) < sp) ? (val / sp + std::sin(std::numbers::pi * val / sp) / std::numbers::pi) : (val > 0. ? 1. : -1.); }
	// inline double getSignDiff(double val, double sp) { return (std::abs(val) < sp) ? (1. + std::cos(std::numbers::pi * val / sp)) / sp : 0.; }
	inline double getSign(double val, double sp, double absg) { return val / std::sqrt(val * val + sp * sp * absg * absg); }
	inline double getSignDiff(double val, double sp, double absg) { return sp * sp * absg * absg / (val * val + sp * sp * absg * absg) * std::sqrt(val * val + sp * sp * absg * absg); }

	void extrapolate(const int maxSteps)
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
							break;
						}
					}
				}
				else _mark[coord] = true;
			});
		}
	}

	void linintg1x1(const GridBasedData<Dim, double> &field, GridBasedData<Dim, double> &result)
	{
		parallelForEach(field.grid, [&](const VectorDi &coord) {
			double sum = 0.0;
			const auto v = field.template stencil<2>(coord - VectorDi::Ones());
			for (double x : v) {
				sum += x;
			}
			if constexpr (Dim == 2) {
				sum += 15. * field[coord];
				result[coord] = sum / 24.;
			}
			else {
				sum += 51. * field[coord];
				result[coord] = sum / 78.; // written down with imagination -_-
			}
			result[coord] *= field.grid.spacing * field.grid.spacing;
		});
	}


public:
	LambdaCorr(const Grid<Dim> &grid) : _mark(grid) { }

	void getCorrection(const GridBasedData<Dim, double> &phi0, const GridBasedData<Dim, double> &absgradphi, const GridBasedData<Dim, double> &phi, GridBasedData<Dim, double> &lambda)
	{
		GridBasedData<Dim, double> orig(lambda.grid), intg(lambda.grid);
		parallelForEach(_mark.grid, [&](const VectorDi &coord) {
			orig[coord] = getSignDiff(phi0[coord], phi.grid.spacing, absgradphi[coord]) * getSignDiff(phi0[coord], phi.grid.spacing, absgradphi[coord]) * absgradphi[coord];
		});
		linintg1x1(orig, intg);
		parallelForEach(_mark.grid, [&](const VectorDi &coord) {
			orig[coord] = getSignDiff(phi0[coord], phi.grid.spacing, absgradphi[coord]) * (phi[coord] - phi0[coord]);
		});
		linintg1x1(orig, lambda);
		parallelForEach(_mark.grid, [&](const VectorDi &coord) {
			if (!_mark[coord]) return;
			lambda[coord] = -2 * lambda[coord] / intg[coord];
		});
	}

	template <int Len>
	void getAbsGrad(const GridBasedData<Dim, double> &field, GridBasedData<Dim, double> &result, std::pair<double, double>(*getDrv)(const std::array<double, Len> &, const double))
	{
		auto square = [](double x) -> double { return x * x; };
		parallelForEach(field.grid, [&](const VectorDi &coord) {
			std::pair<double, double> drv[Dim];
			for (int axis = 0; axis < Dim; axis++)
				drv[axis] = getDrv(field.template stencil<Len - 1>(coord - VectorDi::Unit(axis) * (Len >> 1), axis), field.grid.invSpacing);
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
			result[coord] = std::sqrt(sum);
		});
	}

	void perform(GridBasedData<Dim, double> &phi, const int maxSteps, const double cfl = .5)
	{
		static_assert(Order == 1 || Order == 3 || Order == 5, "Order must be 1, 3 or 5.");
		GridBasedData<Dim, double> phi0(phi), absgrad0(phi.grid), absgrad(phi.grid), lambda(phi.grid);
		if constexpr (Order == 1) {
			getAbsGrad(phi, absgrad0, Derivatives::upwind1);
		}
		else if constexpr (Order == 3) {
			getAbsGrad(phi, absgrad0, Derivatives::eno3);
		}
		else {
			getAbsGrad(phi, absgrad0, Derivatives::weno5);
		}

		// parallelForEach(_mark.grid, [&](const VectorDi &coord) {
		// 	_mark[coord] = true;
		// });
		parallelForEach(_mark.grid, [&](const VectorDi &coord) {
			_mark[coord] = std::abs(phi[coord]) <= _mark.grid.spacing * absgrad[coord] ? true : false;
		});
		// parallelForEach(phi.grid, [&](const VectorDi &coord) {
		// 	for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
		// 		const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
		// 		if (phi.grid.isValid(nbCoord) && phi[coord] * phi[nbCoord] <= 0) {
		// 			_mark[coord] = true;
		// 		}
		// 		else {
		// 			_mark[coord] = false;
		// 		}
		// 	}
		// });
		// extrapolate(1);

		const double dt = phi.grid.spacing * cfl;
		for (int iter = 0; iter < maxSteps; iter++) {
			if constexpr (Order == 1) {
				GridBasedData<Dim, double> newPhi(phi.grid);
				getAbsGrad(phi, newPhi, Derivatives::upwind1);
				parallelForEach(phi.grid, [&](const VectorDi &coord) { phi[coord] = phi[coord] + dt * getSign(phi[coord], phi.grid.spacing, newPhi[coord]) * (1 - newPhi[coord]); });
			}
			else if constexpr (Order == 3) {
				GridBasedData<Dim, double> newPhi1(phi.grid);
				GridBasedData<Dim, double> newPhi2(phi.grid);
				getAbsGrad(phi, newPhi1, Derivatives::eno3);
				parallelForEach(phi.grid, [&](const VectorDi &coord) { newPhi1[coord] = phi[coord] + dt * getSign(phi[coord], phi.grid.spacing, newPhi1[coord]) * (1 - newPhi1[coord]); });
				getAbsGrad(newPhi1, newPhi2, Derivatives::eno3);
				parallelForEach(phi.grid, [&](const VectorDi &coord) { newPhi2[coord] = newPhi1[coord] + dt * getSign(phi[coord], phi.grid.spacing, newPhi2[coord]) * (1 - newPhi2[coord]); });
				parallelForEach(phi.grid, [&](const VectorDi &coord) { phi[coord] = (phi[coord] + newPhi2[coord]) / 2; });
			}
			else {
				GridBasedData<Dim, double> newPhi1(phi.grid);
				GridBasedData<Dim, double> newPhi2(phi.grid);
				getAbsGrad(phi, newPhi1, Derivatives::weno5);
				parallelForEach(phi.grid, [&](const VectorDi &coord) { newPhi1[coord] = phi[coord] + dt * getSign(phi[coord], phi.grid.spacing, newPhi1[coord]) * (1 - newPhi1[coord]); });
				getAbsGrad(newPhi1, newPhi2, Derivatives::weno5);
				parallelForEach(phi.grid, [&](const VectorDi &coord) { newPhi2[coord] = newPhi1[coord] + dt * getSign(phi[coord], phi.grid.spacing, newPhi2[coord]) * (1 - newPhi2[coord]); });
				parallelForEach(phi.grid, [&](const VectorDi &coord) { newPhi1[coord] = (phi[coord] * 3 + newPhi2[coord]) / 4; });
				getAbsGrad(newPhi1, newPhi2, Derivatives::weno5);
				parallelForEach(phi.grid, [&](const VectorDi &coord) { newPhi2[coord] = newPhi1[coord] + dt * getSign(phi[coord], phi.grid.spacing, newPhi2[coord]) * (1 - newPhi2[coord]); });
				parallelForEach(phi.grid, [&](const VectorDi &coord) { phi[coord] = (phi[coord] + newPhi2[coord] * 2) / 3; });
			}

			if constexpr (Order == 1) {
				getAbsGrad(phi, absgrad, Derivatives::upwind1);
			}
			else if constexpr (Order == 3) {
				getAbsGrad(phi, absgrad, Derivatives::eno3);
			}
			else {
				getAbsGrad(phi, absgrad, Derivatives::weno5);
			}
			getCorrection(phi0, absgrad, phi, lambda);

			// GridBasedData<Dim, double> orig(phi.grid), intg0(phi.grid), intg(phi.grid);
			// parallelForEach(_mark.grid, [&](const VectorDi &coord) {
			// 	orig[coord] = (getSign(phi0[coord], phi.grid.spacing, absgrad0[coord]) + 1) / 2;
			// });
			// linintg1x1(orig, intg0);
			// parallelForEach(_mark.grid, [&](const VectorDi &coord) {
			// 	orig[coord] = (getSign(phi[coord], phi.grid.spacing, absgrad[coord]) + 1) / 2;
			// });
			// linintg1x1(orig, intg);
			// double dva = 0.0, sum = 0.0;
			// forEach(_mark.grid, [&](const VectorDi &coord) {
			// 	dva += std::abs(intg[coord] - intg0[coord]);
			// 	sum += std::abs(intg[coord]);
			// });
			// dva /= sum;
			// std::cout << std::format("Deviation: {:.2f}\t", dva);

			parallelForEach(phi.grid, [&](const VectorDi &coord) {
				if (!_mark[coord]) return;
				phi[coord] += lambda[coord] * absgrad[coord] * getSignDiff(phi[coord], phi.grid.spacing, absgrad[coord]);
			});

			// if constexpr (Order == 1) {
			// 	getAbsGrad(phi, absgrad, Derivatives::upwind1);
			// }
			// else if constexpr (Order == 3) {
			// 	getAbsGrad(phi, absgrad, Derivatives::eno3);
			// }
			// else {
			// 	getAbsGrad(phi, absgrad, Derivatives::weno5);
			// }
			// parallelForEach(_mark.grid, [&](const VectorDi &coord) {
			// 	orig[coord] = (getSign(phi[coord], phi.grid.spacing, absgrad[coord]) + 1) / 2;
			// });
			// linintg1x1(orig, intg);
			// dva = 0.0, sum = 0.0;
			// forEach(_mark.grid, [&](const VectorDi &coord) {
			// 	dva += std::abs(intg[coord] - intg0[coord]);
			// 	sum += std::abs(intg[coord]);
			// });
			// dva /= sum;
			// std::cout << std::format("after correction: {:.2f}\n", dva);
		}
	}
};

} // namespace PhysX
