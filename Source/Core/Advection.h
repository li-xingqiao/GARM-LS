#pragma once

#include "VectorField.h"

namespace PhysX::Advection {

template <typename Scheme, int Dim, typename Type>
inline void solve(GridBasedData<Dim, Type> &gbd, const VectorField<Dim> &flow, const double dt)
{
	GridBasedData<Dim, Type> newGbd(gbd.grid);
	Scheme::advect(gbd, newGbd, flow, dt);
	gbd = newGbd;
}

template <typename Scheme, int Dim, typename Type>
inline void solve(SGridBasedData<Dim, Type> &sgbd, const VectorField<Dim> &flow, const double dt)
{
	SGridBasedData<Dim, Type> newSgbd(sgbd.grids);
	for (int axis = 0; axis < Dim; axis++)
		Scheme::advect(sgbd[axis], newSgbd[axis], flow, dt);
	sgbd = newSgbd;
}

template <typename Scheme, int Dim, typename Type>
inline void solve(GridBasedData<Dim, Type> &gbd, GridBasedData<Dim, Upgrade<Type, Dim>> &drv, const VectorField<Dim> &flow, const double dt)
{
	GridBasedData<Dim, Type> newGbd(gbd.grid);
	GridBasedData<Dim, Upgrade<Type, Dim>> newDrv(drv.grid);
	Scheme::advect(gbd, drv, newGbd, newDrv, flow, dt);
	gbd = newGbd, drv = newDrv;
}

template <int RkOrder, int Dim>
inline Vector<Dim, double> trace(const Vector<Dim, double> &startPos, const VectorField<Dim> &flow, const double dt)
{
	DECLARE_DIM_TYPES(Dim)

	if constexpr (RkOrder == 1) return startPos + flow(startPos) * dt;
	else if constexpr (RkOrder == 2) { // the TVD Runge-Kutta scheme of second order
		const VectorDd vel0 = flow(startPos);
		const VectorDd vel1 = flow(startPos + vel0 * dt);
		return startPos + (vel0 + vel1) * dt / 2;
	}
	else if constexpr (RkOrder == 3) { // the TVD Runge-Kutta scheme of third order
		const VectorDd vel0 = flow(startPos);
		const VectorDd vel1 = flow(startPos + vel0 * dt);
		const VectorDd vel2 = flow(startPos + (vel0 + vel1) * dt / 4);
		return startPos + (vel0 + vel1 + 4 * vel2) * dt / 6;
	}
	else if constexpr (RkOrder == 4) { // the classical fourth order Runge-Kutta scheme
		const VectorDd vel0 = flow(startPos);
		const VectorDd vel1 = flow(startPos + vel0 * dt / 2);
		const VectorDd vel2 = flow(startPos + vel1 * dt / 2);
		const VectorDd vel3 = flow(startPos + vel2 * dt);
		return startPos + (vel0 + 2 * vel1 + 2 * vel2 + vel3) * dt / 6;
	}
}

template <int RkOrder, int Dim>
inline std::pair<Vector<Dim, double>, Matrix<Dim, double>> traceDrv(const Vector<Dim, double> &startPos, const VectorField<Dim> &flow, const double dt)
{
	DECLARE_DIM_TYPES(Dim)

	if constexpr (RkOrder == 1) {
		const VectorDd vel0 = flow(startPos);
		const MatrixDd grd0 = flow.gradient(startPos);
		const VectorDd pos1 = startPos + vel0 * dt;
		const MatrixDd rot1 = MatrixDd::Identity() + grd0 * dt;
		return std::pair(pos1, rot1);
	}
	else if constexpr (RkOrder == 2) { // the TVD Runge-Kutta scheme of second order
		const VectorDd vel0 = flow(startPos);
		const MatrixDd grd0 = flow.gradient(startPos);
		const VectorDd pos1 = startPos + vel0 * dt;
		const MatrixDd rot1 = MatrixDd::Identity() + grd0 * dt;

		const VectorDd vel1 = flow(pos1);
		const MatrixDd grd1 = flow.gradient(pos1);
		const VectorDd pos2 = startPos + (vel0 + vel1) * dt / 2;
		const MatrixDd rot2 = MatrixDd::Identity() + (grd0 + grd1 * rot1) * dt / 2;
		return std::pair(pos2, rot2);
	}
	else if constexpr (RkOrder == 3) { // the TVD Runge-Kutta scheme of third order
		const VectorDd vel0 = flow(startPos);
		const MatrixDd grd0 = flow.gradient(startPos);
		const VectorDd pos1 = startPos + vel0 * dt;
		const MatrixDd rot1 = MatrixDd::Identity() + grd0 * dt;

		const VectorDd vel1 = flow(pos1);
		const MatrixDd grd1 = flow.gradient(pos1);
		const VectorDd pos2 = startPos + (vel0 + vel1) * dt / 4;
		const MatrixDd rot2 = MatrixDd::Identity() + (grd0 + grd1 * rot1) * dt / 4;

		const VectorDd vel2 = flow(pos2);
		const MatrixDd grd2 = flow.gradient(pos2);
		const VectorDd pos3 = startPos + (vel0 + vel1 + 4 * vel2) * dt / 6;
		const MatrixDd rot3 = MatrixDd::Identity() + (grd0 + grd1 * rot1 + 4 * grd2 * rot2) * dt / 6;
		return std::pair(pos3, rot3);
	}
	else { // the classical fourth order Runge-Kutta scheme
		const VectorDd vel0 = flow(startPos);
		const MatrixDd grd0 = flow.gradient(startPos);
		const VectorDd pos1 = startPos + vel0 * dt / 2;
		const MatrixDd rot1 = MatrixDd::Identity() + grd0 * dt / 2;

		const VectorDd vel1 = flow(pos1);
		const MatrixDd grd1 = flow.gradient(pos1);
		const VectorDd pos2 = startPos + vel1 * dt / 2;
		const MatrixDd rot2 = MatrixDd::Identity() + grd1 * rot1 * dt / 2;

		const VectorDd vel2 = flow(pos2);
		const MatrixDd grd2 = flow.gradient(pos2);
		const VectorDd pos3 = startPos + vel2 * dt;
		const MatrixDd rot3 = MatrixDd::Identity() + grd2 * rot2 * dt;

		const VectorDd vel3 = flow(pos3);
		const MatrixDd grd3 = flow.gradient(pos3);
		const VectorDd pos4 = startPos + (vel0 + 2 * vel1 + 2 * vel2 + vel3) * dt / 6;
		const MatrixDd rot4 = MatrixDd::Identity() + (grd0 + 2 * grd1 * rot1 + 2 * grd2 * rot2 + grd3 * rot3) * dt / 6;
		return std::pair(pos4, rot4);
	}
}

}

namespace PhysX {

template <int Dim, int RkOrder, typename Intrpl>
class SemiLagrangian
{
	DECLARE_DIM_TYPES(Dim)

	static_assert(1 <= RkOrder && RkOrder <= 4, "Runge-Kutta order must be 1, 2, 3 or 4.");

public:

	template <typename Type>
	static void advect(const GridBasedData<Dim, Type> &gbd, GridBasedData<Dim, Type> &newGbd, const VectorField<Dim> &flow, const double dt)
	{
		const auto intrpl = makeIntrpl<Intrpl>(gbd);
		parallelForEach(gbd.grid, [&](const VectorDi &coord) {
			const VectorDd pos = gbd.grid.position(coord);
			newGbd[coord] = intrpl->interpolate(gbd, Advection::trace<RkOrder>(pos, flow, -dt));
		});
	}

	template <typename Type>
	static void advect(const GridBasedData<Dim, Type> &gbd, const GridBasedData<Dim, Upgrade<Type, Dim>> &drv, GridBasedData<Dim, Type> &newGbd, GridBasedData<Dim, Upgrade<Type, Dim>> &newDrv, const VectorField<Dim> &flow, const double dt)
	{
		const auto intrpl = makeIntrpl<Intrpl>(gbd, drv);
		parallelForEach(gbd.grid, [&](const VectorDi &coord) {
			const VectorDd pos = gbd.grid.position(coord);
			const auto [targetPos, targetGrd] = Advection::traceDrv<RkOrder>(pos, flow, -dt);
			newGbd[coord] = intrpl->interpolate(gbd, drv, targetPos);
			if constexpr (std::is_scalar_v<Type>)
				newDrv[coord] = targetGrd.transpose() * intrpl->interpolateGradient(gbd, drv, targetPos);
			else
				newDrv[coord] = intrpl->interpolateGradient(gbd, drv, targetPos) * targetGrd;
		});
	}
};

template <int Dim, int RkOrder>
class MacCormack
{
	DECLARE_DIM_TYPES(Dim)

	static_assert(1 <= RkOrder && RkOrder <= 4, "Runge-Kutta order must be 1, 2, 3 or 4.");

public:

	template <typename Type>
	static void advect(const GridBasedData<Dim, Type> &gbd, GridBasedData<Dim, Type> &newGbd, const VectorField<Dim> &flow, const double dt)
	{
		GridBasedData<Dim, Type> forwardGbd(gbd.grid);
		SemiLagrangian<Dim, RkOrder, LinearIntrpl<Dim>>::advect(gbd, forwardGbd, flow, dt);
		SemiLagrangian<Dim, RkOrder, LinearIntrpl<Dim>>::advect(forwardGbd, newGbd, flow, -dt);
		parallelForEach(gbd.grid, [&](const VectorDi &coord) {
			auto &val = newGbd[coord];
			val = forwardGbd[coord] + (gbd[coord] - val) * .5;
		});
	}
};

}
