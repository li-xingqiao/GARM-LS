#pragma once

#include "Derivatives.h"

namespace PhysX {

template <int Dim, typename Type> class HermiteIntrpl
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::array<GridBasedData<Dim, Type>, (1 << Dim) - Dim - 1> mixedDrvs;

public:

	HermiteIntrpl(const GridBasedData<Dim, Type> &gbd, const GridBasedData<Dim, Upgrade<Type, Dim>> &drv) : mixedDrvs(makeMixedDrvArray(gbd))
	{
		const double invDx = gbd.grid.invSpacing;
		parallelForEach(gbd.grid, [&](const VectorDi &coord) {
			if constexpr (Dim == 2)
				mixedDrvs[0][coord] = (getC(drv.at(coord + VectorDi::Unit(0)), 1) - getC(drv.at(coord - VectorDi::Unit(0)), 1)) / 2 * invDx;
			else {
				mixedDrvs[0][coord] = (getC(drv.at(coord + VectorDi::Unit(0)), 1) - getC(drv.at(coord - VectorDi::Unit(0)), 1)) / 2 * invDx;
				mixedDrvs[1][coord] = (getC(drv.at(coord + VectorDi::Unit(0)), 2) - getC(drv.at(coord - VectorDi::Unit(0)), 2)) / 2 * invDx;
				mixedDrvs[2][coord] = (getC(drv.at(coord + VectorDi::Unit(1)), 2) - getC(drv.at(coord - VectorDi::Unit(1)), 2)) / 2 * invDx;
				mixedDrvs[3][coord] = (getC(drv.at(coord + VectorDi(1, 1, 0)), 2) - getC(drv.at(coord + VectorDi(-1, 1, 0)), 2) - getC(drv.at(coord + VectorDi(1, -1, 0)), 2) + getC(drv.at(coord + VectorDi(-1, -1, 0)), 2)) / 4 * invDx * invDx;
			}
		});
	}

	Type interpolate(const GridBasedData<Dim, Type> &gbd, const GridBasedData<Dim, Upgrade<Type, Dim>> &drv, const VectorDd &pos) const
	{
		const VectorDi lower = gbd.grid.getLinearLower(pos);
		const VectorDd frac = gbd.grid.getLowerFrac(pos, lower);
		const auto v = getValues(gbd, drv, lower);
		Type value = Zero<Type>();
		for (int i = 0; i < v.size(); i++)
			value += v[i] * weight(frac, i >> Dim, i & (1 << Dim) - 1);
		return value;
	}

	Upgrade<Type, Dim> interpolateGradient(const GridBasedData<Dim, Type> &gbd, const GridBasedData<Dim, Upgrade<Type, Dim>> &drv, const VectorDd &pos) const
	{
		const VectorDi lower = gbd.grid.getLinearLower(pos);
		const VectorDd frac = gbd.grid.getLowerFrac(pos, lower);
		const auto v = getValues(gbd, drv, lower);
		Upgrade<Type, Dim> grad = Zero<Upgrade<Type, Dim>>();
		for (int i = 0; i < v.size(); i++) {
			const int vtx = i >> Dim;
			const int ord = i & (1 << Dim) - 1;
			if constexpr (std::is_scalar_v<Type>)
				grad += v[i] * weightGradient(frac, vtx, ord);
			else
				grad += v[i] * weightGradient(frac, vtx, ord).transpose();
		}
		return grad * gbd.grid.invSpacing;
	}

private:

	auto getValues(const GridBasedData<Dim, Type> &gbd, const GridBasedData<Dim, Upgrade<Type, Dim>> &drv, const VectorDi &lower) const
	{
		const double dx = gbd.grid.spacing;
		std::array<Type, 1 << (Dim << 1)> v;
		for (int i = 0; i < (1 << Dim); i++) {
			if constexpr (Dim == 2) {
				VectorDi coord = lower + VectorDi(i & 1, i >> 1);
				v[i << Dim | 0] = gbd.at(coord);
				v[i << Dim | 1] = getC(drv.at(coord), 0) * dx;
				v[i << Dim | 2] = getC(drv.at(coord), 1) * dx;
				v[i << Dim | 3] = mixedDrvs[0].at(coord) * dx * dx;
			}
			else {
				VectorDi coord = lower + VectorDi(i & 1, i >> 1 & 1, i >> 2);
				v[i << Dim | 0] = gbd.at(coord);
				v[i << Dim | 1] = getC(drv.at(coord), 0) * dx;
				v[i << Dim | 2] = getC(drv.at(coord), 1) * dx;
				v[i << Dim | 3] = mixedDrvs[0].at(coord) * dx * dx;
				v[i << Dim | 4] = getC(drv.at(coord), 2) * dx;
				v[i << Dim | 5] = mixedDrvs[1].at(coord) * dx * dx;
				v[i << Dim | 6] = mixedDrvs[2].at(coord) * dx * dx;
				v[i << Dim | 7] = mixedDrvs[3].at(coord) * dx * dx * dx;
			}
		}
		return v;
	}

	double weight(const VectorDd &frac, const int vtx, const int ord) const
	{
		double wt = 1;
		for (int axis = 0; axis < Dim; axis++) {
			const double x = frac[axis];
			if (ord >> axis & 1)
				wt *= (vtx >> axis & 1) ? x * x * (x - 1) : x * (x - 1) * (x - 1);
			else
				wt *= (vtx >> axis & 1) ? x * x * (3 - x - x) : (x + x + 1) * (x - 1) * (x - 1);
		}
		return wt;
	}

	VectorDd weightGradient(const VectorDd &frac, const int vtx, const int ord) const
	{
		VectorDd wt, wtg;
		for (int axis = 0; axis < Dim; axis++) {
			const double x = frac[axis];
			if (ord >> axis & 1) {
				if (vtx >> axis & 1)
					wt[axis] = x * x * (x - 1), wtg[axis] = x * (3 * x - 2);
				else
					wt[axis] = x * (x - 1) * (x - 1), wtg[axis] = (3 * x - 1) * (x - 1);
			}
			else {
				if (vtx >> axis & 1)
					wt[axis] = x * x * (3 - x - x), wtg[axis] = 6 * x * (1 - x);
				else
					wt[axis] = (x + x + 1) * (x - 1) * (x - 1), wtg[axis] = 6 * x * (x - 1);
			}
		}
		if constexpr (Dim == 2)
			return VectorDd(wtg[0] * wt[1], wt[0] * wtg[1]);
		else
			return VectorDd(wtg[0] * wt[1] * wt[2], wt[0] * wtg[1] * wt[2], wt[0] * wt[1] * wtg[2]);
	}
};

}
