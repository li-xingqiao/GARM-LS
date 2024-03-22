#pragma once

#include "Derivatives.h"
#include "LagrangeIntrpl.h"

namespace PhysX {

template <int Dim, typename Type>
class EnhancedLinearIntrpl
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::array<GridBasedData<Dim, Type>, Dim> drvs;

public:

	EnhancedLinearIntrpl(const GridBasedData<Dim, Type> &gbd) : drvs(makeDrvArray(gbd))
	{
		GridBasedData<Dim, Type> drv(gbd.grid);
		for (int axis = 0; axis < Dim; axis++) {
			Derivatives::computeSecond<2>(gbd, drv, axis);
			parallelForEach(gbd.grid, [&](const VectorDi &coord) {
				const auto v = drv.template stencil<1>(coord);
				int minIdx = 0;
				double minMod = norm(v[0]);
				for (int i = 1; i < v.size(); i++)
					if (double t = norm(v[i]); t < minMod)
						minIdx = i, minMod = t;
				drvs[axis][coord] = v[minIdx];
			});
		}
	}

	Type interpolate(const GridBasedData<Dim, Type> &gbd, const VectorDd &pos) const
	{
		const VectorDi lower = gbd.grid.getLinearLower(pos);
		const double dx = gbd.grid.spacing;
		Type val = LinearIntrpl<Dim>::interpolate(gbd, pos);
		const VectorDd delta = pos - gbd.grid.position(lower);
		for (int axis = 0; axis < Dim; axis++)
			val -= drvs[axis].at(lower) * delta[axis] * (dx - delta[axis]) / 2;
		return val;
	}
};

}
