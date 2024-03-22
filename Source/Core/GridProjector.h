#pragma once

#include "SGridBasedData.h"
#include "StaggeredGrid.h"

namespace PhysX {

template <int Dim>
class GridProjector
{
	DECLARE_DIM_TYPES(Dim)

protected:

	GridBasedData<Dim, double> _reducedPressure;
	GridBasedData<Dim, double> _velocityDiv;

	SparseMatrixd _matLaplacian;

	GridBasedData<Dim, uchar> _mark;

public:

	GridProjector(const StaggeredGrid<Dim> &sGrid);

	void project(
		SGridBasedData<Dim, double> &velocity,
		const GridBasedData<Dim, double> &levelSet,
		const SGridBasedData<Dim, double> &boundaryFraction,
		const SGridBasedData<Dim, double> &boundaryVelocity,
		const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump = nullptr);

	void reproject(
		SGridBasedData<Dim, double> &velocity,
		const GridBasedData<Dim, double> &levelSet,
		const SGridBasedData<Dim, double> &boundaryFraction,
		const SGridBasedData<Dim, double> &boundaryVelocity);

	void buildLinearSystem(
		SGridBasedData<Dim, double> &velocity,
		const GridBasedData<Dim, double> &levelSet,
		const SGridBasedData<Dim, double> &boundaryFraction,
		const SGridBasedData<Dim, double> &boundaryVelocity,
		const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump = nullptr);

	void setRightHandSide(
		SGridBasedData<Dim, double> &velocity,
		const SGridBasedData<Dim, double> &boundaryFraction,
		const SGridBasedData<Dim, double> &boundaryVelocity);

	void applyPressureGradient(
		SGridBasedData<Dim, double> &velocity,
		const GridBasedData<Dim, double> &levelSet,
		const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump = nullptr);

	void setUnknowns(const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction);
	void solveLinearSystem();
};

}
