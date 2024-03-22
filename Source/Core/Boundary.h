#pragma once

#include "ImplicitSurface.h"
#include "SGridBasedData.h"
#include "Surface.h"
#include "StaggeredGrid.h"

namespace PhysX {

template <int Dim>
class Boundary
{
	DECLARE_DIM_TYPES(Dim)

public:

	const ImplicitBox<Dim> domainBox;
	GridBasedData<Dim, double> nodeDist;
	GridBasedData<Dim, double> cellDist;

	SGridBasedData<Dim, double> fraction;
	SGridBasedData<Dim, double> velocity;
	SGridBasedData<Dim, double> normal;

public:

	Boundary(const StaggeredGrid<Dim> &sGrid);

	void reset();
	void unions(const Surface<Dim> &surface);
	void intersects(const Surface<Dim> &surface);
	void finish(const StaggeredGrid<Dim> &sGrid);

	void enforce(SGridBasedData<Dim, double> &fluidVelocity) const;
	void enforce(GridBasedData<Dim, Vector<Dim, double>> &fluidRefMap) const;

	double getFaceFraction(const int axis, const VectorDi &face) const;
};

}
