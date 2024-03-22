#pragma once

#include "GridBasedData.h"
#include "StaggeredGrid.h"
#include "SGridBasedData.h"

#include <queue>
#include <map>

namespace PhysX {

template <int Dim>
struct RotatedCoord
{
	std::array<Vector<Dim, double>, Dim> dirs;
	void orthogonalize()
	{
		dirs[0] = dirs[0].normalized();
		dirs[1] = (dirs[1] - dirs[1].dot(dirs[0]) * dirs[0]).normalized();
		if constexpr (Dim == 3) {
			dirs[2] = dirs[0].cross(dirs[1]);
		}
	}
};

template <int Dim>
class FMMExtrapolator
{
	DECLARE_DIM_TYPES(Dim)

	using HeapElement = std::pair<double, int>;

protected:

	StaggeredGrid<Dim> _sGrid;
	GridBasedData<Dim, int> _mark;
	GridBasedData<Dim, int> _cellType;
	GridBasedData<Dim, int> _closestNode;
	std::map<int, RotatedCoord<Dim>> _nodeNormals;

	GridBasedData<Dim, double> _virtualPhi;
	std::priority_queue<HeapElement, std::vector<HeapElement>, std::greater<HeapElement>> _heap;

public:

	FMMExtrapolator(const StaggeredGrid<Dim> &sGrid) : _sGrid(sGrid), _mark(sGrid.cellGrid), _cellType(sGrid.cellGrid), _closestNode(sGrid.cellGrid), _virtualPhi(sGrid.cellGrid) { }

	void extrapolate(GridBasedData<Dim, double> &phi, const GridBasedData<Dim, double> &rphi, const GridBasedData<Dim, double> &sphi, const SGridBasedData<Dim, double> &spsi, double theta);

	void setupMark(const GridBasedData<Dim, double> &phi, GridBasedData<Dim, double> &sphi, int maxSteps);
};

}
