#pragma once

#include "GridBasedData.h"

#include <queue>

namespace PhysX {

template <int Dim>
class FastMarching
{
	DECLARE_DIM_TYPES(Dim)

	using HeapElement = std::pair<double, int>;

protected:

	GridBasedData<Dim, double> _tent; // tentative signed distance
	GridBasedData<Dim, uchar> _visited;

	std::vector<int> _intfIndices;
	std::priority_queue<HeapElement, std::vector<HeapElement>, std::greater<HeapElement>> _heap;


public:

	FastMarching(const Grid<Dim> &grid) : _tent(grid), _visited(grid) { }

	void perform(GridBasedData<Dim, double> &phi, int maxSteps);
	void perform(GridBasedData<Dim, double> &phi, int maxSteps, const GridBasedData<Dim, uchar> &mark);

protected:

	void updateNeighbors(const VectorDi &coord);

	double solveEikonalEquation(const VectorDi &coord) const;

	static double solveQuadratic(double p0, const double dx) { return p0 + dx; }
	static double solveQuadratic(double p0, double p1, const double dx);
	static double solveQuadratic(double p0, double p1, double p2, const double dx);
};

}
