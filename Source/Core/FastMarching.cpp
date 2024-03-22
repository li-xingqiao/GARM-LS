#include "FastMarching.h"

#include <algorithm>
#include <limits>

#include <cmath>

namespace PhysX {

template <int Dim>
void FastMarching<Dim>::perform(GridBasedData<Dim, double> &phi, int maxSteps)
{
	const double bandWidth = maxSteps * phi.grid.spacing;
	_tent.setConstant(bandWidth > 0 ? bandWidth : std::numeric_limits<double>::infinity());
	_visited.setZero();
	// Initialize interface cells.
	forEach(phi.grid, [&](const VectorDi &coord) {
		VectorDd tempPhi = VectorDd::Ones() * std::numeric_limits<double>::infinity();
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
			if (phi.grid.isValid(nbCoord) && phi[coord] * phi[nbCoord] <= 0) {
				const int axis = Grid<Dim>::neighborAxis(i);
				tempPhi[axis] = std::min(tempPhi[axis], phi[coord] / (phi[coord] - phi[nbCoord]) * phi.grid.spacing);
			}
		}
		if (tempPhi.array().isFinite().any()) {
			_tent[coord] = 1 / tempPhi.cwiseInverse().norm();
			_visited[coord] = true;
			_intfIndices.push_back(int(phi.index(coord)));
		}
	});
	// Perform the algorithm.
	for (const auto index : _intfIndices)
		updateNeighbors(_tent.coordinate(index));
	while (!_heap.empty()) {
		const double val = _heap.top().first;
		const VectorDi coord = _tent.coordinate(_heap.top().second);
		_heap.pop();
		if (_tent[coord] != val) continue;
		_visited[coord] = true;
		updateNeighbors(coord);
	}
	parallelForEach(phi.grid, [&](const VectorDi &coord) {
		phi[coord] = (phi[coord] <= 0 ? -1 : 1) * _tent[coord];
	});
}

template <int Dim>
void FastMarching<Dim>::perform(GridBasedData<Dim, double> &phi, int maxSteps, const GridBasedData<Dim, uchar> &mark)
{
	const double bandWidth = maxSteps * phi.grid.spacing;
	_tent.setConstant(bandWidth >= 0 ? bandWidth : std::numeric_limits<double>::infinity());
	_visited.setZero();
	// Initialize interface cells.
	forEach(phi.grid, [&](const VectorDi &coord) {
		if (mark[coord]) {
			_tent[coord] = std::abs(phi[coord]);
			_visited[coord] = true;
			_intfIndices.push_back(int(phi.index(coord)));
		}
	});
	// Perform the algorithm.
	for (const auto index : _intfIndices)
		updateNeighbors(_tent.coordinate(index));
	while (!_heap.empty()) {
		const double val = _heap.top().first;
		const VectorDi coord = _tent.coordinate(_heap.top().second);
		_heap.pop();
		if (_tent[coord] != val) continue;
		_visited[coord] = true;
		updateNeighbors(coord);
	}
	parallelForEach(phi.grid, [&](const VectorDi &coord) {
		phi[coord] = (phi[coord] <= 0 ? -1 : 1) * _tent[coord];
	});
}

template <int Dim>
void FastMarching<Dim>::updateNeighbors(const VectorDi &coord)
{
	for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
		const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
		if (!_tent.grid.isValid(nbCoord) || _visited[nbCoord]) continue;
		if (const double temp = solveEikonalEquation(nbCoord); temp < _tent[nbCoord]) {
			_tent[nbCoord] = temp;
			_heap.push(HeapElement(temp, int(_tent.index(nbCoord))));
		}
	}
}

template <int Dim>
double FastMarching<Dim>::solveEikonalEquation(const VectorDi &coord) const
{
	VectorDd tempPhi = VectorDd::Ones() * std::numeric_limits<double>::infinity();
	VectorDi mark = VectorDi::Zero();
	for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
		const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
		if (_tent.grid.isValid(nbCoord) && _visited[nbCoord]) {
			const int axis = Grid<Dim>::neighborAxis(i);
			tempPhi[axis] = std::min(tempPhi[axis], _tent[nbCoord]);
		}
	}
	double newPhi;
	if constexpr (Dim == 2) newPhi = solveQuadratic(tempPhi.x(), tempPhi.y(), _tent.grid.spacing);
	else newPhi = solveQuadratic(tempPhi.x(), tempPhi.y(), tempPhi.z(), _tent.grid.spacing);
	if (!std::isfinite(newPhi)) {
		std::cerr << "Error: [FastMarchingReinitializer] failed to solve Eikonal equation." << std::endl;
		std::exit(-1);
	}
	return newPhi;
}

template <int Dim>
double FastMarching<Dim>::solveQuadratic(double p0, double p1, const double dx)
{
	if (p0 > p1) std::swap(p0, p1);
	if (std::isinf(p1) || p1 - p0 > dx) return solveQuadratic(p0, dx);
	else return ((p0 + p1) + std::sqrt(2 * dx * dx - (p0 - p1) * (p0 - p1))) * .5;
}

template <int Dim>
double FastMarching<Dim>::solveQuadratic(double p0, double p1, double p2, const double dx)
{
	if (p0 > p1) std::swap(p0, p1);
	if (p1 > p2) std::swap(p1, p2);
	if (std::isinf(p2) || (p2 - p0) * (p2 - p0) + (p2 - p1) * (p2 - p1) > dx * dx) return solveQuadratic(p0, p1, dx);
	else return (p0 + p1 + p2 + std::sqrt((p0 + p1 + p2) * (p0 + p1 + p2) - 3 * (p0 * p0 + p1 * p1 + p2 * p2 - dx * dx))) / 3;
}

template class FastMarching<2>;
template class FastMarching<3>;

}
