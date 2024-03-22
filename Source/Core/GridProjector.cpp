#include "GridProjector.h"

#include "LinearSystem.h"

namespace PhysX {

template <int Dim>
GridProjector<Dim>::GridProjector(const StaggeredGrid<Dim> &sGrid) :
	_reducedPressure(sGrid.cellGrid),
	_velocityDiv(sGrid.cellGrid),
	_matLaplacian(sGrid.cellCount(), sGrid.cellCount()),
	_mark(sGrid.cellGrid)
{
	std::cout << "** Using internal ICPCG solver..." << std::endl;
}

template <int Dim>
void GridProjector<Dim>::project(SGridBasedData<Dim, double> &velocity, const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction, const SGridBasedData<Dim, double> &boundaryVelocity, const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump)
{
	setUnknowns(levelSet, boundaryFraction);
	buildLinearSystem(velocity, levelSet, boundaryFraction, boundaryVelocity, pressureJump);
	solveLinearSystem();
	applyPressureGradient(velocity, levelSet, pressureJump);
}

template <int Dim>
void GridProjector<Dim>::reproject(SGridBasedData<Dim, double> &velocity, const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction, const SGridBasedData<Dim, double> &boundaryVelocity)
{
	setRightHandSide(velocity, boundaryFraction, boundaryVelocity);
	solveLinearSystem();
	applyPressureGradient(velocity, levelSet);
}

template <int Dim>
void GridProjector<Dim>::buildLinearSystem(SGridBasedData<Dim, double> &velocity, const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction, const SGridBasedData<Dim, double> &boundaryVelocity, const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump)
{
	std::vector<Tripletd> elements;

	forEach(levelSet.grid, [&](const VectorDi &cell) {
		const int idx = int(levelSet.grid.index(cell));
		double diagCoeff = 0;
		double div = 0;
		if (_mark[cell]) {
			for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
				const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
				const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
				const int side = StaggeredGrid<Dim>::cellFaceSide(i);
				const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
				const double weight = 1 - boundaryFraction[axis][face];
				if (weight > 0) {
					if (_mark[nbCell]) {
						const int nbIdx = int(levelSet.grid.index(nbCell));
						diagCoeff += weight;
						elements.push_back(Tripletd(idx, nbIdx, -weight));
					}
					else {
						double theta = levelSet[cell] / (levelSet[cell] - levelSet[nbCell]);
						if (pressureJump)
							div -= pressureJump(cell, nbCell, theta);
						const double intfCoef = 1 / std::max(theta, .001);
						diagCoeff += weight * intfCoef;
					}
					div += side * weight * velocity[axis][face];
				}
				if (weight < 1)
					div += side * (1 - weight) * boundaryVelocity[axis][face];
			}
		}
		else
			diagCoeff = 1;

		_velocityDiv[cell] = div;
		elements.push_back(Tripletd(idx, idx, diagCoeff));
	});

	_matLaplacian.setFromTriplets(elements.begin(), elements.end());
}

template <int Dim>
void GridProjector<Dim>::setRightHandSide(SGridBasedData<Dim, double> &velocity, const SGridBasedData<Dim, double> &boundaryFraction, const SGridBasedData<Dim, double> &boundaryVelocity)
{
	forEach(_mark.grid, [&](const VectorDi &cell) {
		if (!_mark[cell]) return;
		double div = 0;
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
			const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
			const int side = StaggeredGrid<Dim>::cellFaceSide(i);
			const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
			const double weight = 1 - boundaryFraction[axis][face];
			if (weight > 0)
				div += side * weight * velocity[axis][face];
			if (weight < 1)
				div += side * (1 - weight) * boundaryVelocity[axis][face];
		}
		_velocityDiv[cell] = div;
	});
}

template <int Dim>
void GridProjector<Dim>::applyPressureGradient(SGridBasedData<Dim, double> &velocity, const GridBasedData<Dim, double> &levelSet, const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump)
{
	parallelForEach(velocity.grids, [&](const int axis, const VectorDi &face) {
		const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
		const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
		if (!levelSet.grid.isValid(cell0) || !levelSet.grid.isValid(cell1)) return;

		if (!_mark[cell0] && !_mark[cell1]) return;

		const double p0 = _mark[cell0] ? _reducedPressure[cell0] : 0;
		const double p1 = _mark[cell1] ? _reducedPressure[cell1] : 0;

		if (_mark[cell0] && _mark[cell1])
			velocity[axis][face] += p1 - p0;
		else {
			const double phi0 = levelSet[cell0];
			const double phi1 = levelSet[cell1];

			if (phi0 <= 0 && phi1 <= 0) return;

			double theta = phi0 / (phi0 - phi1);

			if (pressureJump)
				velocity[axis][face] += (phi0 <= 0 ? -1 : 1) * pressureJump(cell0, cell1, theta);

			const double intfCoef = 1 / std::max((phi0 <= 0 ? theta : 1 - theta), .001);
			velocity[axis][face] += (p1 - p0) * intfCoef;
		}
	});
}

template<int Dim>
void GridProjector<Dim>::setUnknowns(const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction)
{
	forEach(levelSet.grid, [&](const VectorDi &cell) {
		if (levelSet[cell] <= 0) {
			bool flag = false;
			for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
				const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
				const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
				if (boundaryFraction[axis][face] < 1) {
					flag = true;
					break;
				}
			}
			_mark[cell] = flag;
		}
		else
			_mark[cell] = false;
	});
}

template <int Dim>
void GridProjector<Dim>::solveLinearSystem()
{
	LinearSystem::solve(_matLaplacian, _reducedPressure.asVectorXd(), _velocityDiv.asVectorXd(), 2000, 1e-5);
	std::cout << " (Pressure)";
}

template class GridProjector<2>;
template class GridProjector<3>;
}
