#include "Projector.h"

#include "LinearSystem.h"

namespace PhysX {

template <int Dim>
Projector<Dim>::Projector(const StaggeredGrid<Dim> &sGrid) : _grid2mat(sGrid.cellGrid)
{
	std::cout << "** Using internal ICPCG solver..." << std::endl;
}

template <int Dim>
void Projector<Dim>::project(SGridBasedData<Dim, double> &velocity, const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction, const SGridBasedData<Dim, double> &boundaryVelocity, const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump)
{
	setUnknowns(levelSet, boundaryFraction);
	buildLinearSystem(velocity, levelSet, boundaryFraction, boundaryVelocity, pressureJump);
	solveLinearSystem();
	applyPressureGradient(velocity, levelSet, pressureJump);
}

template <int Dim>
void Projector<Dim>::reproject(SGridBasedData<Dim, double> &velocity, const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction, const SGridBasedData<Dim, double> &boundaryVelocity)
{
	setRightHandSide(velocity, boundaryFraction, boundaryVelocity);
	solveLinearSystem();
	applyPressureGradient(velocity, levelSet);
}

template <int Dim>
void Projector<Dim>::buildLinearSystem(SGridBasedData<Dim, double> &velocity, const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction, const SGridBasedData<Dim, double> &boundaryVelocity, const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump)
{
	std::vector<Tripletd> elements;

	for (int r = 0; r < _rdP.size(); r++) {
		const VectorDi cell = levelSet.grid.coordinate(_mat2grid[r]);
		double diagCoeff = 0;
		double div = 0;
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
			const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
			const int side = StaggeredGrid<Dim>::cellFaceSide(i);
			const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
			const double weight = 1 - boundaryFraction[axis][face];
			if (weight > 0) {
				const int c = _grid2mat[nbCell];
				if (c >= 0) {
					diagCoeff += weight;
					elements.push_back(Tripletd(r, c, -weight));
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
		_rhs[r] = div;
		elements.push_back(Tripletd(r, r, diagCoeff));
	}

	_matLaplacian.setFromTriplets(elements.begin(), elements.end());
}

template <int Dim>
void Projector<Dim>::setRightHandSide(SGridBasedData<Dim, double> &velocity, const SGridBasedData<Dim, double> &boundaryFraction, const SGridBasedData<Dim, double> &boundaryVelocity)
{
	for (int r = 0; r < _rhs.size(); r++) {
		const VectorDi cell = _grid2mat.grid.coordinate(_mat2grid[r]);
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
		_rhs[r] = div;
	}
}

template <int Dim>
void Projector<Dim>::applyPressureGradient(SGridBasedData<Dim, double> &velocity, const GridBasedData<Dim, double> &levelSet, const std::function<double(const VectorDi &, const VectorDi &, double &)> pressureJump)
{
	parallelForEach(velocity.grids, [&](const int axis, const VectorDi &face) {
		const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
		const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
		if (!levelSet.grid.isValid(cell0) || !levelSet.grid.isValid(cell1)) return;

		const int id0 = _grid2mat[cell0];
		const int id1 = _grid2mat[cell1];
		if (id0 < 0 && id1 < 0) return;

		const double p0 = id0 >= 0 ? _rdP[_grid2mat[cell0]] : 0;
		const double p1 = id1 >= 0 ? _rdP[_grid2mat[cell1]] : 0;

		if (id0 >= 0 && id1 >= 0)
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
void Projector<Dim>::setUnknowns(const GridBasedData<Dim, double> &levelSet, const SGridBasedData<Dim, double> &boundaryFraction)
{
	_grid2mat.setConstant(-1);
	_mat2grid.clear();

	int n = 0;

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
			if (flag) {
				_grid2mat[cell] = n++;
				_mat2grid.push_back(int(levelSet.grid.index(cell)));
			}
		}
	});

	_matLaplacian.resize(n, n);
	_rdP.resize(n);
	_rhs.resize(n);

	_rdP.setZero();
}

template <int Dim>
void Projector<Dim>::solveLinearSystem()
{
	LinearSystem::solve(_matLaplacian, _rdP, _rhs, 2000, 1e-5);
	std::cout << " (Pressure)";
}

template class Projector<2>;
template class Projector<3>;
}
