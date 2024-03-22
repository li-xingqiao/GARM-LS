#include "Boundary.h"

#include "VectorField.h"

namespace PhysX {

template <int Dim>
Boundary<Dim>::Boundary(const StaggeredGrid<Dim> &sGrid) :
	domainBox(sGrid.domainOrigin(), sGrid.domainLengths()),
	nodeDist(sGrid.nodeGrid),
	cellDist(sGrid.cellGrid),
	fraction(sGrid.faceGrids),
	velocity(sGrid.faceGrids),
	normal(sGrid.faceGrids)
{
	reset();
}

template <int Dim>
void Boundary<Dim>::reset()
{
	parallelForEach(nodeDist.grid, [&](const VectorDi &node) {
		const VectorDd pos = nodeDist.grid.position(node);
		nodeDist[node] = -domainBox.signedDistance(pos);
	});
	parallelForEach(cellDist.grid, [&](const VectorDi &cell) {
		const VectorDd pos = cellDist.grid.position(cell);
		cellDist[cell] = -domainBox.signedDistance(pos);
	});
	parallelForEach(fraction.grids, [&](const int axis, const VectorDi &face) {
		const VectorDd pos = fraction[axis].grid.position(face);
		fraction[axis][face] = -domainBox.signedDistance(pos); // fraction as an sdf currently
		normal[axis][face] = -domainBox.closestNormal(pos)[axis];
	});
}

template <int Dim>
void Boundary<Dim>::unions(const Surface<Dim> &surface)
{
	parallelForEach(nodeDist.grid, [&](const VectorDi &node) {
		const VectorDd pos = nodeDist.grid.position(node);
		nodeDist[node] = std::min(nodeDist[node], surface.signedDistance(pos));
	});
	parallelForEach(cellDist.grid, [&](const VectorDi &cell) {
		const VectorDd pos = cellDist.grid.position(cell);
		cellDist[cell] = std::min(cellDist[cell], surface.signedDistance(pos));
	});
	parallelForEach(fraction.grids, [&](const int axis, const VectorDi &face) {
		const VectorDd pos = fraction[axis].grid.position(face);
		const double phi = surface.signedDistance(pos);
		if (phi < fraction[axis][face]) { // fraction as an sdf currently
			fraction[axis][face] = phi;
			normal[axis][face] = surface.closestNormal(pos)[axis];
		}
	});
}

template <int Dim>
void Boundary<Dim>::intersects(const Surface<Dim> &surface)
{
	parallelForEach(nodeDist.grid, [&](const VectorDi &node) {
		const VectorDd pos = nodeDist.grid.position(node);
		nodeDist[node] = std::max(nodeDist[node], surface.signedDistance(pos));
	});
	parallelForEach(cellDist.grid, [&](const VectorDi &cell) {
		const VectorDd pos = cellDist.grid.position(cell);
		cellDist[cell] = std::max(cellDist[cell], surface.signedDistance(pos));
	});
	parallelForEach(fraction.grids, [&](const int axis, const VectorDi &face) {
		const VectorDd pos = fraction[axis].grid.position(face);
		const double phi = surface.signedDistance(pos);
		if (phi > fraction[axis][face]) { // fraction as an sdf currently
			fraction[axis][face] = phi;
			normal[axis][face] = surface.closestNormal(pos)[axis];
		}
	});
}

template <int Dim>
void Boundary<Dim>::finish(const StaggeredGrid<Dim> &sGrid)
{
	// Convert fraction to be boundary fraction.
	parallelForEach(fraction.grids, [&](const int axis, const VectorDi &face) {
		if (sGrid.isBoundaryFace(axis, face)) fraction[axis][face] = 1;
		else fraction[axis][face] = getFaceFraction(axis, face);
	});
}

template <int Dim>
void Boundary<Dim>::enforce(SGridBasedData<Dim, double> &fluidVelocity) const
{
	const auto normalView = VectorFieldView<LinearIntrpl<Dim>>(normal);
	const auto velocityView = VectorFieldView<LinearIntrpl<Dim>>(velocity);
	const auto fluidVelocityView = VectorFieldView<LinearIntrpl<Dim>>(fluidVelocity);

	auto newFluidVelocity = fluidVelocity;
	parallelForEach(fluidVelocity.grids, [&](const int axis, const VectorDi &face) {
		if (fraction[axis][face] > 0.1 && fraction[axis][face] <= 1) {
			const VectorDd pos = fluidVelocity.grids[axis].position(face);
			const VectorDd n = normalView(pos).normalized();
			if (n.any()) {
				const VectorDd vel = fluidVelocityView(pos);
				const VectorDd vel0 = velocityView(pos);
				const double vn = (vel - vel0).dot(n);
				const VectorDd vt = (vel - vel0) - (vel - vel0).dot(n) * n;
				if (vn < 0) {
					newFluidVelocity[axis][face] = (vel - (vel - vel0).dot(n) * n)[axis] - 0.8 * vt[axis];
				} else {
					newFluidVelocity[axis][face] = fluidVelocity[axis][face];
				}
			}
			else
				newFluidVelocity[axis][face] = fluidVelocity[axis][face];
		}
		else if (fraction[axis][face] == 1) {
			newFluidVelocity[axis][face] = velocity[axis][face];
		}
	});
	fluidVelocity = newFluidVelocity;
}

template <int Dim>
void Boundary<Dim>::enforce(GridBasedData<Dim, Vector<Dim, double>> &fluidRefMap) const
{
	auto newFluidRefMap = fluidRefMap;
	parallelForEach(fluidRefMap.grid, [&](const VectorDi &cell) {
		if (cellDist[cell] < 0) {
			const VectorDd pos = fluidRefMap.grid.position(cell);
			newFluidRefMap[cell] = pos;
		}
	});
	fluidRefMap = newFluidRefMap;
}

template <int Dim>
double Boundary<Dim>::getFaceFraction(const int axis, const VectorDi &face) const
{
	const auto theta = [](const double phi0, const double phi1)->double { return phi0 / (phi0 - phi1); };

	double faceFraction;

	if constexpr (Dim == 2) {
		const auto fraction = [&theta](const double phi0, const double phi1)->double {
			if (phi0 <= 0 && phi1 <= 0) return 1;
			else if (phi0 <= 0 && phi1 > 0) return theta(phi0, phi1);
			else if (phi0 > 0 && phi1 <= 0) return theta(phi1, phi0);
			else return 0;
		};
		const double phi0 = nodeDist[StaggeredGrid<Dim>::faceNode(axis, face, 0)];
		const double phi1 = nodeDist[StaggeredGrid<Dim>::faceNode(axis, face, 1)];
		faceFraction = fraction(phi0, phi1);
	}
	else {
		const auto fraction = [&](double phi0, double phi1, double phi2)->double {
			if (phi0 > phi1) std::swap(phi0, phi1);
			if (phi1 > phi2) std::swap(phi1, phi2);
			if (phi0 > phi1) std::swap(phi0, phi1);
			if (phi2 <= 0) return 1;
			else if (phi1 <= 0) return 1 - theta(phi2, phi0) * theta(phi2, phi1);
			else if (phi0 <= 0) return theta(phi0, phi1) * theta(phi0, phi2);
			else return 0;
		};
		const double phi0 = nodeDist[StaggeredGrid<Dim>::faceNode(axis, face, 0)];
		const double phi1 = nodeDist[StaggeredGrid<Dim>::faceNode(axis, face, 1)];
		const double phi2 = nodeDist[StaggeredGrid<Dim>::faceNode(axis, face, 2)];
		const double phi3 = nodeDist[StaggeredGrid<Dim>::faceNode(axis, face, 3)];
		const double centerPhi = (phi0 + phi1 + phi2 + phi3) * .25;
		faceFraction = (fraction(centerPhi, phi0, phi1) + fraction(centerPhi, phi1, phi3) + fraction(centerPhi, phi3, phi2) + fraction(centerPhi, phi2, phi0)) * .25;
	}

	return faceFraction > .9 ? 1 : faceFraction;
	// return faceFraction > .1 ? 1 : 0;
	// return faceFraction > .9 ? 1 : (faceFraction < .1 ? 0 : 0.5);
}

template class Boundary<2>;
template class Boundary<3>;

}
