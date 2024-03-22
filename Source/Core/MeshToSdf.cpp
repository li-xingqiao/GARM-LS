#include "MeshToSdf.h"
#include "FastMarching.h"

namespace PhysX {

static constexpr double MTS_EPS = 1e-10;

template <int Dim>
static double getSegmentDistance(const Vector<Dim, double> &p, const Vector<Dim, double> &a, const Vector<Dim, double> &b)
{
	double pmap = (p - a).dot(p - b);
	double length = (b - a).norm();
	if (pmap < 0) return (p - a).norm();
	else if (pmap > length) return (p - b).norm();
	else return (p - a - pmap * (b - a) / length).norm();
}

template <int Dim>
static bool computeBarycentric(const Vector<Dim, double> &p, const Vector<Dim, double> &a, const Vector<Dim, double> &b, const Vector<Dim, double> &c, Vector3d &result)
{
	Vector<Dim, double> n = (b - a).cross(c - a);
	double nn = n.squaredNorm();
	if (nn < MTS_EPS * MTS_EPS) {
		return false; // degenerated
	}
	Vector<Dim, double> na = (p - b).cross(p - c);
	Vector<Dim, double> nb = (p - c).cross(p - a);
	Vector<Dim, double> nc = (p - a).cross(p - b);
	result(0) = n.dot(na) / nn;
	result(1) = n.dot(nb) / nn;
	result(2) = n.dot(nc) / nn;
	double sum = result(0) + result(1) + result(2);
	result /= sum;
	return true;
}

template <int Dim>
static double getTriangleDistance(const Vector<Dim, double> &p, const Vector<Dim, double> &a, const Vector<Dim, double> &b, const Vector<Dim, double> &c)
{
	Vector<Dim, double> ab(b - a), bc(c - b), ca(a - c);
	Vector<Dim, double> nab((c - a) - (c - a).dot(ab) * ab / ab.squaredNorm()), nbc((a - b) - (a - b).dot(bc) * bc / bc.squaredNorm()), nca((b - c) - (b - c).dot(ca) * ca / ca.squaredNorm());
	if (nab.norm() < MTS_EPS) {
		// degenerated
		if (ab.dot(bc) >= 0) {
			return getSegmentDistance(p, a, c);
		}
		else if (bc.dot(ca) >= 0) {
			return getSegmentDistance(p, b, a);
		}
		else {
			return getSegmentDistance(p, c, b);
		}
	}
	if ((p - a).dot(nab) < 0) {
		return getSegmentDistance(p, a, b);
	}
	else if ((p - b).dot(nbc) < 0) {
		return getSegmentDistance(p, b, c);
	}
	else if ((p - c).dot(nca) < 0) {
		return getSegmentDistance(p, c, a);
	}
	else {
		if constexpr (Dim == 2) {
			return 0.0;
		}
		else {
			bc = ab.cross(ca); // bc now represent the normal
			return std::abs((p - a).dot(bc)) / bc.norm();
		}
	}
}

// This implementation can't deal with p==a / p on ab circumstances
// static bool intersectTriangle3(const Vector3d &p, const Vector3d &a, const Vector3d &b, const Vector3d &c, const int axis, Vector3d &result) {
// 	Vector3d proj_p(p), proj_a(a), proj_b(b), proj_c(c);
// 	proj_p(axis) = 0.0, proj_a(axis) = 0.0, proj_b(axis) = 0.0, proj_c(axis) = 0.0;
// 	Vector3d bcCoord;
// 	if (!computeBarycentric(proj_p, proj_a, proj_b, proj_c, bcCoord)) return false;
// 	if (bcCoord(0) < 0.0 || bcCoord(1) < 0.0 || bcCoord(2) < 0.0) return false;
// 	result = proj_p;
// 	result(axis) = bcCoord(0) * a(axis) + bcCoord(1) * b(axis) + bcCoord(2) * c(axis);
// 	return true;
// }

static int orientation(double x1, double y1, double x2, double y2, double &twice_signed_area)
{
	twice_signed_area = y1 * x2 - x1 * y2;
	if (twice_signed_area > 0) return 1;
	else if (twice_signed_area < 0) return -1;
	else if (y2 > y1) return 1;
	else if (y2 < y1) return -1;
	else if (x1 > x2) return 1;
	else if (x1 < x2) return -1;
	else return 0; // only true when x1==x2 and y1==y2
}

static bool intersectTriangle3(const Vector3d &p, const Vector3d &a, const Vector3d &b, const Vector3d &c, const int axis, Vector3d &result)
{
	int ax1 = (axis + 1) % 3;
	int ax2 = (axis + 2) % 3;
	double ra, rb, rc;
	int signa = orientation(p(ax1) - b(ax1), p(ax2) - b(ax2), p(ax1) - c(ax1), p(ax2) - c(ax2), ra);
	if (signa == 0) return false;
	int signb = orientation(p(ax1) - c(ax1), p(ax2) - c(ax2), p(ax1) - a(ax1), p(ax2) - a(ax2), rb);
	if (signb != signa) return false;
	int signc = orientation(p(ax1) - a(ax1), p(ax2) - a(ax2), p(ax1) - b(ax1), p(ax2) - b(ax2), rc);
	if (signc != signa) return false;
	double sum = ra + rb + rc;
	assert(sum != 0.0);
	result = p;
	result(axis) = (ra * a(axis) + rb * b(axis) + rc * c(axis)) / sum;
	return true;
}

template <int Dim>
static void spreadSign(GridBasedData<Dim, int> &sign)
{
	DECLARE_DIM_TYPES(Dim)
	std::queue<VectorDi> queue;
	forEach(sign.grid, [&](const VectorDi &cell) {
		if (sign[cell] != 0)
			queue.push(cell);
	});
	VectorDi coord;
	while (!queue.empty()) {
		coord = queue.front();
		queue.pop();
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
			if (sign.grid.isValid(nbCoord) && sign[nbCoord] == 0) {
				sign[nbCoord] = sign[coord];
				queue.push(nbCoord);
			}
		}
	}
}

template <int Dim>
void convertMeshToSdf(const SurfaceMesh<Dim> &mesh, GridBasedData<Dim, double> &sdf, const int maxSteps, const int bandWidth)
{
	DECLARE_DIM_TYPES(Dim)

	// GridBasedData<Dim, int> closestTriangle(sdf.grid);
	GridBasedData<Dim, int> outsideMesh(sdf.grid);
	GridBasedData<Dim, int> trianglesUpside(sdf.grid); // along axis 2
	parallelForEach(sdf.grid, [&](const VectorDi &cell) {
		sdf[cell] = 1e9;
	});

	for (int i = 0; i < mesh.indices.size(); i += Dim) {
		if constexpr (Dim == 2) {
			const VectorDd a = mesh.positions[mesh.indices[i]];
			const VectorDd b = mesh.positions[mesh.indices[i + 1]];

			VectorDd lCorner = a.cwiseMin(b);
			VectorDd rCorner = a.cwiseMax(b);
			VectorDi minCoord = sdf.grid.getLinearLower(lCorner) - bandWidth * VectorDi::Ones();
			VectorDi maxCoord = sdf.grid.getLinearLower(rCorner) + (bandWidth + 1) * VectorDi::Ones();
			minCoord = sdf.grid.clamp(minCoord);
			maxCoord = sdf.grid.clamp(maxCoord);

			for (int x0 = minCoord[0]; x0 <= maxCoord[0]; ++x0) {
				for (int x1 = minCoord[1]; x1 <= maxCoord[1]; ++x1) {
					VectorDi cell(x0, x1);
					VectorDd cellCenter(sdf.grid.position(cell));
					double dist = getSegmentDistance(cellCenter, a, b);
					if (dist < sdf[cell] + MTS_EPS) {
						sdf[cell] = dist;
						// closestTriangle[cell] = i;
						if (mesh.faceNormals[i / Dim].dot(cellCenter - a) >= 0) outsideMesh[cell] = 1;
						else if (outsideMesh[cell] == 0) outsideMesh[cell] = -1;
					}
				}
				Vector2d intersection;
				Vector2i iterCoord = Vector2i(x0, 0);
				intersection = sdf.grid.position(iterCoord);
				if (std::abs(a(0) - b(0)) < MTS_EPS) continue;
				// double c = std::max(a(0), b(0));
				if ((a(0) - intersection(0)) * (intersection(0) - b(0)) > 0 && std::abs(b(0) - intersection(0)) > MTS_EPS) {
					intersection(1) = ((a(0) - intersection(0)) * b(1) + (intersection(0) - b(0)) * a(1)) / (a(0) - b(0));
					iterCoord(1) = sdf.grid.getLinearLower(intersection)(1) + 1;
					while (iterCoord(1) < sdf.grid.size(1)) {
						++trianglesUpside[iterCoord];
						iterCoord(1) += 1;
					}
				}
			}
		}
		else if constexpr (Dim == 3) {
			const VectorDd a = mesh.positions[mesh.indices[i]];
			const VectorDd b = mesh.positions[mesh.indices[i + 1]];
			const VectorDd c = mesh.positions[mesh.indices[i + 2]];

			VectorDd lCorner = a.cwiseMin(b).cwiseMin(c);
			VectorDd rCorner = a.cwiseMax(b).cwiseMax(c);
			VectorDi minCoord = sdf.grid.getLinearLower(lCorner) - bandWidth * VectorDi::Ones();
			VectorDi maxCoord = sdf.grid.getLinearLower(rCorner) + (bandWidth + 1) * VectorDi::Ones();
			minCoord = sdf.grid.clamp(minCoord);
			maxCoord = sdf.grid.clamp(maxCoord);

			for (int x0 = minCoord[0]; x0 <= maxCoord[0]; ++x0)
			for (int x1 = minCoord[1]; x1 <= maxCoord[1]; ++x1) {
				for (int x2 = minCoord[2]; x2 <= maxCoord[2]; ++x2) {
					VectorDi cell(x0, x1, x2);
					VectorDd cellCenter(sdf.grid.position(cell));
					double dist = getTriangleDistance(cellCenter, a, b, c);
					if (dist < sdf[cell] + MTS_EPS) {
						sdf[cell] = dist;
						// closestTriangle[cell] = i;
						if (mesh.faceNormals[i / Dim].dot(cellCenter - a) >= 0) outsideMesh[cell] = 1;
						else if (outsideMesh[cell] == 0) outsideMesh[cell] = -1;
					}
				}
				Vector3d intersection;
				Vector3i iterCoord(x0, x1, 0);
				if (intersectTriangle3(sdf.grid.position(iterCoord), a, b, c, /*axis=*/2, intersection)) {
					iterCoord(2) = sdf.grid.getLinearLower(intersection)(2) + 1;
					while (iterCoord(2) < sdf.grid.size(2)) {
						if (iterCoord(2) > 0)
							++trianglesUpside[iterCoord];
						iterCoord(2) += 1;
					}
				}
			}
		}
	}

	// spreadSign(outsideMesh);
	parallelForEach(sdf.grid, [&](const VectorDi &cell) {
		// sdf[cell] = sdf[cell] * outsideMesh[cell];
		sdf[cell] = sdf[cell] * ((trianglesUpside[cell] % 2 == 0) ? 1 : -1);
	});

	FastMarching<Dim> fmm(sdf.grid);
	fmm.perform(sdf, maxSteps);
}

template void convertMeshToSdf<2>(const SurfaceMesh<2> &mesh, GridBasedData<2, double> &sdf, const int maxSteps, const int bandWidth);
template void convertMeshToSdf<3>(const SurfaceMesh<3> &mesh, GridBasedData<3, double> &sdf, const int maxSteps, const int bandWidth);

}
