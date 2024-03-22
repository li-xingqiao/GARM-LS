#include "SurfaceMesh.h"

#include <algorithm>
#include <format>
#include <limits>
#include <sstream>

namespace PhysX {

template<int Dim>
Vector<Dim, double> SurfaceMesh<Dim>::closestPosition(const VectorDd& pos) const
{
	VectorDd normal;
	return closestPosition(pos, normal);
}

template<int Dim>
Vector<Dim, double> SurfaceMesh<Dim>::closestNormal(const VectorDd& pos) const
{
	VectorDd normal;
	const VectorDd closestPos = closestPosition(pos, normal);
	return normal;
}

template <int Dim>
double SurfaceMesh<Dim>::signedDistance(const VectorDd& pos) const
{
	VectorDd normal;
	const VectorDd closestPos = closestPosition(pos, normal);
	return (pos - closestPos).norm() * (normal.dot(pos - closestPos) <= 0 ? -1 : 1);
}

template<int Dim>
bool SurfaceMesh<Dim>::isInside(const VectorDd& pos) const
{
	VectorDd normal;
	const VectorDd closestPos = closestPosition(pos, normal);
	return normal.dot(pos - closestPos) <= 0;
}

template <int Dim>
Vector<Dim, double> SurfaceMesh<Dim>::closestPosition(const VectorDd& pos, VectorDd& normal) const
{
	VectorDd closestPos;
	double minDist = std::numeric_limits<double>::infinity();

	for (size_t i = 0; i < positions.size(); i++) {
		const double tempDist = (positions[i] - pos).norm();
		if (tempDist < minDist) {
			minDist = tempDist;
			closestPos = positions[i];
			normal = normals[i];
		}
	}

	for (size_t i = 0; i < indices.size(); i += Dim) {
		const VectorDd n = faceNormals[i / Dim];
		if constexpr (Dim == 2) {
			const VectorDd a = positions[indices[i]];
			const VectorDd b = positions[indices[i + 1]];

			const double tempDist = (pos - a).dot(n);
			const VectorDd p = pos - tempDist * n;

			if ((p - a).dot(p - b) <= 0 && std::abs(tempDist) < minDist) {
				minDist = std::abs(tempDist);
				closestPos = p;
				normal = n;
			}
		}
		else {
			const VectorDd a = positions[indices[i]];
			const VectorDd b = positions[indices[i + 1]];
			const VectorDd c = positions[indices[i + 2]];

			const double tempDist = (pos - a).dot(n);
			const VectorDd p = pos - tempDist * n;

			const VectorDd v0 = c - a;
			const VectorDd v1 = b - a;
			const VectorDd v2 = p - a;

			const double dot00 = v0.dot(v0);
			const double dot01 = v0.dot(v1);
			const double dot02 = v0.dot(v2);
			const double dot11 = v1.dot(v1);
			const double dot12 = v1.dot(v2);

			const double deno = 1 / (dot00 * dot11 - dot01 * dot01);

			const double u = (dot11 * dot02 - dot01 * dot12) * deno;
			const double v = (dot00 * dot12 - dot01 * dot02) * deno;
			if (0 <= u && 0 <= v && u + v <= 1 && std::abs(tempDist) < minDist) {
				minDist = std::abs(tempDist);
				closestPos = p;
				normal = n;
			}
		}
	}

	return closestPos;
}

template <int Dim>
inline void SurfaceMesh<Dim>::write(std::ostream &out) const
{
	IO::writeValue(out, uint(positions.size()));
	for (const auto &pos : positions)
		IO::writeValue(out, pos.template cast<float>().eval());
	if constexpr (Dim == 3) {
		if (normals.size() != positions.size()) {
			std::vector<VectorDd> normals_;
			computeNormals(normals_);
			for (const auto &normal : normals_)
				IO::writeValue(out, normal.template cast<float>().eval());
		}
		else {
			for (const auto &normal : normals)
				IO::writeValue(out, normal.template cast<float>().eval());
		}
	}
	IO::writeValue(out, uint(indices.size()));
	IO::writeArray(out, indices.data(), indices.size());
}

template <int Dim>
void SurfaceMesh<Dim>::writeOBJ(const std::string &fileName) const
{
	std::ofstream fout(fileName);
	if constexpr (Dim == 2) {
		for (const auto &pos : positions)
			fout << std::format("v {:.15e} {:.15e} 0", pos.x(), pos.y()) << std::endl;
		for (size_t i = 0; i < indices.size(); i += 3)
			fout << std::format("f {} {} {}", indices[i] + 1, indices[i + 1] + 1, indices[i + 2] + 1) << std::endl;
	}
	else {
		for (const auto &pos : positions)
			fout << std::format("v {:.15e} {:.15e} {:.15e}", pos.x(), pos.y(), pos.z()) << std::endl;
		for (const auto &normal : normals)
			fout << std::format("vn {:.15e} {:.15e} {:.15e}", normal.x(), normal.y(), normal.z()) << std::endl;
		for (size_t i = 0; i < indices.size(); i += 3)
			fout << std::format("f {} {} {}", indices[i] + 1, indices[i + 1] + 1, indices[i + 2] + 1) << std::endl;
	}
}

template <int Dim>
void SurfaceMesh<Dim>::loadOBJ(const std::string& fileName)
{
	clear();

	std::ifstream fin(fileName);
	std::cout << std::format("** Loading {}... ", fileName) << std::flush;

	std::string line;
	while (std::getline(fin, line)) {
		if (line.substr(0, 2) == "v ") {
			std::istringstream iss(line.substr(2));
			VectorDd pos;
			iss >> pos[0] >> pos[1];
			if constexpr (Dim == 3) iss >> pos[2];
			positions.push_back(pos);
		}
		else if (line.substr(0, 2) == "f ") {
			std::istringstream iss(line.substr(2));
			std::vector<std::string> vtx(3);
			iss >> vtx[0];
			iss >> vtx[1];
			int i = 0;
			while (iss >> vtx[2]) {
				for (const std::string &s : vtx) {
					indices.push_back(uint(std::stoi(s.substr(0, s.find_first_of('/'))) - 1));
				}
				vtx[1] = vtx[2];
			}
		}
	}

	std::cout << std::format("{} vertices, {} {}.", positions.size(), indices.size() / Dim, Dim == 2 ? "lines" : "faces") << std::endl;
}

template <int Dim>
void SurfaceMesh<Dim>::computeFaceNormals(std::vector<VectorDd>& faceNormals_) const
{
	const size_t numFaces = indices.size() / Dim;
	faceNormals_.resize(numFaces);
	for (size_t i = 0; i < indices.size(); i += Dim)
		faceNormals_[i / Dim] = getFaceNormal(i);
}

template <int Dim>
void SurfaceMesh<Dim>::computeNormals(std::vector<VectorDd> &normals_) const
{
	normals_.resize(positions.size());
	std::fill(normals_.begin(), normals_.end(), VectorDd::Zero());

	const bool dirty = faceNormals.size() * Dim != indices.size();

	for (size_t i = 0; i < indices.size(); i += Dim) {
		const VectorDd n = dirty ? getFaceNormal(i) : faceNormals[i / Dim];
		for (size_t j = i; j < i + Dim; j++) normals_[indices[j]] += n;
	}
	for (auto &normal : normals_) normal.normalize();
}

template <int Dim>
Vector<Dim, double> SurfaceMesh<Dim>::getFaceNormal(const size_t i) const
{
	if constexpr (Dim == 2) {
		const VectorDd a = positions[indices[i + 1]] - positions[indices[i]];
		return VectorDd(a[1], -a[0]).normalized();
	}
	else {
		const VectorDd a = positions[indices[i + 1]] - positions[indices[i]];
		const VectorDd b = positions[indices[i + 2]] - positions[indices[i]];
		return a.cross(b).normalized();
	}
}

template class SurfaceMesh<2>;
template class SurfaceMesh<3>;

}
