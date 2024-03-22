#pragma once

#include "IO.h"
#include "Surface.h"

#include <vector>

namespace PhysX {

template <int Dim>
class SurfaceMesh : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	std::vector<VectorDd> positions;
	std::vector<VectorDd> normals;
	std::vector<uint> indices;

	std::vector<VectorDd> faceNormals;

public:

	SurfaceMesh() = default;

	virtual VectorDd closestPosition(const VectorDd &pos) const override;
	virtual VectorDd closestNormal(const VectorDd &pos) const override;
	virtual double distance(const VectorDd &pos) const override { return (pos - closestPosition(pos)).norm(); }
	virtual double signedDistance(const VectorDd &pos) const override;
	virtual bool isInside(const VectorDd &pos) const override;

	VectorDd closestPosition(const VectorDd &pos, VectorDd &normal) const;

	void scale(const double factor) { for (auto &pos : positions) pos *= factor; }
	void translate(const VectorDd &delta) { for (auto &pos : positions) pos += delta; }

	void write(std::ostream &out) const;

	void writeOBJ(const std::string &fileName) const;

	void loadOBJ(const std::string &fileName);

	// TODO: Convert to SDF.
	// see https://github.com/AcademySoftwareFoundation/openvdb/blob/v9.0.0/openvdb/openvdb/tools/MeshToVolume.h#L3124

	void clear()
	{
		positions.clear();
		normals.clear();
		indices.clear();
		faceNormals.clear();
	}

	void computeFaceNormals(std::vector<VectorDd> &faceNormals_) const;
	void computeFaceNormals() { computeFaceNormals(faceNormals); }

	void computeNormals(std::vector<VectorDd> &normals_) const;
	void computeNormals() { computeNormals(normals); }

private:

	VectorDd getFaceNormal(const size_t i) const;
};

}
