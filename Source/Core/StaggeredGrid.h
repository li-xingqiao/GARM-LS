#pragma once

#include "Grid.h"

namespace PhysX {

template <int Dim>
class StaggeredGrid
{
	DECLARE_DIM_TYPES(Dim)

public:

	const int boundaryWidth;

	const double spacing;
	const double invSpacing;
	const VectorDi resolution;
	const VectorDd origin;

	const Grid<Dim> nodeGrid;
	const Grid<Dim> cellGrid;
	const std::array<Grid<Dim>, Dim> faceGrids;

public:

	StaggeredGrid(const int boundaryWidth_, const double spacing_, const VectorDi &resolution_, const VectorDd &center_ = VectorDd::Zero()) :
		boundaryWidth(boundaryWidth_),
		spacing(spacing_),
		invSpacing(1 / spacing),
		resolution(resolution_),
		origin(center_ - resolution.template cast<double>() * spacing / 2),
		nodeGrid(spacing, resolution + VectorDi::Ones(), origin),
		cellGrid(spacing, resolution, origin + VectorDd::Ones() * spacing / 2),
		faceGrids(makeFaceGrids(spacing, resolution, origin))
	{ }

	bool operator==(const StaggeredGrid &other) const { return boundaryWidth == other.boundaryWidth && spacing == other.spacing && resolution == other.resolution && origin == other.origin; }

	VectorDd domainOrigin() const { return origin + VectorDd::Ones() * boundaryWidth * spacing; }
	VectorDd domainLengths() const { return (resolution - VectorDi::Ones() * boundaryWidth * 2).template cast<double>() * spacing; }
	double radius() const { return .5 * (Dim == 2 ? domainLengths().y() : domainLengths().norm()); }

	size_t nodeCount() const { return nodeGrid.count(); }
	size_t cellCount() const { return cellGrid.count(); }
	size_t faceCount(const int axis) const { return faceGrids[axis].count(); }

	VectorDd nodeCenter(const VectorDi &node) const { return nodeGrid.position(node); }
	VectorDd cellCenter(const VectorDi &cell) const { return cellGrid.position(cell); }
	VectorDd faceCenter(const int axis, const VectorDi &face) const { return faceGrids[axis].position(face); }

	bool isInsideNode(const VectorDi &node) const { return nodeGrid.isInside(node, boundaryWidth); }
	bool isBoundaryNode(const VectorDi &node) const { return !nodeGrid.isInside(node, boundaryWidth + 1); }
	bool isInsideCell(const VectorDi &cell) const { return cellGrid.isInside(cell, boundaryWidth); }
	bool isBoundaryCell(const VectorDi &cell) const { return !isInsideCell(cell); }
	bool isInsideFace(const int axis, const VectorDi &face) const { return faceGrids[axis].isInside(face, boundaryWidth); }
	bool isBoundaryFace(const int axis, const VectorDi &face) const { return face[axis] <= boundaryWidth || face[axis] >= resolution[axis] - boundaryWidth || !isInsideFace(axis, face); }

	static constexpr int numberOfCellNodes() { return 1 << Dim; }
	static constexpr int numberOfCellFaces() { return Dim << 1; }
	static constexpr int numberOfCellEdges() { return Dim * (1 << (Dim - 1)); }
	static constexpr int numberOfFaceNodes() { return 1 << (Dim - 1); }

	static VectorDi cellNode(const VectorDi &cell, const int ord)
	{
		if constexpr (Dim == 2) return cell + VectorDi(ord & 1, ord >> 1 & 1);
		else return cell + VectorDi(ord & 1, ord >> 1 & 1, ord >> 2 & 1);
	};

	static VectorDi nodeCell(const VectorDi &node, const int ord)
	{
		if constexpr (Dim == 2) return node - VectorDi((~ord) & 1, (~ord) >> 1 & 1);
		else return node - VectorDi((~ord) & 1, (~ord) >> 1 & 1, (~ord) >> 2 & 1);
	};

	static VectorDi faceNode(const int axis, const VectorDi &face, const int ord)
	{
		if constexpr (Dim == 2) return face + VectorDi::Unit(axis ^ 1) * (ord & 1);
		else return face + VectorDi::Unit((axis + 1) % 3) * (ord & 1) + VectorDi::Unit((axis + 2) % 3) * (ord >> 1 & 1);
	}

	static constexpr int cellFaceAxis(const int ord) { return ord >> 1; }
	static constexpr int cellFaceSide(const int ord) { return ord & 1 ? 1 : -1; }
	static VectorDi cellFace(const VectorDi &cell, const int ord) { return cell + VectorDi::Unit(cellFaceAxis(ord)) * (ord & 1); }
	static VectorDi faceAdjacentCell(const int axis, const VectorDi &face, const int ord) { return face - VectorDi::Unit(axis) * (ord & 1 ^ 1); }

	static constexpr int cellEdgeAxis(const int ord) { return ord >> (Dim - 1); }
	static VectorDi edgeAdjacentNode(const int axis, const VectorDi &edge, const int ord) { return edge + VectorDi::Unit(axis) * (ord & 1); }

	static VectorDi cellEdge(const VectorDi &cell, const int ord)
	{
		if constexpr (Dim == 2) return cell + VectorDi::Unit(cellEdgeAxis(ord) ^ 1) * (ord & 1);
		else return cell + VectorDi::Unit((cellEdgeAxis(ord) + 1) % 3) * (ord & 1) + VectorDi::Unit((cellEdgeAxis(ord) + 2) % 3) * (ord >> 1 & 1);
	}

	static VectorDi nodeFace(const VectorDi &node, const int ord)
	{
		if constexpr (Dim == 2) return node - VectorDi::Unit(cellEdgeAxis(ord) ^ 1) * ((ord & 1) ^ 1);
		else return node - VectorDi::Unit((cellEdgeAxis(ord) + 1) % 3) * ((ord & 1) ^ 1) - VectorDi::Unit((cellEdgeAxis(ord) + 2) % 3) * ((ord >> 1 & 1) ^ 1);
	}

private:

	static auto makeFaceGrids(const double spacing_, const VectorDi &resolution_, const VectorDd &origin_)
	{
		if constexpr (Dim == 2) {
			return std::array {
				Grid<Dim>(spacing_, resolution_ + VectorDi(1, 0), origin_ + VectorDd(0, 1) * spacing_ / 2),
				Grid<Dim>(spacing_, resolution_ + VectorDi(0, 1), origin_ + VectorDd(1, 0) * spacing_ / 2)
			};
		}
		else {
			return std::array {
				Grid<Dim>(spacing_, resolution_ + VectorDi(1, 0, 0), origin_ + VectorDd(0, 1, 1) * spacing_ / 2),
				Grid<Dim>(spacing_, resolution_ + VectorDi(0, 1, 0), origin_ + VectorDd(1, 0, 1) * spacing_ / 2),
				Grid<Dim>(spacing_, resolution_ + VectorDi(0, 0, 1), origin_ + VectorDd(1, 1, 0) * spacing_ / 2)
			};
		}
	}
};

}
