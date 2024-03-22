#include "Contourer.h"

#include "StaggeredGrid.h"

#include <memory>

namespace PhysX {

template <int Dim, bool Filled>
SurfaceMesh<Dim> Contourer<Dim, Filled>::contour(const GridBasedData<Dim, double> &gbd, const double value)
{
	SurfaceMesh<Dim> mesh;
	// Initialize auxiliary grids and data structures.
	const Grid<Dim> &nodeGrid = gbd.grid;
	const Grid<Dim> cellGrid(nodeGrid.spacing, nodeGrid.size - VectorDi::Ones(), nodeGrid.origin + VectorDd::Ones() * nodeGrid.spacing / 2);
	const auto edgeGrids = makeEdgeGrids(nodeGrid);
	auto edgeMark = makeEdgeMark(edgeGrids);
	// Initialize for filling.
	std::unique_ptr<GridBasedData<Dim, int>> nodeMarkPtr;
	if constexpr (Filled) {
		nodeMarkPtr = std::make_unique<GridBasedData<Dim, int>>(nodeGrid, -1);
		forEach(nodeGrid, [&](const VectorDi &node) {
			if (gbd[node] <= value) {
				nodeMarkPtr->operator[](node) = int(mesh.positions.size());
				mesh.positions.push_back(nodeGrid.position(node));
			}
		});
	}
	// Perform the marching-cubes algorithm.
	forEach(cellGrid, [&](const VectorDi &cell) {
		const uint cellType = getCellType(gbd, cell, value);
		uint edgeState;
		if constexpr (Dim == 2) edgeState = _kEdgeStateTable2[cellType];
		else edgeState = _kEdgeStateTable3[cellType];
		for (int i = 0; i < StaggeredGrid<Dim>::numberOfCellEdges(); i++) {
			if (edgeState >> i & 1) {
				const int axis = StaggeredGrid<Dim>::cellEdgeAxis(i);
				const VectorDi edge = StaggeredGrid<Dim>::cellEdge(cell, i);
				if (edgeMark[axis][edge] < 0) {
					const VectorDi node0 = StaggeredGrid<Dim>::edgeAdjacentNode(axis, edge, 0);
					const VectorDi node1 = StaggeredGrid<Dim>::edgeAdjacentNode(axis, edge, 1);
					const double theta = (gbd[node0] - value) / (gbd[node0] - gbd[node1]);
					const VectorDd pos = (1 - theta) * nodeGrid.position(node0) + theta * nodeGrid.position(node1);
					edgeMark[axis][edge] = int(mesh.positions.size());
					mesh.positions.push_back(pos);
				}
			}
		}
		if constexpr (!Filled) {
			const int *edgeOrds;
			if constexpr (Dim == 2) edgeOrds = _kEdgeOrdsTable2[cellType];
			else edgeOrds = _kEdgeOrdsTable3[cellType];
			for (int i = 0; edgeOrds[i] != -1; i++) {
				const int axis = StaggeredGrid<Dim>::cellEdgeAxis(edgeOrds[i]);
				const VectorDi edge = StaggeredGrid<Dim>::cellEdge(cell, edgeOrds[i]);
				mesh.indices.push_back(uint(edgeMark[axis][edge]));
			}
		}
		else {
			const int *vertexOrds = _kVertexOrdsTable2[cellType];
			for (int i = 0; vertexOrds[i] != -1; i++) {
				if (vertexOrds[i] < 4) {
					const int axis = StaggeredGrid<Dim>::cellEdgeAxis(vertexOrds[i]);
					const VectorDi edge = StaggeredGrid<Dim>::cellEdge(cell, vertexOrds[i]);
					mesh.indices.push_back(uint(edgeMark[axis][edge]));
				}
				else {
					const VectorDi node = StaggeredGrid<Dim>::cellNode(cell, vertexOrds[i] & 3);
					mesh.indices.push_back(uint(nodeMarkPtr->operator[](node)));
				}
			}
		}
	});
	return mesh;
}

template <int Dim, bool Filled>
uint Contourer<Dim, Filled>::getCellType(const GridBasedData<Dim, double> &gbd, const VectorDi &cell, const double value)
{
	uint type = 0;
	for (int i = 0; i < StaggeredGrid<Dim>::numberOfCellNodes(); i++) {
		const VectorDi node = StaggeredGrid<Dim>::cellNode(cell, i);
		if (gbd[node] <= value) type |= 1 << i;
	}
	// Handle diagonal cases for 2D.
	if constexpr (Dim == 2) {
		if (type == 6 || type == 9) {
			const double phi0 = gbd[StaggeredGrid<Dim>::cellNode(cell, 0)];
			const double phi1 = gbd[StaggeredGrid<Dim>::cellNode(cell, 1)];
			const double phi2 = gbd[StaggeredGrid<Dim>::cellNode(cell, 2)];
			const double phi3 = gbd[StaggeredGrid<Dim>::cellNode(cell, 3)];
			const double centerPhi = (phi0 + phi1 + phi2 + phi3) / 4;
			if (centerPhi <= value) type = 16 + (type > 8);
		}
	}
	return type;
}

template class Contourer<2>;
template class Contourer<2, true>;
template class Contourer<3>;

}
