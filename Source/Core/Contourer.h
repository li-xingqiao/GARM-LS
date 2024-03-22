#pragma once

#include "GridBasedData.h"
#include "SurfaceMesh.h"

#include <array>

namespace PhysX {

template <int Dim, bool Filled = false>
class Contourer
{
	DECLARE_DIM_TYPES(Dim)
	static_assert(Dim == 2 || !Filled, "Dimension must be 2 if filled.");

protected:

#include "MarchingCubesTables.inc"

public:

	static SurfaceMesh<Dim> contour(const GridBasedData<Dim, double> &gbd, const double value = 0);

protected:

	static auto makeEdgeGrids(const Grid<Dim> &nodeGrid)
	{
		if constexpr (Dim == 2) {
			return std::array {
				Grid<2>(nodeGrid.spacing, nodeGrid.size - VectorDi::Unit(0), nodeGrid.origin + VectorDd::Unit(0) * nodeGrid.spacing / 2),
				Grid<2>(nodeGrid.spacing, nodeGrid.size - VectorDi::Unit(1), nodeGrid.origin + VectorDd::Unit(1) * nodeGrid.spacing / 2)
			};
		}
		else {
			return std::array {
				Grid<3>(nodeGrid.spacing, nodeGrid.size - VectorDi::Unit(0), nodeGrid.origin + VectorDd::Unit(0) * nodeGrid.spacing / 2),
				Grid<3>(nodeGrid.spacing, nodeGrid.size - VectorDi::Unit(1), nodeGrid.origin + VectorDd::Unit(1) * nodeGrid.spacing / 2),
				Grid<3>(nodeGrid.spacing, nodeGrid.size - VectorDi::Unit(2), nodeGrid.origin + VectorDd::Unit(2) * nodeGrid.spacing / 2)
			};
		}
	}

	static auto makeEdgeMark(const std::array<Grid<Dim>, Dim> &edgeGrids)
	{
		if constexpr (Dim == 2) {
			return std::array {
				GridBasedData<2, int>(edgeGrids[0], -1),
				GridBasedData<2, int>(edgeGrids[1], -1)
			};
		}
		else {
			return std::array {
				GridBasedData<3, int>(edgeGrids[0], -1),
				GridBasedData<3, int>(edgeGrids[1], -1),
				GridBasedData<3, int>(edgeGrids[2], -1)
			};
		}
	}

	static uint getCellType(const GridBasedData<Dim, double> &gbd, const VectorDi &cell, const double value);
};

}
