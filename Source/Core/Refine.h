#pragma once

#include "Grid.h"

namespace PhysX::Refine {

template <int Dim>
inline Grid<Dim> get(const Grid<Dim> &original, int scale)
{
	return Grid<Dim>(original.spacing / scale, original.size * scale, original.origin);
}

template <int Dim, typename Data, typename Field>
void fill(GridBasedData<Dim, Data> &result, const Field &field)
{
	DECLARE_DIM_TYPES(Dim)

	parallelForEach(result.grid, [&](const VectorDi &cell) {
		const VectorDd pos = result.grid.position(cell);
		result[cell] = field(pos);
	});
}

}
