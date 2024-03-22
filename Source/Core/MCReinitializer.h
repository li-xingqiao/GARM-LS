#pragma once

#include "Contourer.h"
#include "MeshToSdf.h"
#include "Refine.h"

namespace PhysX {

template <int Dim>
class MCReinitializer
{
	DECLARE_DIM_TYPES(Dim)

public:
	MCReinitializer() { }

	void perform(GridBasedData<Dim, double> &phi, const int maxSteps)
	{
		SurfaceMesh<Dim> liquidMesh = Contourer<Dim, Dim == 2>::contour(phi);
		liquidMesh.computeFaceNormals();
		convertMeshToSdf(liquidMesh, phi, maxSteps);
	}

	template <typename Field>
	void perform(GridBasedData<Dim, double> &phi, const Field &phiView, const int maxSteps)
	{
		GridBasedData<Dim, double> refinedLevelSet(Refine::get(phi.grid, 4));
		Refine::fill(refinedLevelSet, phiView);
		SurfaceMesh<Dim> liquidMesh = Contourer<Dim, Dim == 2>::contour(refinedLevelSet);
		liquidMesh.computeFaceNormals();
		convertMeshToSdf(liquidMesh, phi, maxSteps);
	}
};

}
