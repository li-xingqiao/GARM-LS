#pragma once

#include "SurfaceMesh.h"
#include "GridBasedData.h"

namespace PhysX {

template <int Dim>
void convertMeshToSdf(const SurfaceMesh<Dim> &mesh, GridBasedData<Dim, double> &sdf, const int maxSteps, const int bandWidth = 4);

}
