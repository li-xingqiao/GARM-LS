#include "LevelSetLiquid.h"
#include "LevelSetSurface.h"
#include "TransformedSurface.h"
#include "Contourer.h"

#include <numbers>

using namespace PhysX;

int main() {
	constexpr int Dim = 3;
	DECLARE_DIM_TYPES(Dim)

	const int scale = 256;
	const double length = 1;
	const int boundaryWidth = 2;
	VectorDi resolution(scale, scale, static_cast<int>(scale * 0.2));

	StaggeredGrid<Dim> sGrid(boundaryWidth, length / (scale - 2 * boundaryWidth), resolution);
	GridBasedData<Dim, double> boundarylevelset(sGrid.cellGrid);
	parallelForEach(sGrid.cellGrid, [&](const VectorDi &cell) {
		boundarylevelset[cell] = 1e5;
	});

	for (int i = 0; i < 5; ++i)
	for (int j = 0; j < 8 + (i % 2); ++j) {
		double h = 0.1 * (i - 1.5);
		double x = 0.1 * (j - 0.5 * (i % 2)) - 0.35;
		if constexpr (Dim == 2) {
			ImplicitSphere<Dim> rotatedCylinder(VectorDd(x, -h) * length, length * 0.035);
			unionSurfaceLevelSet(boundarylevelset, rotatedCylinder);
		} else {
			ImplicitCylinder<Dim> cylinder(VectorDd(x, -0.5, -h) * length, length * 0.035, length);
			RotatedSurface<Dim> rotatedCylinder(cylinder, VectorDd::Zero(), AngleAxisd(0.5 * std::numbers::pi, VectorDd::Unit(0)));
			unionSurfaceLevelSet(boundarylevelset, rotatedCylinder);
		}
	}

	const auto liquidMesh = Contourer<Dim, Dim == 2>::contour(boundarylevelset);
	// GridBasedData<Dim, double> refinedLevelSet(Refine::get(boundarylevelset.grid, 2));
	// Refine::fill(refinedLevelSet, _levelSetView);
	// const auto liquidMesh = Contourer<Dim, Dim == 2>::contour(refinedLevelSet);
	liquidMesh.writeOBJ("./boundary.obj");
}