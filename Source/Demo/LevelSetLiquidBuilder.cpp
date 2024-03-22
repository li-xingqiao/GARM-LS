#include "LevelSetLiquidBuilder.h"

#include "SurfaceMesh.h"
#include "MeshToSdf.h"
#include "LevelSetSurface.h"
#include "VectorField.h"
#include "TransformedSurface.h"

#include <numbers>

namespace PhysX {

template <int Dim>
std::unique_ptr<LevelSetLiquid<Dim>> LevelSetLiquidBuilder<Dim>::build(const int scale, const int option, std::string parameters)
{
	switch (option) {
	case 0:
		return buildCase0(scale, parameters);
	case 1:
		return buildCase1(scale);
	default:
		reportError("invalid option");
		return nullptr;
	}
}

template <int Dim>
std::unique_ptr<LevelSetLiquid<Dim>> LevelSetLiquidBuilder<Dim>::buildCase0(int scale, std::string parameters)
{
	if (scale < 0) scale = (Dim == 2) ? 256 : 128;
	const double length = 20;
	const int boundaryWidth = 2;
	const VectorDi resolution = scale * VectorDi::Ones() - int(scale * 0.375) * VectorDi::Unit(1);

	double st, v;
	if (parameters.length() > 0) {
		std::istringstream sin(parameters);
		std::string tmp;
		getline(sin, tmp, ',');
		st = std::stod(tmp);
		getline(sin, tmp, ',');
		v = std::stod(tmp);
		std::cout << "param: " << parameters << " s=" << st << " v=" << v << std::endl;
	} else {
		st = (Dim == 2) ? 300 : 5000;
		v = (Dim == 2) ? 0.5 : 2.0;
	}

	StaggeredGrid<Dim> sGrid(boundaryWidth, length / (scale - 2 * boundaryWidth), resolution);
	auto liquid = std::make_unique<LevelSetLiquid<Dim>>(sGrid);
	// liquid->_enableSurfaceTension = true;
	liquid->_surfaceTensionCoeff = st;
	liquid->_surfaceTensionDeltaMax = 8.0;
	liquid->_surfaceTensionDeltaMin = 8.0;
	// liquid->_surfaceTensionDeltaCritTime = 0.2;

	liquid->unionLevelSet(ImplicitPlane<Dim>(-VectorDd::Unit(1) * length * 0.15, VectorDd::Unit(1)));
	liquid->unionLevelSet(ImplicitSphere<Dim>(VectorDd::Unit(1) * (-length * 0.15 + length * 0.02 + length * 0.03), length * 0.03));

	parallelForEach(sGrid.faceGrids[1], [&](const VectorDi &face) {
		liquid->_velocity[1][face] = -v * length;
	});
	return liquid;
}

template <int Dim>
std::unique_ptr<LevelSetLiquid<Dim>> LevelSetLiquidBuilder<Dim>::buildCase1(int scale)
{
	if (scale < 0) scale = (Dim == 2) ? 256 : 64;
	const double length = 1;
	const int boundaryWidth = 2;
	const VectorDi resolution = scale * (VectorDi::Ones() + VectorDi::Unit(0)) + int(0.2 * scale) * VectorDi::Unit(1);

	StaggeredGrid<Dim> sGrid(boundaryWidth, length / (scale - 2 * boundaryWidth), resolution);
	auto liquid = std::make_unique<LevelSetLiquid<Dim>>(sGrid);
	liquid->_enableSurfaceTension = true;
	liquid->_surfaceTensionCoeff = 30;
	liquid->_surfaceTensionDeltaMax = 5.5;
	liquid->_surfaceTensionDeltaMin = 5.5;

	liquid->unionLevelSet(ImplicitBox<Dim>(sGrid.domainOrigin(), (VectorDd::Ones() - VectorDd::Unit(0) * 0.6 - VectorDd::Unit(1) * 0.6) * length));
	// liquid->unionLevelSet(ImplicitBox<Dim>(sGrid.domainOrigin() + VectorDd::Unit(0) * (sGrid.domainLengths()[0] - length / 3), (VectorDd::Ones() - VectorDd::Unit(0) * 2 / 3 - VectorDd::Unit(1) * 0.9) * length));

	liquid->_boundary.unions(ImplicitCylinder<Dim>(-VectorDd::Unit(1) * length * (0.5 + 0.2 / 2), length * 0.1, length * (0.5 + 0.2 / 2)));

	return liquid;
}

template <int Dim>
void LevelSetLiquidBuilder<Dim>::reportError(const std::string &msg)
{
	std::cerr << std::format("Error: [LevelSetLiquidBuilder] encountered {}.", msg) << std::endl;
	std::exit(-1);
}

template class LevelSetLiquidBuilder<2>;
template class LevelSetLiquidBuilder<3>;

}
