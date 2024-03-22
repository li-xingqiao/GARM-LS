#include "LevelSetLiquid.h"

#include "Advection.h"
#include "Constants.h"
#include "Contourer.h"
#include "Extrapolation.h"
#include "Reinitialization.h"
#include "LambdaCorr.h"
#include "HybridNewton.h"
#include "MCReinitializer.h"
#include "SurfaceTension.h"
#include "Refine.h"
#include "SdfToVDB.h"
#include "RefMapView.h"
#include "FMMExtrapolator.h"
#include "LevelSetSurface.h"
#include "TransformedSurface.h"

#include "../../MC-style-vol-eval/fraction.hpp"

#include <filesystem>
#include <limits>
#include <chrono>

namespace PhysX {

static std::mt19937 engine(static_cast<std::mt19937::result_type>(19260817));
static std::uniform_real_distribution<double> distribution(0.0, 1.0);
static inline double getRandomNumber(double num1, double num2) {
	double ratio = distribution(engine);
	return ratio * num1 + (1 - ratio) * num2;
}
template <int Dim>
static inline Vector<Dim, double> getRandomPosition(const Vector<Dim, double> &pos1, const Vector<Dim, double> &pos2)
{
	Vector<Dim, double> result;
	for (int i = 0; i < Dim; ++i) {
		double ratio = distribution(engine);
		result[i] = ratio * pos1[i] + (1 - ratio) * pos2[i];
	}
	return result;
}

// #define CORRECTEDADVECTION

template <int Dim>
LevelSetLiquid<Dim>::LevelSetLiquid(const StaggeredGrid<Dim> &sGrid) :
	_sGrid(sGrid),
	_velocity(_sGrid.faceGrids),
	_levelSet(_sGrid.cellGrid, std::numeric_limits<double>::infinity()),
#ifdef _ENABLE_NEWTON_
	_levelSetGrad(_sGrid.cellGrid),
#endif
	_levelSetView(_levelSet),
#ifndef _DISABLE_RM_
	_refMap(_sGrid.cellGrid),
	_refMapGrad(_sGrid.cellGrid),
	_tmpLevelSet(_sGrid.cellGrid),
	_advRegion(_sGrid.cellGrid),
	_oriLevelSet(_sGrid.cellGrid),
	_refLevelSet(_sGrid.cellGrid),
#ifndef _DISABLE_GA_
	_refLevelSetGrad(_sGrid.cellGrid),
	_refLevelSetView(_refLevelSet, _refLevelSetGrad),
#else
	_refLevelSetView(_refLevelSet),
#endif
#endif
	_velocityDiv(_sGrid.cellGrid),
	_boundary(_sGrid),
	_projector(_sGrid),
	_reinitmark(_sGrid.cellGrid),
	_refColor(_sGrid.cellGrid),
	_refColorView(_refColor),
	_color(_sGrid.cellGrid),
	_initialBoundary(_sGrid.cellGrid)
{ }

template <int Dim>
void LevelSetLiquid<Dim>::writeDescription(YAML::Node &root) const
{
	root["radius"] = float(_sGrid.radius());
	{ // Description of neumann
		YAML::Node node;
		node["name"] = "neumann";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["indexed"] = false;
		node["material"]["diffuse_albedo"] = Vector4f(0.5f, 0.5f, 0.5, 1.0f);
		root["objects"].push_back(node);
	}
	// { // Description of drops
	// 	YAML::Node node;
	// 	node["name"] = "drops";
	// 	node["data_mode"] = "dynamic";
	// 	node["primitive_type"] = "point_list";
	// 	node["indexed"] = false;
	// 	node["material"]["diffuse_albedo"] = Vector4f(1.0f, 0.0f, 0.0f, 1.0f);
	// 	root["objects"].push_back(node);
	// }
	// if constexpr (Dim == 2) { // Description of velocity.
	// 	YAML::Node node;
	// 	node["name"] = "velocity";
	// 	node["data_mode"] = "dynamic";
	// 	node["primitive_type"] = "line_list";
	// 	node["indexed"] = false;
	// 	node["color_map"]["enabled"] = true;
	// 	root["objects"].push_back(node);
	// }
	// if constexpr (Dim == 2) { // Description of velocity at its point.
	// 	YAML::Node node;
	// 	node["name"] = "velocityact";
	// 	node["data_mode"] = "dynamic";
	// 	node["primitive_type"] = "line_list";
	// 	node["indexed"] = false;
	// 	node["color_map"]["enabled"] = true;
	// 	root["objects"].push_back(node);
	// }
#ifdef _ENABLE_NEWTON_
	// if constexpr (Dim == 2) { // Description of levelSetGrad.
	// 	YAML::Node node;
	// 	node["name"] = "levelSetGrad";
	// 	node["data_mode"] = "dynamic";
	// 	node["primitive_type"] = "line_list";
	// 	node["indexed"] = false;
	// 	node["color_map"]["enabled"] = true;
	// 	root["objects"].push_back(node);
	// }
#endif
	{ // Description of liquid.
		YAML::Node node;
		node["name"] = "liquid";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "triangle_list";
		node["material"]["diffuse_albedo"] = Vector4f(0, 0, 1, 1);
		node["indexed"] = true;
		root["objects"].push_back(node);
	}
#ifdef _ENABLE_PARTICLE_
	// {
	// 	YAML::Node temp;
	// 	temp["name"] = "particle";
	// 	temp["data_mode"] = "dynamic";
	// 	temp["primitive_type"] = "point_list";
	// 	temp["indexed"] = false;
	// 	temp["material"]["diffuse_albedo"] = Vector4f(0.0f, 1.0f, 1.0f, 1.0f);
	// 	root["objects"].push_back(temp);
	// }
#endif
	if constexpr (Dim == 2) { // Description of level set.
#ifndef _DISABLE_RM_
		// {
		// 	YAML::Node temp;
		// 	temp["name"] = "rmdelta";
		// 	temp["data_mode"] = "dynamic";
		// 	temp["primitive_type"] = "point_list";
		// 	temp["indexed"] = false;
		// 	temp["color_map"]["enabled"] = true;
		// 	root["objects"].push_back(temp);
		// }
#endif
#ifdef _ENABLE_NEWTON_
		// {
		// 	YAML::Node temp;
		// 	temp["name"] = "mark";
		// 	temp["data_mode"] = "dynamic";
		// 	temp["primitive_type"] = "point_list";
		// 	temp["indexed"] = false;
		// 	temp["color_map"]["enabled"] = true;
		// 	root["objects"].push_back(temp);
		// }
#endif
		// {
		// 	YAML::Node node;
		// 	node["name"] = "levelSet";
		// 	node["data_mode"] = "dynamic";
		// 	node["primitive_type"] = "point_list";
		// 	node["indexed"] = false;
		// 	node["color_map"]["enabled"] = true;
		// 	root["objects"].push_back(node);
		// }
		// {
		// 	YAML::Node node;
		// 	node["name"] = "velocityDiv";
		// 	node["data_mode"] = "dynamic";
		// 	node["primitive_type"] = "point_list";
		// 	node["indexed"] = false;
		// 	node["color_map"]["enabled"] = true;
		// 	root["objects"].push_back(node);
		// }
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{
		// Write volume.
		std::ofstream fout(frameDir + "/volume.csv");
		size_t steps = _volume.size();
		for (int i = 0; i < steps; ++i) {
			fout << std::format("{},{:.6f}\n", _restarted.at(i), _volume.at(i));
		}
		fout.close();
	}
	{
		// Write profile.
		std::ofstream fout(frameDir + "/profile.csv");
		size_t steps = _profile.size();
		for (int i = 0; i < steps; ++i) {
			fout << std::format("{:.6f},{:.6f},{:.6f},{:.6f}\n",
			 	_profile.at(i).tadvect, _profile.at(i).trestart, _profile.at(i).tproject, _profile.at(i).textrap);
		}
		fout.close();
	}
	{ // Write neumann.
		std::ofstream fout(frameDir + "/neumann.mesh", std::ios::binary);
		uint cnt = 0;
		forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
			const VectorDd pos = _sGrid.cellCenter(cell);
			if (_boundary.domainBox.isInside(pos) && _boundary.cellDist[cell] <= 0)
				cnt++;
		});
		IO::writeValue(fout, cnt);
		forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
			const VectorDd pos = _sGrid.cellCenter(cell);
			if (_boundary.domainBox.isInside(pos) && _boundary.cellDist[cell] <= 0)
				IO::writeValue(fout, pos.template cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
				const VectorDd pos = _sGrid.cellCenter(cell);
				if (_boundary.domainBox.isInside(pos) && _boundary.cellDist[cell] <= 0)
					IO::writeValue(fout, VectorFieldView<LinearIntrpl<Dim>>(_boundary.normal)(pos).normalized().template cast<float>().eval());
			});
		}
	}
	// { // Write drops.
	// 	std::ofstream fout(frameDir + "/drops.mesh", std::ios::binary);
	// 	uint cnt = 0;
	// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 		if (checkDrop(cell)) cnt++;
	// 	});
	// 	IO::writeValue(fout, cnt);
	// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 		const VectorDd pos = _sGrid.cellCenter(cell);
	// 		if (checkDrop(cell))
	// 			IO::writeValue(fout, pos.template cast<float>().eval());
	// 	});
	// 	if constexpr (Dim == 3) {
	// 		forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 			if (checkDrop(cell))
	// 				IO::writeValue(fout, Vector3d::Unit(1).eval());
	// 		});
	// 	}
	// }
	// if constexpr (Dim == 2) { // Write velocity.
	// 	std::ofstream fout(frameDir + "/velocity.mesh", std::ios::binary);
	// 	IO::writeValue(fout, uint(2 * _sGrid.cellCount()));
	// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 		const VectorDd pos = _sGrid.cellCenter(cell);
	// 		const VectorDd dir = VectorFieldView<LinearIntrpl<Dim>>(_velocity)(pos).normalized() * _sGrid.spacing * std::sqrt(Dim) / 2;
	// 		IO::writeValue(fout, pos.template cast<float>().eval());
	// 		IO::writeValue(fout, (pos + dir).template cast<float>().eval());
	// 	});
	// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 		const VectorDd pos = _sGrid.cellCenter(cell);
	// 		const float vel = float(VectorFieldView<LinearIntrpl<Dim>>(_velocity)(pos).norm());
	// 		IO::writeValue(fout, vel);
	// 		IO::writeValue(fout, vel);
	// 	});
	// }
	// if constexpr (Dim == 2) { // Write velocity at its point.
	// 	double maxvel = _velocity.normMax();
	// 	std::ofstream fout(frameDir + "/velocityact.mesh", std::ios::binary);
	// 	IO::writeValue(fout, uint(2 * (_sGrid.faceCount(0) + _sGrid.faceCount(1))));
	// 	forEach(_sGrid.faceGrids, [&](const int axis, const VectorDi &face) {
	// 		const VectorDd pos = _sGrid.faceGrids[axis].position(face);
	// 		const VectorDd dir = VectorDd::Unit(axis) * _velocity[axis][face] / maxvel * _sGrid.spacing * 2;
	// 		IO::writeValue(fout, pos.template cast<float>().eval());
	// 		IO::writeValue(fout, (pos + dir).template cast<float>().eval());
	// 	});
	// 	forEach(_sGrid.faceGrids, [&](const int axis, const VectorDi &face) {
	// 		// const float vel = float(_velocity[axis][face]);
	// 		IO::writeValue(fout, 0.0f);
	// 		IO::writeValue(fout, 0.0f);
	// 	});
	// }
#ifdef _ENABLE_NEWTON_
	// if constexpr (Dim == 2) { // Write level set gradient.
	// 	std::ofstream fout(frameDir + "/levelSetGrad.mesh", std::ios::binary);
	// 	IO::writeValue(fout, uint(2 * _sGrid.cellCount()));
	// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 		const VectorDd pos = _sGrid.cellCenter(cell);
	// 		const VectorDd dir = VectorFieldView<LinearIntrpl<Dim>>(_levelSetGrad)(pos).normalized() * _sGrid.spacing * std::sqrt(Dim) / 2;
	// 		IO::writeValue(fout, pos.template cast<float>().eval());
	// 		IO::writeValue(fout, (pos + dir).template cast<float>().eval());
	// 	});
	// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 		const VectorDd pos = _sGrid.cellCenter(cell);
	// 		const float vel = float(VectorFieldView<LinearIntrpl<Dim>>(_levelSetGrad)(pos).norm());
	// 		IO::writeValue(fout, vel);
	// 		IO::writeValue(fout, vel);
	// 	});
	// }
#endif
	{ // Write liquid.
		std::filesystem::path framePath(frameDir);
		std::filesystem::path baseDir = framePath.parent_path().parent_path();
		std::filesystem::path frame = framePath.filename();
		std::ofstream fout(frameDir + "/liquid.mesh", std::ios::binary);
#ifndef _DISABLE_RM_
		const GridBasedData<Dim, double> &sourceLevelSet = _tmpLevelSet;
#else
		const GridBasedData<Dim, double> &sourceLevelSet = _levelSet;
#endif
		const std::filesystem::path renderDir = baseDir / "render";
		if (!std::filesystem::exists(renderDir))
			std::filesystem::create_directories(renderDir);
		std::string objfile = (renderDir / (frame.string() + ".obj")).string();
#ifdef _WITH_OPENVDB
		std::string vdbfile = (renderDir / (frame.string() + ".vdb")).string();
#endif

		SurfaceMesh<Dim> liquidMesh;
		if (_enableRefine) {
#ifdef _ENABLE_NEWTON_
			DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>> tmpLevelSetView = ScalarFieldView<HermiteIntrpl<Dim, double>>(sourceLevelSet, _levelSetGrad);
#else
			ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, CubicIntrpl<Dim>> tmpLevelSetView(sourceLevelSet);
#endif
			Grid<Dim> gridEx = Refine::get(sourceLevelSet.grid, 4);
			GridBasedData<Dim, double> refinedLevelSet(gridEx);
			Refine::fill(refinedLevelSet, tmpLevelSetView);
			liquidMesh = Contourer<Dim, Dim == 2>::contour(refinedLevelSet);
#ifdef _WITH_OPENVDB
			if constexpr (Dim == 3) saveVDB(refinedLevelSet, vdbfile);
#endif
		}
		else {
			liquidMesh = Contourer<Dim, Dim == 2>::contour(sourceLevelSet);
#ifdef _WITH_OPENVDB
			if constexpr (Dim == 3) saveVDB(sourceLevelSet, vdbfile);
#endif
		}
		liquidMesh.write(fout);
		liquidMesh.writeOBJ(objfile);

		// write color for liquidMesh
		if (_enableColorOutput) {
			std::string colorfile = (renderDir / (frame.string() + ".txt")).string();
			std::ofstream fout2(colorfile);
			auto colorView = ScalarFieldView<LinearIntrpl<Dim>>(_color);
			for (const auto &pos : liquidMesh.positions) {
				double color_pos = colorView(pos);
				fout2 << color_pos << "\n";
			}
		}
	}
#ifdef _ENABLE_PARTICLE_
	// { // Write particle.
	// 	std::ofstream fout(frameDir + "/particle.mesh", std::ios::binary);
	// 	uint cnt = _particles.size();
	// 	IO::writeValue(fout, cnt);
	// 	for (const auto &p : _particles) {
	// 		IO::writeValue(fout, p.pos.template cast<float>().eval());
	// 	}
	// 	for (const auto &p : _particles) {
	// 		IO::writeValue(fout, Vector3d::Unit(1).eval());
	// 	}
	// }
#endif
	if constexpr (Dim == 2) { // Write level set.
#ifndef _DISABLE_RM_
		// {
		// 	std::ofstream fout(frameDir + "/rmdelta.mesh", std::ios::binary);
		// 	IO::writeValue(fout, uint(_sGrid.cellCount()));
		// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		// 		const VectorDd pos = _sGrid.cellCenter(cell);
		// 		IO::writeValue(fout, pos.template cast<float>().eval());
		// 	});
		// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		// 		IO::writeValue(fout, float((_refMap[cell] - _sGrid.cellCenter(cell)).norm()));
		// 	});
		// }
#endif
// #ifdef _ENABLE_NEWTON_
// 		{
// 			std::ofstream fout(frameDir + "/mark.mesh", std::ios::binary);
// 			IO::writeValue(fout, uint(_sGrid.cellCount()));
// 			forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
// 				const VectorDd pos = _sGrid.cellCenter(cell);
// 				IO::writeValue(fout, pos.template cast<float>().eval());
// 			});
// 			forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
// 				IO::writeValue(fout, float(_reinitmark[cell] ? 1.f : 0.f));
// 			});
// 		}
// #endif
		// {
		// 	std::ofstream fout(frameDir + "/levelSet.mesh", std::ios::binary);
		// 	IO::writeValue(fout, uint(_sGrid.cellCount()));
		// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		// 		const VectorDd pos = _sGrid.cellCenter(cell);
		// 		IO::writeValue(fout, pos.template cast<float>().eval());
		// 	});
		// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		// 		IO::writeValue(fout, float(_levelSet[cell]));
		// 	});
		// }
		// {
		// 	std::ofstream fout(frameDir + "/velocityDiv.mesh", std::ios::binary);
		// 	IO::writeValue(fout, uint(_sGrid.cellCount()));
		// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		// 		const VectorDd pos = _sGrid.cellCenter(cell);
		// 		IO::writeValue(fout, pos.template cast<float>().eval());
		// 	});
		// 	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		// 		IO::writeValue(fout, float(_velocityDiv[cell]));
		// 	});
		// }
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::saveFrame(const std::string &frameDir) const
{
	if (_enableResume) {
	{ // Save velocity.
		std::ofstream fout(frameDir + "/velocity.sav", std::ios::binary);
		_velocity.save(fout);
	}
#ifndef _DISABLE_RM_
	{ // Save reference map.
		std::ofstream fout(frameDir + "/refMap.sav", std::ios::binary);
		_refMap.save(fout);
	}
	{ // Save reference level set.
		std::ofstream fout(frameDir + "/refLevelSet.sav", std::ios::binary);
		_refLevelSet.save(fout);
	}
	{ // Save advection region.
		std::ofstream fout(frameDir + "/advRegion.sav", std::ios::binary);
		_advRegion.save(fout);
	}
#ifndef _DISABLE_GA_
	{ // Save reference map gradient.
		std::ofstream fout(frameDir + "/refMapGrad.sav", std::ios::binary);
		_refMapGrad.save(fout);
	}
#ifdef _ENABLE_NEWTON_
	{ // Save reference level set gradient.
		std::ofstream fout(frameDir + "/refLevelSetGrad.sav", std::ios::binary);
		_refLevelSetGrad.save(fout);
	}
#endif
#endif
#else
	{ // Save level set.
		std::ofstream fout(frameDir + "/levelSet.sav", std::ios::binary);
		_levelSet.save(fout);
	}
#endif
	if (_enableColorOutput) {
		std::ofstream fout(frameDir + "/color.sav", std::ios::binary);
		_color.save(fout);
		fout.close();
		fout.open(frameDir + "/refColor.sav", std::ios::binary);
		_refColor.save(fout);
	}
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::loadFrame(const std::string &frameDir)
{
	if (_enableResume) {
	_boundary.finish(_sGrid);

	{ // Load velocity.
		std::ifstream fin(frameDir + "/velocity.sav", std::ios::binary);
		_velocity.load(fin);
	}
#ifndef _DISABLE_RM_
	{ // Load reference map.
		std::ifstream fin(frameDir + "/refMap.sav", std::ios::binary);
		_refMap.load(fin);
	}
	{ // Load reference level set.
		std::ifstream fin(frameDir + "/refLevelSet.sav", std::ios::binary);
		_refLevelSet.load(fin);
	}
	{ // Load advection region.
		std::ifstream fin(frameDir + "/advRegion.sav", std::ios::binary);
		_advRegion.load(fin);
	}
#ifndef _DISABLE_GA_
	{ // Load reference map gradient.
		std::ifstream fin(frameDir + "/refMapGrad.sav", std::ios::binary);
		_refMapGrad.load(fin);
	}
#ifdef _ENABLE_NEWTON_
	{ // Load reference level set gradient.
		std::ifstream fin(frameDir + "/refLevelSetGrad.sav", std::ios::binary);
		_refLevelSetGrad.load(fin);
	}
#else
	Derivatives::computeFirst<4>(_refLevelSet, _refLevelSetGrad);
#endif
#endif
	_refLevelSetView.reset();
#else
	{ // Load level set.
		std::ifstream fin(frameDir + "/levelSet.sav", std::ios::binary);
		_levelSet.load(fin);
	}
#endif
	if (_enableColorOutput) {
		std::ifstream fin(frameDir + "/color.sav", std::ios::binary);
		_color.load(fin);
		fin.close();
		fin.open(frameDir + "/color.sav", std::ios::binary);
		_refColor.load(fin);
	}
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::initialize()
{
	transformBoundary();
	_boundary.finish(_sGrid);
	enforceBoundaryLevelSet(_levelSet);

#ifndef _DISABLE_RM_
	std::cout << "Using Reference Map\n";
#ifdef _ENABLE_NEWTON_
	std::cout << "Restart with hybrid newton\n";
#elif defined(_ENABLE_MC_)
	std::cout << "Restart with marching cubes\n";
#elif defined(_ENABLE_FMM_)
	std::cout << "Restart with fast marching\n";
#else
	std::cout << "Restart with weno\n";
#endif
#ifdef _ENABLE_NEWTON_
	Derivatives::computeFirst<4>(_levelSet, _levelSetGrad);
#endif
	reinitializeLevelSet();
	_tmpLevelSet = _levelSet;
#ifndef _DISABLE_RM_
	_refColor = _color;
	_refColorView.reset();
	_refLevelSet = _levelSet;
	_oriLevelSet = _levelSet;
#ifdef _ENABLE_NEWTON_
	_refLevelSetGrad = _levelSetGrad;
#endif
	_refLevelSetView.reset();
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		_refMap[cell] = _sGrid.cellCenter(cell);
		_advRegion[cell] = 1;
#ifndef _DISABLE_GA_
		_refMapGrad[cell] = MatrixDd::Identity();
#endif
	});
#endif
	// restart();
#elif defined(_ENABLE_PARTICLE_)
	std::cout << "Using Particle Level-Set\n";
	reinitializeLevelSet();
#else
	std::cout << "Using Plain Level-Set\n";
	reinitializeLevelSet();
#endif
	_restarted.push_back(1);
	_volume.push_back(calculateVolume());

	applyFluidSources(0.0);
#ifndef _DISABLE_RM_
	_levelSet = _oriLevelSet;
#endif
	projectVelocity(0);
	extrapolateFieldsToAir();
	_boundary.enforce(_velocity);

	if (_enableBounce) {
		parallelForEach(_sGrid.faceGrids, [&](const int &axis, const VectorDi &face) {
			_boundary.fraction[axis][face] = 0;
		});
	}

#ifdef _ENABLE_PARTICLE_
	resampleParticles();
#endif

	_reinitCount = 0;
}

// #define PROFILE(var, code) code
#define PROFILE(var, code)\
    t0 = std::chrono::high_resolution_clock::now();\
    code\
    t1 = std::chrono::high_resolution_clock::now();\
	elapsed = t1 - t0;\
	(var) += elapsed.count();\

template <int Dim>
void LevelSetLiquid<Dim>::advance(const double dt)
{
	// std::cout << "Delta t=" << dt << std::endl;

	TimeInfo tinfo;
	std::chrono::high_resolution_clock::time_point t0, t1;
  	std::chrono::duration<double, std::chrono::seconds::period> elapsed;

	PROFILE(tinfo.tadvect, advectFields(dt);)

	_elapsedTime += dt * 0.5 * std::numbers::pi;
	transformBoundary();

#ifndef _DISABLE_RM_
	PROFILE(tinfo.tadvect, generateCurrentFields();)
	PROFILE(tinfo.trestart, 
		reinitializeLevelSet();
		applyFluidSources(dt);
	)
	PROFILE(tinfo.textrap, extrapolateFieldsToCollider();)
	PROFILE(tinfo.trestart, 
	// if (true) {
	if (!checkReferenceMap(0.883022221559489, 4, 1e-5)) {
	// if (_reinitCount >= INT_MAX) {
	// if (_reinitCount >= 3) {
		_reinitCount = 0;
		applyFluidSources(dt);
		restart();
		generateCurrentFields();
		extrapolateFieldsToCollider();
		// reinitializeLevelSet();
		_restarted.push_back(1);
	}
	else {
		std::cout << "            " << std::flush;
		_rainingCount += _rainingRate * dt * std::sin(_elapsedTime) * std::sin(_elapsedTime);
		_restarted.push_back(0);
	}
	)
#else
	PROFILE(tinfo.trestart, 
		reinitializeLevelSet();
		applyFluidSources(dt);
	)
	_restarted.push_back(1);
#endif

#ifdef _ENABLE_PARTICLE_
	PROFILE(tinfo.tadvect, 
		correctLevelSet();
		reradiusParticles();
	)
#endif

	_volume.push_back(calculateVolume());

	applyBodyForces(dt);

#ifndef _DISABLE_RM_
	// Extrapolation::solve(_levelSet, 4, [&](const VectorDi &cell) {
	// 	return _boundary.cellDist[cell] > 0;
	// });
#endif

	PROFILE(tinfo.tproject,
	projectVelocity(dt);
	if (_enableSurfaceTension) {
		applySurfaceTension(dt);
		trickLittleDrop(dt);
	}
	)

	PROFILE(tinfo.textrap, extrapolateFieldsToAir();)

	if (!_enableBounce) {
		enforceBoundaryConditions();
	}

#ifndef _DISABLE_RM_
	// enforceBoundaryLevelSet(_levelSet, _sGrid.spacing);
	// _levelSet = _oriLevelSet;
	// _levelSetView.reset();
#endif

	// double maxdiv = calculateVelocityDiv();
	// std::cout << "\nmax(nabla·v) (max volume change rate): " << maxdiv << "\n";

	_profile.push_back(tinfo);
	_reinitCount++;
}

template <int Dim>
void LevelSetLiquid<Dim>::unionLevelSet(const Surface<Dim> &surface)
{
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		const VectorDd pos = _sGrid.cellCenter(cell);
		_levelSet[cell] = std::min(_levelSet[cell], surface.signedDistance(pos));
	});
}

template <int Dim>
void LevelSetLiquid<Dim>::intersectLevelSet(const Surface<Dim> &surface)
{
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		const VectorDd pos = _sGrid.cellCenter(cell);
		_levelSet[cell] = std::max(_levelSet[cell], surface.signedDistance(pos));
	});
}

template <int Dim>
void LevelSetLiquid<Dim>::advectFields(const double dt)
{
#ifndef _ENABLE_PARTICLE_

#ifndef _DISABLE_RM_
#ifdef CORRECTEDADVECTION
		auto oriRefMap = _refMap;
		auto oriRefMapGrad = _refMapGrad;
		auto newRefMap = _refMap;
		auto newRefMapGrad = _refMapGrad;
		// do not use extrapolated refmap
		parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
			if (_boundary.cellDist[cell] < 0) {
				newRefMap[cell] = _sGrid.cellGrid.position(cell);
				newRefMapGrad[cell] = MatrixDd::Identity();
			}
		});
#endif
#endif

#ifndef _DISABLE_GA_
	SGridBasedData<Dim, VectorDd> velocityGrad(_sGrid.faceGrids);
	for (int axis = 0; axis < Dim; axis++)
		Derivatives::computeFirst<2>(_velocity[axis], velocityGrad[axis]);
	auto velocityView = VectorFieldView<LinearIntrpl<Dim>>(_velocity, velocityGrad);
#ifdef CORRECTEDADVECTION
	Advection::solve<SemiLagrangian<Dim, 2, HermiteIntrpl<Dim, VectorDd>>>(newRefMap, newRefMapGrad, velocityView, dt);
	Advection::solve<SemiLagrangian<Dim, 2, HermiteIntrpl<Dim, VectorDd>>>(oriRefMap, oriRefMapGrad, velocityView, dt);
#else
	Advection::solve<SemiLagrangian<Dim, 2, HermiteIntrpl<Dim, VectorDd>>>(_refMap, _refMapGrad, velocityView, dt);
#endif

#else
	auto velocityView = VectorFieldView<LinearIntrpl<Dim>>(_velocity);

#if defined(_USE_BFECC_)
#define ADVECT(x) Advection::solve<MacCormack<Dim, 2>>(x, velocityView, dt)
#elif defined(_USE_LERP_)
#define ADVECT(x) Advection::solve<SemiLagrangian<Dim, 2, LinearIntrpl<Dim>>>(x, velocityView, dt)
#else
#define ADVECT(x) Advection::solve<SemiLagrangian<Dim, 2, CubicIntrpl<Dim>>>(x, velocityView, dt)
#endif

#ifndef _DISABLE_RM_
#ifdef CORRECTEDADVECTION
		ADVECT(oriRefMap);
		ADVECT(newRefMap);
#else
		ADVECT(_refMap);
#endif
#else
	ADVECT(_levelSet);
#endif

#undef ADVECT
#endif

	// parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 	if (_levelSet[cell] < 0 || _oriLevelSet[cell] || velocityView(_sGrid.cellGrid.position(cell)).norm() > 1e-6) {
	// 		_advRegion[cell] = 1;
	// 	} else {
	// 		_advRegion[cell] = 0;
	// 	}
	// });

#ifndef _DISABLE_RM_
#ifdef CORRECTEDADVECTION
		parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
			// if (_boundary.cellDist[cell] < 0) {
				_refMap[cell] = oriRefMap[cell];
			// }
			// else {
			// 	_refMap[cell] = newRefMap[cell];
			// }
		});
#ifndef _DISABLE_GA_
		parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
			// if (_boundary.cellDist[cell] < 0) {
				_refMapGrad[cell] = oriRefMapGrad[cell];
			// } else {
			// 	_refMapGrad[cell] = newRefMapGrad[cell];
			// }
		});
#endif
#endif
#endif

#if defined(_DISABLE_GA_) && !defined(_DISABLE_RM_)
	Derivatives::computeFirst<4>(_refMap, _refMapGrad);
#endif

#else

	Advection::solve<SemiLagrangian<Dim, 2, CubicIntrpl<Dim>>>(_levelSet, VectorFieldView<LinearIntrpl<Dim>>(_velocity), dt);
	advectParticle(dt);
	correctLevelSet();

	// reinitializeLevelSet();
	// correctLevelSet();
	// reradiusParticles();

#endif

// std::chrono::high_resolution_clock::time_point t0, t1;
// t0 = std::chrono::high_resolution_clock::now();
#ifndef _USE_LERP_
	Advection::solve<SemiLagrangian<Dim, 2, CubicIntrpl<Dim>>>(_velocity, VectorFieldView<LinearIntrpl<Dim>>(_velocity), dt);
#else
	Advection::solve<SemiLagrangian<Dim, 2, LinearIntrpl<Dim>>>(_velocity, VectorFieldView<LinearIntrpl<Dim>>(_velocity), dt);
#endif
// t1 = std::chrono::high_resolution_clock::now();
// std::cout << "tv: " << std::chrono::duration<double, std::chrono::seconds::period>(t1 - t0).count() << std::endl;

}

template <int Dim>
void LevelSetLiquid<Dim>::applyBodyForces(const double dt)
{
	if (_enableGravity) {
		parallelForEach(_sGrid.faceGrids[1], [&](const VectorDi &face) {
			_velocity[1][face] -= kGravity * dt;
		});
	}
	if (_enableRotation) {
		if constexpr (Dim == 2) {
			constexpr double OMEGA = 0.3;
			auto velocityView = VectorFieldView<LinearIntrpl<Dim>, Dim>(_velocity);
			parallelForEach(_sGrid.faceGrids, [&](const int &axis, const VectorDi &face) {
				const VectorDd pos = _sGrid.faceGrids[axis].position(face);
				const double forceComp = -OMEGA * OMEGA * pos[axis] - 2 * OMEGA * VectorDd(velocityView(pos)[1], -(velocityView(pos)[0]))[axis];
				_velocity[axis][face] += forceComp * dt;
			});
		}
	}
#ifndef _DISABLE_GA_
	if (_enableBounce) {
		const double K1 = 1. / dt;
		// const double K2 = 0.3 / dt;
		parallelForEach(_sGrid.faceGrids, [&](const int &axis, const VectorDi &face) {
			const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
			const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
			if (_oriLevelSet[cell0] > 0 && _oriLevelSet[cell1] > 0) return;

			const VectorDd pos = _sGrid.faceGrids[axis].position(face);
			const auto normalView = VectorFieldView<LinearIntrpl<Dim>, Dim>(_boundary.normal);
			const auto depthView = ScalarFieldView<LinearIntrpl<Dim>, Dim>(_boundary.cellDist);
			const auto velocityView = VectorFieldView<LinearIntrpl<Dim>, Dim>(_velocity);
			const VectorDd n = normalView(pos).normalized();
			const double depth = depthView(pos);
			const double vn = velocityView(pos).dot(n);
			if (depth < 0 && n.any()) {
				const double forceComp = -n[axis] * K1 * vn * depth / _sGrid.spacing;
				// const double forceComp = -n[axis] * K2 * depth;
				_velocity[axis][face] += forceComp * dt;
			}
		});
	}
#endif
}

template <int Dim>
void LevelSetLiquid<Dim>::applyFluidSources(const double dt)
{
	GridBasedData<Dim, double> droplets(_sGrid.cellGrid);
	if (_enableRaining) {
		parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
			droplets[cell] = 1e5;
		});

		// set parameters
		VectorDd origin = _sGrid.domainOrigin();
		VectorDd lengths = _sGrid.domainLengths();
		VectorDd maxspace = origin + lengths - VectorDd::Ones() * _sGrid.spacing * 8.0;
		VectorDd minspace = origin + VectorDd::Unit(1) * 0.8 * lengths[1] + VectorDd::Ones() * _sGrid.spacing * 8.0;
		double minradius = _sGrid.spacing * 1.5;
		double maxradius = _sGrid.spacing * 3.5;

		// generate droplet
		// _rainingCount += _rainingRate * dt;
		_rainingCount += _rainingRate * dt * std::sin(_elapsedTime) * std::sin(_elapsedTime);
		while (_rainingCount > 0.0) {
			ImplicitSphere<Dim> droplet(getRandomPosition(minspace, maxspace), getRandomNumber(minradius, maxradius));
			unionSurfaceLevelSet(droplets, droplet);
			_rainingCount -= 1.0;
		}

		_source = std::make_unique<LevelSetSurfaceMap<Dim>>(droplets);
	}
	if (_source) {
		parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
			const VectorDd pos = _sGrid.cellGrid.position(cell);
			double res = _source->signedDistance(pos);
#ifndef _DISABLE_RM_
			// if (false) {
			if (res < _levelSet[cell]) {
				_refLevelSet[cell] = res;
				_refMap[cell] = pos;
				_refMapGrad[cell] = MatrixDd::Identity();
#ifdef _ENABLE_NEWTON_
				_levelSetGrad[cell] = _source->closestNormal(pos);
#endif
#ifndef _DISABLE_GA_
				_refLevelSetGrad[cell] = _source->closestNormal(pos);
#endif
			}
			_oriLevelSet[cell] = std::min(_oriLevelSet[cell], res);
#else
			_levelSet[cell] = std::min(_levelSet[cell], res);
#endif
		});

		_levelSetView.reset();

		if (_enableRaining) {
			parallelForEach(_sGrid.faceGrids, [&](const int &axis, const VectorDi &face) {
				const VectorDd pos = _sGrid.faceGrids[axis].position(face);
				if (_source->signedDistance(pos) < 4 * _sGrid.spacing) {
					// if (axis == 1)
					// 	_velocity[axis][face] = -_rainingSpeed * 2 * std::sin(_elapsedTime) * std::sin(_elapsedTime);
					// else
						_velocity[axis][face] = 0.0;
				}
			});
		}
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::transformBoundary()
{
	if (_enableMovingBoundary) {
		LevelSetSurface<Dim> boundary(_initialBoundary);
		MovedSurface<Dim> newboundary(boundary, _boundaryTransformation(_elapsedTime));
		_boundary.reset();
		_boundary.unions(newboundary);
		_boundary.finish(_sGrid);

		parallelForEach(_sGrid.faceGrids, [&](const int axis, const VectorDi &face) {
			VectorDd pos = _sGrid.faceGrids[axis].position(face);
			if (newboundary.signedDistance(pos) < _sGrid.spacing)
				_boundary.velocity[axis][face] = _boundaryVelocity(_elapsedTime, pos)[axis];
			else
				_boundary.velocity[axis][face] = 0;
		});
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::enforceBoundaryLevelSet(GridBasedData<Dim, double> &sdf, const double tolerance)
{
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		if (sdf[cell] < -_boundary.cellDist[cell] - tolerance) {
			VectorDd pos = sdf.grid.position(cell);
			sdf[cell] = -_boundary.cellDist[cell] - tolerance;
		}
	});
}

template <int Dim>
void LevelSetLiquid<Dim>::enforceBoundaryLevelSet(GridBasedData<Dim, double> &sdf, GridBasedData<Dim, Vector<Dim, double>> &drv, const double tolerance)
{
#ifdef _ENABLE_NEWTON_
	const auto norm = VectorFieldView<LinearIntrpl<Dim>, Dim>(_boundary.normal);
#endif
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		if (sdf[cell] < -_boundary.cellDist[cell] - tolerance) {
			VectorDd pos = sdf.grid.position(cell);
			sdf[cell] = -_boundary.cellDist[cell] - tolerance;
#ifdef _ENABLE_NEWTON_
			drv[cell] = -norm(pos);
#endif
		}
	});
}

template <int Dim>
void LevelSetLiquid<Dim>::cutBoundaryLevelSet(GridBasedData<Dim, double> &sdf, GridBasedData<Dim, Vector<Dim, double>> &drv, const double tolerance)
{
#ifdef _ENABLE_NEWTON_
	const auto norm = VectorFieldView<LinearIntrpl<Dim>, Dim>(_boundary.normal);
#endif
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		if (cell[1] > 7) return;
		if (_boundary.cellDist[cell] < tolerance) {
			VectorDd pos = sdf.grid.position(cell);
			sdf[cell] = tolerance - _boundary.cellDist[cell];
#ifdef _ENABLE_NEWTON_
			drv[cell] = -norm(pos);
#endif
		}
		else if (sdf[cell] < -_boundary.cellDist[cell] + tolerance) {
			VectorDd pos = sdf.grid.position(cell);
			sdf[cell] = -_boundary.cellDist[cell] + tolerance;
#ifdef _ENABLE_NEWTON_
			drv[cell] = -norm(pos);
#endif
		}
	});
}

template <int Dim>
void LevelSetLiquid<Dim>::projectVelocity(const double dt)
{
	_projector.project(_velocity, _levelSet, _boundary.fraction, _boundary.velocity);
}

template <int Dim>
void LevelSetLiquid<Dim>::applySurfaceTension(const double dt)
{
	extrapolateFieldsToAir(true);
	enforceBoundaryConditions();

	double surfaceTensionDelta = _surfaceTensionDeltaMin + (_surfaceTensionDeltaMax - _surfaceTensionDeltaMin) * 2 / (1 + std::exp(_elapsedTime / _surfaceTensionDeltaCritTime));
	SurfaceTension::solve(_velocity, _levelSetView, _boundary.fraction, _surfaceTensionCoeff / _density, surfaceTensionDelta, dt);

	_projector.reproject(_velocity, _levelSet, _boundary.fraction, _boundary.velocity);
}

template <int Dim>
void LevelSetLiquid<Dim>::trickLittleDrop(const double dt)
{
	const double v0 = _levelSet.grid.spacing / dt;
	const double vt = v0 * .15;
	forEach(_levelSet.grid, [&](const VectorDi &cell) {
		if (checkDrop(cell)) {
			const VectorDi face0 = cell;
			const VectorDi face1 = cell + VectorDi::Unit(1);
//			std::cout << "## " << cell.transpose() << std::endl;
//			std::cout << "## " << _velocity[1][face0] << ' ' << _velocity[1][face1] << std::endl;
			const double v = (_velocity[1][face0] + _velocity[1][face1]) * .5;
			if (v < vt * 1.25 && v > -vt && _velocity[1][face0] > 0 && _velocity[1][face1] < 0) {
//				std::cout << "found!" << std::endl;
				_velocity[1][face0] -= vt + v;
				_velocity[1][face1] -= vt + v;
			}
//			std::cout << std::endl;
		}
	});
}

template <int Dim>
void LevelSetLiquid<Dim>::extrapolateFieldsToAir(bool onlyv)
{
	if (!onlyv) {
#ifndef _DISABLE_RM_
		Extrapolation::solve(_refMap, _airExtrapolationRadius, [&](const VectorDi &cell) {
			return _advRegion[cell] > 0;
		});
		parallelForEach(_sGrid.cellGrid, [this](const VectorDi &cell) {
			_advRegion[cell] = (_oriLevelSet[cell] <= 0) ? 1 : 0;
		});
		Extrapolation::solve(_advRegion, _airExtrapolationRadius, [&](const VectorDi &cell) {
			return _oriLevelSet[cell] <= 0;
		});
#endif
	}
	Extrapolation::solve(_velocity, _airExtrapolationRadius + 2, [&](const int axis, const VectorDi &face) {
		const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
		const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
#ifndef _DISABLE_RM_
		return _boundary.fraction[axis][face] < 1 && (_oriLevelSet[cell0] <= 0 || _oriLevelSet[cell1] <= 0);
#else
		return _boundary.fraction[axis][face] < 1 && (_levelSet[cell0] <= 0 || _levelSet[cell1] <= 0);
#endif
	});
}

template <int Dim>
void LevelSetLiquid<Dim>::extrapolateFieldsToCollider()
{
	if (!_enableBounce) {
#ifndef _DISABLE_RM_
		Extrapolation::solve(_oriLevelSet, _boundaryExtrapolationRadius, [&](const VectorDi &cell) {
			// return _boundary.cellDist[cell] > 0;
			// should use solved fluid region
			if (_boundary.cellDist[cell] > 0) return true;
			// for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			// 	const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
			// 	const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
			// 	if (_boundary.fraction[axis][face] < 1) {
			// 		return true;
			// 	}
			// }
			return false;
		});
		Extrapolation::solve(_levelSet, _boundaryExtrapolationRadius, [&](const VectorDi &cell) {
			// return _boundary.cellDist[cell] > 0;
			// should use solved fluid region
			if (_boundary.cellDist[cell] > 0) return true;
			// for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			// 	const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
			// 	const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
			// 	if (_boundary.fraction[axis][face] < 1) {
			// 		return true;
			// 	}
			// }
			return false;
		});

		if (_enableContactAngle) {
			reinitializeLevelSet();

			FMMExtrapolator extrapolator(_sGrid);
			extrapolator.setupMark(_oriLevelSet, _boundary.cellDist, 5);
			extrapolator.extrapolate(_oriLevelSet, _levelSet, _boundary.cellDist, _boundary.normal, _contactAngle);
			// _tmpLevelSet = _levelSet;
		}
#else
	// Extrapolation::solve(_levelSet, 4, [&](const VectorDi &cell) {
	// 	return _boundary.cellDist[cell] > 0;
	// });
#endif
// #ifdef _ENABLE_NEWTON_
// 	Extrapolation::solve(_refMapGrad, 3, 3, [&](const VectorDi &cell) {
// 		return _oriLevelSet[cell] < 2.0 * _levelSet.grid.spacing;
// 	});
// #endif
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::enforceBoundaryConditions()
{
	_boundary.enforce(_velocity);
	// #ifndef _DISABLE_RM_
	// 	_boundary.enforce(_refMap);
	// #endif
}

template <int Dim>
void LevelSetLiquid<Dim>::reinitializeLevelSet()
{
#ifdef _ENABLE_FMM_
	Reinitialization::solve<FastMarching<Dim>>(_levelSet, 15);
#elif defined (_ENABLE_CORR_)
	Reinitialization::solve<LambdaCorr<Dim, 5>>(_levelSet, 30);
#else
	Reinitialization::solve<5>(_levelSet, 30, 0.9);
#endif

	_levelSetView.reset();
}

template <int Dim>
double LevelSetLiquid<Dim>::calculateVolume() const
{
#ifdef _DISABLE_RM_
	auto const &_tmpLevelSet = _levelSet;
#endif
	double res = 0.0;
	double dx = _tmpLevelSet.grid.spacing;
	if constexpr (Dim == 2) {
		auto mesh = Contourer<2, true>::contour(_tmpLevelSet);
		for (size_t i = 0; i < mesh.indices.size(); i += 3) {
			const Vector2d a = mesh.positions[mesh.indices[i]];
			const Vector2d b = mesh.positions[mesh.indices[i + 1]];
			const Vector2d c = mesh.positions[mesh.indices[i + 2]];
			const Vector2d ab = b - a, ac = c - a;
			res += (ab[0] * ac[1] - ab[1] * ac[0]) / 2;
		}
	}
	else {
		double dx3 = dx * dx * dx;
		forEach(_tmpLevelSet.grid, [&](const Vector3i &cell) {
			if (_tmpLevelSet.grid.isValid(cell + Vector3i::Ones())) {
				std::array<double, 8> stencil = _tmpLevelSet.template stencil<1>(cell);
				res += Fraction::get_mc_vol(stencil) * dx3;
			}
		});
	}
	return res;
}

template <int Dim>
double LevelSetLiquid<Dim>::calculateVelocityDiv()
{
	double maxdiv = 0.0;
	forEach(_levelSet.grid, [&](const VectorDi &cell) {
		if (_levelSet[cell] > 0) { _velocityDiv[cell] = 0.0; return; }
		double div = 0;
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
			const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
			const int side = StaggeredGrid<Dim>::cellFaceSide(i);
			const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
			const double weight = 1 - _boundary.fraction[axis][face];
			if (weight > 0)
				div += side * weight * _velocity[axis][face];
			if (weight < 1)
				div += side * (1 - weight) * _boundary.velocity[axis][face];
		}
		div = std::abs(div);
		_velocityDiv[cell] = div;
		maxdiv = std::max(div, maxdiv);
	});
	return maxdiv;
}

#ifndef _DISABLE_RM_
template <int Dim>
void LevelSetLiquid<Dim>::generateCurrentFields()
{
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		_levelSet[cell] = _refLevelSetView(_refMap[cell]);
		if (_enableColorOutput) {
			_color[cell] = _refColorView(_refMap[cell]);
		}
#ifdef _ENABLE_NEWTON_
		_levelSetGrad[cell] = _refMapGrad[cell].transpose() * _refLevelSetView.gradient(_refMap[cell]);
#endif
#ifndef _DISABLE_RM_
		if (_advRegion[cell] == 0 && _levelSet[cell] < _airExtrapolationRadius * _sGrid.cellGrid.spacing) {
			_levelSet[cell] = _airExtrapolationRadius * _levelSet.grid.spacing;
		}
#endif
	});

	enforceBoundaryLevelSet(_levelSet, 3.0 * _sGrid.spacing);
	_tmpLevelSet = _levelSet;
	if (_enablePrettifiedOutput) {
		enforceBoundaryLevelSet(_tmpLevelSet, 0.0 * _sGrid.spacing); // prettify the output
	}
	_oriLevelSet = _levelSet;
	// enforceBoundaryLevelSet(_oriLevelSet, 0.5 * _sGrid.spacing);
	_levelSetView.reset();
}

template <int Dim>
bool LevelSetLiquid<Dim>::checkReferenceMap(const double threshold, const double bandStep, const double ratio) const
{
	const auto square = [](const double x) { return x * x; };
	int tot = 0, cnt = 0;
	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		if (std::abs(_oriLevelSet[cell]) < bandStep * _sGrid.spacing && _boundary.cellDist[cell] > 0) {
			tot++;
			const MatrixDd rmgrad = _refMapGrad[cell];
			// Eigen::JacobiSVD<MatrixDd> svd(rmgrad);
			// double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
			// if (cond > 4) {
			// 	cnt++;
			// 	return;
			// }
			// if (std::abs(_levelSet[cell] - _oriLevelSet[cell]) > 0.5 * bandStep * _sGrid.spacing) {
			// 	cnt++;
			// 	return;
			// }
			if constexpr (Dim == 2) {
				const VectorDd a = rmgrad.col(0);
				const VectorDd b = rmgrad.col(1);
				const double vab = a.dot(b) * a.dot(b) / a.squaredNorm() / b.squaredNorm();
				if (vab > threshold) cnt++;
			}
			else {
				const VectorDd a = rmgrad.col(0);
				const VectorDd b = rmgrad.col(1);
				const VectorDd c = rmgrad.col(2);
				const double vab = square(a.dot(b)) / a.squaredNorm() / b.squaredNorm();
				const double vac = square(a.dot(c)) / a.squaredNorm() / c.squaredNorm();
				const double vbc = square(b.dot(c)) / b.squaredNorm() / c.squaredNorm();
				if (vab > threshold || vac > threshold || vbc > threshold) cnt++;
			}
		}
	});
	return 1.0 * cnt <= ratio * tot;
}

template <int Dim>
void LevelSetLiquid<Dim>::restart()
{
	std::cout << "  Restart..." << std::flush;

#ifdef _ENABLE_NEWTON_
	GridBasedData<Dim, VectorDd> points(_levelSet.grid);
	parallelForEach(points.grid, [&](const VectorDi &coord) {
		VectorDd pos = points.grid.position(coord);
		// std::cout << pos.transpose() << " # " << _levelSetView(pos) << " # " << _levelSetView.gradient(pos).transpose() << std::endl;
		points[coord] = pos - _levelSetView(pos) * _levelSetView.gradient(pos);
	});
	HybridNewton<Dim> reinizer(_levelSet.grid);
	RefMapWrapper<Dim, decltype(_refLevelSetView), HermiteIntrpl<Dim, VectorDd>> refMapView(_refLevelSetView, _refMap, _refMapGrad, true);
	// reinizer.perform(_oriLevelSet, _levelSetGrad, points, refMapView, 20, 0.9);
	// reinizer.perform(_oriLevelSet, _levelSetGrad, refMapView, 20, 0.9);
	reinizer.perform(_oriLevelSet, _levelSetGrad, points, 20, 0.9);
	// reinizer.perform(_oriLevelSet, _levelSetGrad, 20, 0.9);

	// cutBoundaryLevelSet(_oriLevelSet, _levelSetGrad, 2.0 * _sGrid.spacing);

	_reinitmark = reinizer._mark;
	_refLevelSet = _oriLevelSet;
	_refLevelSetGrad = _levelSetGrad;
#elif defined(_ENABLE_MC_)
	generateCurrentFields();
	RefMapWrapper<Dim, decltype(_refLevelSetView), HermiteIntrpl<Dim, VectorDd>> refMapView(_refLevelSetView, _refMap, _refMapGrad, true);
	MCReinitializer<Dim>().perform(_levelSet, _levelSetView, 20);

	_refLevelSet = _levelSet;

#ifndef _DISABLE_GA_
	Derivatives::computeFirst<4>(_refLevelSet, _refLevelSetGrad);
#endif

#else
	// const double EPS = 0.5 * _sGrid.spacing;
	// parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 	_levelSet[cell] = _refLevelSetView(_refMap[cell]) - EPS;
	// });
	// reinitializeLevelSet();
	// parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
	// 	_levelSet[cell] += EPS;
	// });

	_refLevelSet = _levelSet;

#ifndef _DISABLE_GA_
	Derivatives::computeFirst<4>(_refLevelSet, _refLevelSetGrad);
#endif

#endif

	_refLevelSetView.reset();

	if (_enableColorOutput) {
		_refColor = _color;
		_refColorView.reset();
	}

	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		_refMap[cell] = _sGrid.cellCenter(cell);

#ifndef _DISABLE_GA_
		_refMapGrad[cell] = MatrixDd::Identity();
#endif
	});
}
#endif

#ifdef _ENABLE_PARTICLE_
template <int Dim>
void LevelSetLiquid<Dim>::advectParticle(double dt) {
	auto velocityView = VectorFieldView<LinearIntrpl<Dim>, Dim>(_velocity);
	for (auto &particle : _particles) {
		const VectorDd vel0 = velocityView(particle.pos);
		const VectorDd vel1 = velocityView(particle.pos + vel0 * dt);
		particle.pos += (vel0 + vel1) * dt / 2;
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::correctLevelSet() {
	GridBasedData<Dim, double> phiPositive = _levelSet;
	GridBasedData<Dim, double> phiNegative = _levelSet;
	for (int i = 0; i < _particles.size(); i++) {
		const Particle &p = _particles[i];
		const double tPhi = _levelSetView(p.pos);
		if (tPhi * p.sign < 0 && std::abs(tPhi) > p.rad) {
			const VectorDi lower = _sGrid.cellGrid.getLinearLower(p.pos);

			if constexpr (Dim == 2) {
				for (int j = 0; j <= 1; j++)
					for (int i = 0; i <= 1; i++) {
						const Vector2i coord = lower + Vector2i(i, j);
						if (p.sign > 0)
							phiPositive[coord] = std::max(phiPositive[coord], p.rad - (p.pos - _levelSet.grid.position(coord)).norm());
						else
							phiNegative[coord] = std::min(phiNegative[coord], (p.pos - _levelSet.grid.position(coord)).norm() - p.rad);
					}
			}
			else {
				for (int k = 0; k <= 1; k++)
					for (int j = 0; j <= 1; j++)
						for (int i = 0; i <= 1; i++) {
							const Vector3i coord = lower + Vector3i(i, j, k);
							if (p.sign > 0)
								phiPositive[coord] = std::max(phiPositive[coord], p.rad - (p.pos - _levelSet.grid.position(coord)).norm());
							else
								phiNegative[coord] = std::min(phiNegative[coord], (p.pos - _levelSet.grid.position(coord)).norm() - p.rad);
						}
			}
		}
	}
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &coord) {
		if (std::abs(phiPositive[coord]) <= std::abs(phiNegative[coord])) _levelSet[coord] = phiPositive[coord];
		else _levelSet[coord] = phiNegative[coord];
	});
	_levelSetView.reset();
}

template <int Dim>
void LevelSetLiquid<Dim>::resampleParticles() {
	constexpr int sampleRate = (Dim == 3) ? 128 : 32;

	_particles.clear();
	const double dx = _sGrid.cellGrid.spacing;
	forEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		if (std::abs(_levelSet[cell]) <= 3 * dx) {
			for (int i = 0; i < sampleRate; i++) {
				const VectorDd pos = _sGrid.cellGrid.position(cell) + VectorDd::Random() * .5 * dx;
				const double phiVal = _levelSetView(pos);
				const int sign = phiVal < 0 ? -1 : 1;
				const double rad = std::clamp(sign * phiVal, .1 * dx, .5 * dx);
				_particles.emplace_back(pos, rad, sign);
			}
		}
	});
}

template <int Dim>
void LevelSetLiquid<Dim>::reradiusParticles() {
	const double dx = _sGrid.cellGrid.spacing;
	for (int i = 0; i < _particles.size(); i++) {
		const VectorDd pos = _particles[i].pos;
		const int sign = _particles[i].sign;
		_particles[i].rad = std::clamp(sign * _levelSetView(pos), .1 * dx, .5 * dx);
	}
}

#endif

template class LevelSetLiquid<2>;
template class LevelSetLiquid<3>;

}
