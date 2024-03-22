#include "Advection.h"
#include "ArgsParser.h"
#include "Contourer.h"
#include "ImplicitSurface.h"
#include "Reinitialization.h"
#include "ScalarField.h"
#include "RefMapView.h"
#include "StaggeredGrid.h"
#include "Yaml.h"
#include "Refine.h"
#include "FastMarching.h"

#define _ENABLE_FMM_
#include "HybridNewton.h"

#include "../../MC-style-vol-eval/fraction.hpp"

#include <omp.h>

#include <chrono>
#include <filesystem>
#include <numbers>

using namespace PhysX;

constexpr double length = 10;

std::ofstream volFout;

inline auto BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("output", 'o', "the output directory", "output");
	parser->addArgument<int>("test", 't', "the test case index", 0);
	parser->addArgument<std::string>("mode", 'm', "the mode of advection");
	parser->addArgument<bool>("refmap", 'r', "enable the reference map", false);
	parser->addArgument<uint>("end", 'e', "the end frame (excluding)", 201);
	parser->addArgument<int>("scale", 's', "the scale of grid", 100);
	return parser;
}

inline void prepareWriter(const std::string outputDir, const double radius)
{
	std::filesystem::create_directory(outputDir);
	{ // Write description.
		YAML::Node root;
		root["Dimension"] = 3;
		root["Radius"] = float(radius);
		// { // Description of level set.
		// 	YAML::Node node;
		// 	node["Name"] = "levelSet";
		// 	node["data_mode"] = "dynamic";
		// 	node["primitive_type"] = "point_list";
		// 	node["indexed"] = false;
		// 	node["color_map"]["enabled"] = true;
		// 	root["objects"].push_back(node);
		// }
		{ // Description of contour.
			YAML::Node node;
			node["Name"] = "contour";
			node["Animated"] = true;
			node["Indexed"] = true;
			// node["primitive_type"] = "triangle_list";
			// node["material"]["diffuse_albedo"] = Vector4f(0, 0, 1, 1);
			// node["indexed"] = true;
			root["Objects"].push_back(node);
		}
		// { // Description of markers.
		// 	YAML::Node node;
		// 	node["name"] = "markers";
		// 	node["data_mode"] = "dynamic";
		// 	node["primitive_type"] = "point_list";
		// 	node["material"]["diffuse_albedo"] = Vector4f(1, 0, 1, 1);
		// 	node["indexed"] = false;
		// 	root["objects"].push_back(node);
		// }
		std::ofstream fout(outputDir + "/description.yaml");
		fout << root;
	}
}

double calculateVolume(const GridBasedData<3, double> &phi)
{
	double volume = 0.0;
	double dx3 = phi.grid.spacing * phi.grid.spacing * phi.grid.spacing;
	forEach(phi.grid, [&](const Vector3i &cell) {
		if (phi.grid.isValid(cell + Vector3i::Ones())) {
			std::array<double, 8> stencil = phi.template stencil<1>(cell);
			volume += Fraction::get_mc_vol(stencil) * dx3;
		}
	});
	return volume;
}

inline double writeFrame(const GridBasedData<3, double> &phi, const std::vector<Vector3d> &pPos, const std::string outputDir, const int frame)
{
	{ // Write the last frame.
		std::ofstream fout(outputDir + "/frame_count.txt");
		fout << frame + 1 << std::endl;
	}
	std::filesystem::create_directory(std::format("{}/frames", outputDir));
	std::filesystem::create_directory(std::format("{}/render", outputDir));
	const std::string frameDir = std::format("{}/frames/{}", outputDir, frame);
	std::filesystem::create_directory(frameDir);
	{ // Write level set.
		// std::ofstream fout(frameDir + "/levelSet.mesh", std::ios::binary);
		// IO::writeValue(fout, uint(phi.grid.count()));
		// forEach(phi.grid, [&](const Vector3i &cell) {
		// 	const Vector3d pos = phi.grid.position(cell);
		// 	IO::writeValue(fout, pos.template cast<float>().eval());
		// });
		// forEach(phi.grid, [&](const Vector3i &cell) {
		// 	IO::writeValue(fout, float(phi[cell]));
		// });
	}
	{ // Write contour.
		std::ofstream fout(frameDir + "/contour.out", std::ios::binary);
		auto mesh = Contourer<3>::contour(phi);
		mesh.write(fout);
		mesh.writeOBJ(std::format("{}/render/{}.obj", outputDir, frame));
		// for (size_t i = 0; i < mesh.indices.size(); i += 3) {
		// 	const Vector3d a = mesh.positions[mesh.indices[i]];
		// 	const Vector3d b = mesh.positions[mesh.indices[i + 1]];
		// 	const Vector3d c = mesh.positions[mesh.indices[i + 2]];
		// 	const Vector3d ab = b - a, ac = c - a;
		// 	volume += (ab[0] * ac[1] - ab[1] * ac[0]) / 2;
		// }
	}
	// { // Write markers.
	// 	std::ofstream fout(frameDir + "/markers.mesh", std::ios::binary);
	// 	IO::writeValue(fout, uint(pPos.size()));
	// 	for (const auto &pos : pPos) {
	// 		IO::writeValue(fout, pos.template cast<float>().eval());
	// 	}
	// }
	double volume = calculateVolume(phi);
	volFout << std::format("{:.6e}", volume) << std::endl;
	return volume;
}

class TestField : public VectorField<3>
{
public:
	double coeff = 1.0;
	TestField() = default;
	TestField(double c) : coeff(c) { }
	virtual Vector3d operator()(const Vector3d &pos) const override = 0;
	virtual Matrix3d gradient(const Vector3d &pos) const override = 0;
	inline void setCoeff(double c) { coeff = c; }
};

class Enright : public TestField
{
public:
	using TestField::coeff;
	Enright() = default;
	Enright(double c) : TestField(c) { }

	virtual Vector3d operator()(const Vector3d &pos) const override
	{
		double x = 0.5 + pos(0) / length;
		double y = 0.5 + pos(1) / length;
		double z = 0.5 + pos(2) / length;
		return Vector3d(
			2 * std::sin(x * std::numbers::pi) * std::sin(x * std::numbers::pi) * std::sin(2 * y * std::numbers::pi) * std::sin(2 * z * std::numbers::pi),
			-std::sin(y * std::numbers::pi) * std::sin(y * std::numbers::pi) * std::sin(2 * x * std::numbers::pi) * std::sin(2 * z * std::numbers::pi),
			-std::sin(z * std::numbers::pi) * std::sin(z * std::numbers::pi) * std::sin(2 * x * std::numbers::pi) * std::sin(2 * y * std::numbers::pi)
		) * 2.5 * coeff;
	}

	virtual Matrix3d gradient(const Vector3d &pos) const override
	{
		double x = 0.5 + pos(0) / length;
		double y = 0.5 + pos(1) / length;
		double z = 0.5 + pos(2) / length;
		Matrix3d ret;
		ret << 2 * std::sin(2 * x * std::numbers::pi) * std::sin(2 * y * std::numbers::pi) * std::sin(2 * z * std::numbers::pi),
			4 * std::sin(x * std::numbers::pi) * std::sin(x * std::numbers::pi) * std::cos(2 * y * std::numbers::pi) * std::sin(2 * z * std::numbers::pi),
			4 * std::sin(x * std::numbers::pi) * std::sin(x * std::numbers::pi) * std::sin(2 * y * std::numbers::pi) * std::cos(2 * z * std::numbers::pi),

			-2 * std::sin(y * std::numbers::pi) * std::sin(y * std::numbers::pi) * std::cos(2 * x * std::numbers::pi) * std::sin(2 * z * std::numbers::pi),
			-std::sin(2 * y * std::numbers::pi) * std::sin(2 * x * std::numbers::pi) * std::sin(2 * z * std::numbers::pi),
			-2 * std::sin(y * std::numbers::pi) * std::sin(y * std::numbers::pi) * std::sin(2 * x * std::numbers::pi) * std::cos(2 * z * std::numbers::pi),

			-2 * std::sin(z * std::numbers::pi) * std::sin(z * std::numbers::pi) * std::cos(2 * x * std::numbers::pi) * std::sin(2 * y * std::numbers::pi),
			-2 * std::sin(z * std::numbers::pi) * std::sin(z * std::numbers::pi) * std::sin(2 * x * std::numbers::pi) * std::cos(2 * y * std::numbers::pi),
			-std::sin(2 * z * std::numbers::pi) * std::sin(2 * x * std::numbers::pi) * std::sin(2 * y * std::numbers::pi);
		return ret * std::numbers::pi / length * 2.5 * coeff;
	}
};

inline std::unique_ptr<Enright> buildCase0(GridBasedData<3, double> &phi, GridBasedData<3, Vector3d> &psi, GridBasedData<3, double> &rphi)
{
	const ImplicitSphere<3> sphere(-Vector3d(1., 1., 1.) * length * .15, .15 * length);
	parallelForEach(phi.grid, [&](const Vector3i &cell) {
		const Vector3d pos = phi.grid.position(cell);
		phi[cell] = sphere.signedDistance(pos);
		psi[cell] = sphere.closestNormal(pos);
	});
	parallelForEach(rphi.grid, [&](const Vector3i &cell) {
		const Vector3d pos = rphi.grid.position(cell);
		rphi[cell] = sphere.signedDistance(pos);
	});
	return std::make_unique<Enright>();
}

class Twist3D : public TestField
{
public:
	virtual Vector3d operator()(const Vector3d &pos) const override { return 4 / (pos.template segment<2>(0).norm() + .1) * Vector3d(pos[1], -pos[0], 0); }

	virtual Matrix3d gradient(const Vector3d &pos) const override
	{
		const double r = pos.template segment<2>(0).norm();
		Matrix3d ret;
		ret << -pos[0] * pos[1], pos[0] * pos[0] + r * .1, 0, -pos[1] * pos[1] - r * .1, pos[0] * pos[1], 0, 0, 0, 0;
		return ret * 4 / r / (r + .1) / (r + .1);
	}
};

inline std::unique_ptr<Twist3D> buildCase1(GridBasedData<3, double> &phi, GridBasedData<3, Vector3d> &psi, GridBasedData<3, double> &rphi)
{
	const ImplicitSphere<3> sphere(Vector3d::Unit(1) * length * .25, .2 * length);
	parallelForEach(phi.grid, [&](const Vector3i &cell) {
		const Vector3d pos = phi.grid.position(cell);
		phi[cell] = sphere.signedDistance(pos);
		psi[cell] = sphere.closestNormal(pos);
	});
	parallelForEach(rphi.grid, [&](const Vector3i &cell) {
		const Vector3d pos = rphi.grid.position(cell);
		rphi[cell] = sphere.signedDistance(pos);
	});
	return std::make_unique<Twist3D>();
}

void advectParticle(std::vector<Vector3d> &pPos, const VectorField<3> &flow, const double dt)
{
	for (auto &pos : pPos) {
		const Vector3d vel0 = flow(pos);
		const Vector3d vel1 = flow(pos + vel0 * dt);
		pos += (vel0 + vel1) * dt / 2;
	}
}

void correctLevelSet(GridBasedData<3, double> &phi, const std::vector<Vector3d> &pPos, const std::vector<double> &pRad, const std::vector<int> &pSign)
{
	GridBasedData<3, double> phiPositive = phi;
	GridBasedData<3, double> phiNegative = phi;
	for (int i = 0; i < pPos.size(); i++) {
		const Vector3d pos = pPos[i];
		const double rad = pRad[i];
		const int sign = pSign[i];
		const double tPhi = ScalarFieldView<LinearIntrpl<3>>(phi)(pos);
		if (tPhi * sign < 0 && std::abs(tPhi) > rad) {
			const Vector3i lower = phi.grid.getLinearLower(pos);
			for (int j = 0; j <= 1; j++)
				for (int k = 0; k <= 1; k++)
					for (int l = 0; l <= 1; l++) {
						const Vector3i coord = lower + Vector3i(j, k, l);
						if (sign > 0)
							phiPositive[coord] = std::max(phiPositive[coord], rad - (pos - phi.grid.position(coord)).norm());
						else
							phiNegative[coord] = std::min(phiNegative[coord], (pos - phi.grid.position(coord)).norm() - rad);
					}
		}
	}
	parallelForEach(phi.grid, [&](const Vector3i &coord) {
		if (std::abs(phiPositive[coord]) <= std::abs(phiNegative[coord])) phi[coord] = phiPositive[coord];
		else phi[coord] = phiNegative[coord];
	});
}

void resample(const GridBasedData<3, double> &phi, std::vector<Vector3d> &pPos, std::vector<double> &pRad, std::vector<int> &pSign)
{
	pPos.clear();
	pRad.clear();
	pSign.clear();
	const double dx = phi.grid.spacing;
	forEach(phi.grid, [&](const Vector3i &cell) {
		if (std::abs(phi[cell]) <= 3 * dx) {
			for (int i = 0; i < 256; i++) {
				const Vector3d pos = phi.grid.position(cell) + Vector3d::Random() * .5 * dx;
				const double phiVal = ScalarFieldView<LinearIntrpl<3>>(phi)(pos);
				const int sign = phiVal < 0 ? -1 : 1;
				const double rad = std::clamp(sign * phiVal, .1 * dx, .5 * dx);
				pPos.push_back(pos);
				pRad.push_back(rad);
				pSign.push_back(sign);
			}
		}
	});
}

void reradius(const GridBasedData<3, double> &phi, std::vector<Vector3d> &pPos, std::vector<double> &pRad, std::vector<int> &pSign)
{
	const double dx = phi.grid.spacing;
	for (int i = 0; i < pPos.size(); i++) {
		const Vector3d pos = pPos[i];
		double &rad = pRad[i];
		const int sign = pSign[i];
		rad = std::clamp(sign * ScalarFieldView<LinearIntrpl<3>>(phi)(pos), .1 * dx, .5 * dx);
	}
}

int main(int argc, char *argv[])
{
#ifdef _OPENMP
	omp_set_num_threads(std::max(omp_get_num_procs() / 3, 1));
#endif

	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

	const auto output = std::any_cast<std::string>(parser->getValueByName("output"));
	const auto test = std::any_cast<int>(parser->getValueByName("test"));
	const auto mode = std::any_cast<std::string>(parser->getValueByName("mode"));
	const auto refmap = std::any_cast<bool>(parser->getValueByName("refmap"));
	const auto end = std::any_cast<uint>(parser->getValueByName("end"));
	const auto scale = std::any_cast<int>(parser->getValueByName("scale"));

	if (mode != "linear" && mode != "quadratic" && mode != "cubic" && mode != "bfecc" && mode != "hermite" && mode != "hnewton" && mode != "hfmm" && mode != "particle") {
		std::cerr << "Error: [main] encountered invalid mode." << std::endl;
		std::exit(-1);
	}

	const double dt = .02;
	const Vector3i resolution = scale * Vector3i::Ones();
	StaggeredGrid<3> sGrid(0, length / scale, resolution);
	GridBasedData<3, Vector3d> points(sGrid.cellGrid);
	GridBasedData<3, double> phi(sGrid.cellGrid);
	GridBasedData<3, Vector3d> psi(sGrid.cellGrid);
	GridBasedData<3, double> curPhi(sGrid.cellGrid);
	GridBasedData<3, Vector3d> curPsi(sGrid.cellGrid);
	GridBasedData<3, Vector3d> xi(sGrid.cellGrid);
	GridBasedData<3, Matrix3d> pi(sGrid.cellGrid);
	std::unique_ptr<TestField> velocityPtr;

	Grid<3> gridEx(Refine::get(sGrid.cellGrid, 8));
	GridBasedData<3, double> refinePhi(gridEx);

	switch (test) {
	case 0:
		velocityPtr = buildCase0(phi, psi, refinePhi);
		break;
	case 1:
		velocityPtr = buildCase1(phi, psi, refinePhi);
		break;
	default:
		std::cerr << "Error: [main] encountered invalid test case." << std::endl;
		std::exit(-1);
	}

	std::vector<Vector3d> pPos;
	std::vector<double> pRad;
	std::vector<int> pSign;
	if (mode == "particle") {
		// Reinitialization::solve<FastMarching<3>>(phi, -1);
		resample(phi, pPos, pRad, pSign);
	}

	prepareWriter(output, sGrid.radius());
	volFout.open(output + "/volume.txt");
	// std::cout << std::format("Initial volume: {:.6e}", writeFrame(phi, pPos, output, 0)) << std::endl;
	std::cout << std::format("Initial volume: {:.6e}", writeFrame(refinePhi, pPos, output, 0)) << std::endl;

	parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
		curPhi[cell] = phi[cell];
		curPsi[cell] = psi[cell];
	});
	auto phiViewLinear = ScalarFieldView<LinearIntrpl<3>>(phi);
	auto phiViewQuad = ScalarFieldView<EnhancedLinearIntrpl<3, double>>(phi);
	auto phiViewCubic = ScalarFieldView<CubicIntrpl<3>>(phi);
	auto phiViewBfecc = ScalarFieldView<QuadraticIntrpl<3>>(phi);
	auto phiViewHermite = ScalarFieldView<HermiteIntrpl<3, double>>(phi, psi);
	// auto phiView1Linear = ScalarFieldView<LinearIntrpl<3>>(curPhi);
	// auto phiView1Quad = ScalarFieldView<EnhancedLinearIntrpl<3, double>>(curPhi);
	// auto phiView1Cubic = ScalarFieldView<CubicIntrpl<3>>(curPhi);
	// auto phiView1Bfecc = ScalarFieldView<QuadraticIntrpl<3>>(curPhi);
	// auto phiView1Hermite = ScalarFieldView<HermiteIntrpl<3, double>>(curPhi, curPsi);
	auto phiView1Linear = RefMapFieldView<LinearIntrpl<3>>(phiViewLinear, xi);
	auto phiView1Quad = RefMapFieldView<EnhancedLinearIntrpl<3, Vector3d>>(phiViewQuad, xi);
	auto phiView1Cubic = RefMapFieldView<CubicIntrpl<3>>(phiViewCubic, xi);
	auto phiView1Bfecc = RefMapFieldView<QuadraticIntrpl<3>>(phiViewBfecc, xi);
	auto phiView1Hermite = RefMapFieldView<HermiteIntrpl<3, Vector3d>>(phiViewHermite, xi, pi);

	if (refmap) { // Initialize reference map.
		parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
			xi[cell] = sGrid.cellCenter(cell);
			pi[cell] = Matrix3d::Identity();
		});
	}

	const auto beginTime = std::chrono::steady_clock::now();
	constexpr int RkOrder = 3;

	for (uint i = 1; i < end; i++) {
		velocityPtr->setCoeff(std::cos(std::numbers::pi * i / end));
		std::cout << std::format("[{:>4}] Advect by {} scheme, {} reference map... ", i, mode, refmap ? "with" : "without") << std::flush;

		if (refmap) {
			if (mode == "linear") {
				Advection::solve<SemiLagrangian<3, RkOrder, LinearIntrpl<3>>>(xi, *velocityPtr, dt);
				parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
					curPhi[cell] = phiViewLinear(xi[cell]);
				});
				phiView1Linear.reset();
				Refine::fill(refinePhi, phiView1Linear);
			}
			else if (mode == "quadratic") {
				Advection::solve<SemiLagrangian<3, RkOrder, EnhancedLinearIntrpl<3, Vector3d>>>(xi, *velocityPtr, dt);
				parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
					curPhi[cell] = phiViewQuad(xi[cell]);
				});
				phiView1Quad.reset();
				Refine::fill(refinePhi, phiView1Quad);
			}
			else if (mode == "cubic") {
				Advection::solve<SemiLagrangian<3, RkOrder, CubicIntrpl<3>>>(xi, *velocityPtr, dt);
				parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
					curPhi[cell] = phiViewCubic(xi[cell]);
				});
				phiView1Cubic.reset();
				Refine::fill(refinePhi, phiView1Cubic);
			}
			else if (mode == "bfecc") {
				Advection::solve<MacCormack<3, RkOrder>>(xi, *velocityPtr, dt);
				parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
					curPhi[cell] = phiViewBfecc(xi[cell]);
				});
				phiView1Bfecc.reset();
				Refine::fill(refinePhi, phiView1Bfecc);
			}
			else {
				Advection::solve<SemiLagrangian<3, RkOrder, HermiteIntrpl<3, Vector3d>>>(xi, pi, *velocityPtr, dt);
				parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
					curPhi[cell] = phiViewHermite(xi[cell]);
					curPsi[cell] = pi[cell].transpose() * phiViewHermite.gradient(xi[cell]);
				});
				phiView1Hermite.reset();
				if (mode == "hnewton" && i % 10 == 1) {
					parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
						phi[cell] = curPhi[cell];
					});
					Reinitialization::solve<FastMarching<3>>(phi, -1);
					phiViewLinear.reset();
					parallelForEach(points.grid, [&](const Vector3i &cell) {
						Vector3d pos = points.grid.position(cell);
						points[cell] = pos - phiViewLinear(pos) * phiViewLinear.gradient(pos);
					});
					HybridNewton<3> reinizer(sGrid.cellGrid);
					reinizer.perform(curPhi, curPsi, points, 20, 0.9);
					parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
						phi[cell] = curPhi[cell];
						psi[cell] = curPsi[cell];
						xi[cell] = sGrid.cellCenter(cell);
						pi[cell] = Matrix3d::Identity();
					});
					phiViewHermite.reset();
					phiView1Hermite.reset();
				}
				if (mode == "hfmm" && i % 10 == 1) {
					Reinitialization::solve<5>(curPhi, 50);
					Derivatives::computeFirst<2>(curPhi, curPsi);
					parallelForEach(sGrid.cellGrid, [&](const Vector3i &cell) {
						phi[cell] = curPhi[cell];
						psi[cell] = curPsi[cell];
						xi[cell] = sGrid.cellCenter(cell);
						pi[cell] = Matrix3d::Identity();
					});
					phiViewHermite.reset();
					phiView1Hermite.reset();
				}
				Refine::fill(refinePhi, phiView1Hermite);
			}

		}
		else {
			if (mode == "linear") {
				Advection::solve<SemiLagrangian<3, RkOrder, LinearIntrpl<3>>>(phi, *velocityPtr, dt);
				phiViewLinear.reset();
				Refine::fill(refinePhi, phiViewLinear);
			}
			else if (mode == "quadratic") {
				Advection::solve<SemiLagrangian<3, RkOrder, EnhancedLinearIntrpl<3, double>>>(phi, *velocityPtr, dt);
				phiViewQuad.reset();
				Refine::fill(refinePhi, phiViewQuad);
			}
			else if (mode == "cubic") {
				Advection::solve<SemiLagrangian<3, RkOrder, CubicIntrpl<3>>>(phi, *velocityPtr, dt);
				phiViewCubic.reset();
				Refine::fill(refinePhi, phiViewCubic);
			}
			else if (mode == "bfecc") {
				Advection::solve<MacCormack<3, RkOrder>>(phi, *velocityPtr, dt);
				phiViewBfecc.reset();
				Refine::fill(refinePhi, phiViewBfecc);
			}
			else if (mode == "particle") {
				Advection::solve<SemiLagrangian<3, 2, LinearIntrpl<3>>>(phi, *velocityPtr, dt);
				advectParticle(pPos, *velocityPtr, dt);
				correctLevelSet(phi, pPos, pRad, pSign);
				Reinitialization::solve<FastMarching<3>>(phi, -1);
				correctLevelSet(phi, pPos, pRad, pSign);
				reradius(phi, pPos, pRad, pSign);

				phiViewLinear.reset();
				Refine::fill(refinePhi, phiViewLinear);
			}
			else {
				Advection::solve<SemiLagrangian<3, RkOrder, HermiteIntrpl<3, double>>>(phi, psi, *velocityPtr, dt);
				phiViewHermite.reset();
				if (mode == "hnewton") {
					HybridNewton<3> reinizer(sGrid.cellGrid);
					reinizer.perform(phi, psi, 20, 0.9);
					phiViewHermite.reset();
				}
				if (mode == "hfmm") {
					Reinitialization::solve<FastMarching<3>>(phi, 50);
					Derivatives::computeFirst<2>(phi, psi);
					phiViewHermite.reset();
				}
				Refine::fill(refinePhi, phiViewHermite);
			}
		}

		// std::cout << std::format("Volume: {:.6e}", writeFrame(refmap ? curPhi : phi, pPos, output, i)) << std::endl;
		std::cout << std::format("Volume: {:.6e}", writeFrame(refinePhi, pPos, output, i)) << std::endl;
	}

	const auto endTime = std::chrono::steady_clock::now();
	std::cout << std::format("Time: {:>10.3f}s used", std::chrono::duration<double>(endTime - beginTime).count()) << std::endl;

	return 0;
}
