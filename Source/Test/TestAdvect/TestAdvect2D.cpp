#include "Advection.h"
#include "ArgsParser.h"
#include "Contourer.h"
#include "ImplicitSurface.h"
#include "ParticleLevelSet.h"
#include "Reinitialization.h"
#include "Refine.h"
#include "ScalarField.h"
#include "StaggeredGrid.h"
#include "Yaml.h"

#define _ENABLE_FMM_
#include "HybridNewton.h"

#include <chrono>
#include <filesystem>
#include <numbers>

using namespace PhysX;

constexpr double c_Length = 10.;
constexpr double c_Dt     = .02;

double g_Time = 0;

class TestField : public VectorField<2> {
public:
	virtual std::pair<double, Vector2d> getLevelSetInit(const Vector2d &pos) const = 0;
	virtual double getLevelSetGt(const Vector2d &pos) const = 0;
};

class Rotation : public TestField {
public:
	virtual std::pair<double, Vector2d> getLevelSetInit(const Vector2d &pos) const override {
		const double phi0 = s_Sphere.signedDistance(pos);
		const double phi1 = s_Box.signedDistance(pos);
		if (phi0 >= -phi1) {
			return { phi0, s_Sphere.closestNormal(pos) };
		} else {
			return { -phi1, -s_Box.closestNormal(pos) };
		}
	}

	virtual double getLevelSetGt(const Vector2d &pos) const override {
		const double alpha = g_Time / c_Dt / 157. * std::numbers::pi;
		const double cosa = std::cos(alpha);
		const double sina = std::sin(alpha);
		const Vector2d pos0 = { cosa * pos[0] + sina * pos[1], -sina * pos[0] + cosa * pos[1] };
		return getLevelSetInit(pos0).first;
	}

	virtual Vector2d operator()(const Vector2d &pos) const override { return Vector2d(-pos[1], pos[0]) * std::numbers::pi / c_Dt / 157; }

	virtual Matrix2d gradient(const Vector2d &pos) const override {
		Matrix2d ret;
		ret << 0, -1, 1, 0;
		return ret * std::numbers::pi / c_Dt / 157; // 314 frames per rotation
	}

private:
	inline static ImplicitSphere<2> s_Sphere = ImplicitSphere<2>(Vector2d::Unit(1) * c_Length * .25, .15 * c_Length);
	inline static ImplicitBox<2> s_Box = ImplicitBox<2>(Vector2d::Unit(0) * c_Length * -.025, (Vector2d::Ones() * .05 + Vector2d::Unit(1) * .3) * c_Length);
};

class Twist : public TestField
{
public:
	virtual std::pair<double, Vector2d> getLevelSetInit(const Vector2d &pos) const override {
		return { s_Sphere.signedDistance(pos), s_Sphere.closestNormal(pos) };
	}

	virtual double getLevelSetGt(const Vector2d &pos) const override {
		const double omega = 4 / (pos.norm() + .1);
		const double alpha = omega * g_Time;
		const double cosa = std::cos(alpha);
		const double sina = std::sin(alpha);
		const Vector2d pos0 = { cosa * pos[0] - sina * pos[1], sina * pos[0] + cosa * pos[1] };
		return getLevelSetInit(pos0).first;
	}

	virtual Vector2d operator()(const Vector2d &pos) const override { return 4 / (pos.norm() + .1) * Vector2d(pos[1], -pos[0]); }

	virtual Matrix2d gradient(const Vector2d &pos) const override {
		const double r = pos.norm();
		if (r < std::numeric_limits<double>::epsilon()) return Matrix2d::Zero();
		Matrix2d ret;
		ret << -pos[0] * pos[1], pos[0] * pos[0] + r * .1, -pos[1] * pos[1] - r * .1, pos[0] * pos[1];
		return ret * 4 / r / (r + .1) / (r + .1);
	}

private:
	inline static ImplicitSphere<2> s_Sphere = ImplicitSphere<2>(Vector2d::Unit(1) * c_Length * .25, .2 * c_Length);
};

class Bell : public TestField {
public:
	virtual std::pair<double, Vector2d> getLevelSetInit(const Vector2d &pos) const override {
		return { s_Sphere.signedDistance(pos), s_Sphere.closestNormal(pos) };
	}

	virtual double getLevelSetGt(const Vector2d &pos) const override {
		return 0;
	}

	virtual Vector2d operator()(const Vector2d &pos) const override {
		const double x = (0.5 + pos(0) / c_Length) * std::numbers::pi;
		const double y = (0.5 + pos(1) / c_Length) * std::numbers::pi;
		double a = g_Time / c_Dt / 628;
		if (a < .5 + 1e-6) a -= 1. / 628;
		return Vector2d(-std::sin(x) * std::sin(x) * std::sin(2 * y), std::sin(y) * std::sin(y) * std::sin(2 * x)) * c_Length * std::cos(a * std::numbers::pi);
	}

	virtual Matrix2d gradient(const Vector2d &pos) const override {
		const double x = (0.5 + pos(0) / c_Length) * std::numbers::pi;
		const double y = (0.5 + pos(1) / c_Length) * std::numbers::pi;
		double a = g_Time / c_Dt / 628;
		if (a < .5 + 1e-6) a -= 1. / 628;
		Matrix2d ret;
		ret << -std::sin(2 * x) * std::sin(2 * y), -2 * std::sin(x) * std::sin(x) * std::cos(2 * y), 2 * std::sin(y) * std::sin(y) * std::cos(2 * x), std::sin(2 * y) * std::sin(2 * x);
		return ret * std::numbers::pi * std::cos(a * std::numbers::pi);
	}

private:
	inline static ImplicitSphere<2> s_Sphere = ImplicitSphere<2>(Vector2d::Unit(1) * c_Length * .25, .15 * c_Length);
};

std::string   g_OutputDir;
int           g_TestCaseIdx;
std::string   g_Mode;
bool          g_RefMap;
uint          g_EndFrame;
int           g_Scale;
int           g_NPPC;

std::ofstream g_VolFout;

inline double Linear() { return g_Mode == "linear"; }
inline double Quadratic() { return g_Mode == "quadratic"; }
inline double Cubic() { return g_Mode == "cubic"; }
inline double BFECC() { return g_Mode == "bfecc"; }
inline double Hermite() { return g_Mode == "hermite" || g_Mode == "hnewton" || g_Mode == "hfmm"; }
inline double Particle() { return g_Mode == "particle"; }

// inline size_t GetRefineRes() {
// 	if (Linear() || Particle()) {
// 		return 1;
// 	} else if (Quadratic() || BFECC()) {
// 		return 2;
// 	} else {
// 		return 4;
// 	}
// }

inline void PrepareWritter(const float radius) {
	std::filesystem::create_directory(g_OutputDir);
	{ // Write description.
		YAML::Node root;
		root["dimension"] = 2;
		root["radius"] = radius;
		{ // Description of contour.
			YAML::Node node;
			node["name"] = "contour";
			node["data_mode"] = "dynamic";
			node["primitive_type"] = "triangle_list";
			node["material"]["diffuse_albedo"] = Vector4f(0, 0, 1, 1);
			node["indexed"] = true;
			root["objects"].push_back(node);
		}
		std::ofstream fout(g_OutputDir + "/description.yaml");
		fout << root;
	}
	std::filesystem::create_directory(std::format("{}/frames", g_OutputDir));
	std::filesystem::remove_all(std::format("{}/render", g_OutputDir));
	std::filesystem::create_directory(std::format("{}/render", g_OutputDir));
	g_VolFout.open(g_OutputDir + "/volume.txt");
}

inline double WriteFrame(const GridBasedData<2, double> &phi, const int frame)
{
	double volume = 0;
	{ // Write the last frame.
		std::ofstream fout(g_OutputDir + "/end_frame.txt");
		fout << frame + 1 << std::endl;
	}
	const std::string frameDir = std::format("{}/frames/{}", g_OutputDir, frame);
	std::filesystem::create_directory(frameDir);
	{ // Write contour.
		std::ofstream fout(frameDir + "/contour.mesh", std::ios::binary);
		auto mesh = Contourer<2, true>::contour(phi);
		mesh.write(fout);
		mesh.writeOBJ(std::format("{}/render/{}.obj", g_OutputDir, frame));
		for (size_t i = 0; i < mesh.indices.size(); i += 3) {
			const Vector2d a = mesh.positions[mesh.indices[i]];
			const Vector2d b = mesh.positions[mesh.indices[i + 1]];
			const Vector2d c = mesh.positions[mesh.indices[i + 2]];
			const Vector2d ab = b - a, ac = c - a;
			volume += (ab[0] * ac[1] - ab[1] * ac[0]) / 2;
		}
	}
	return volume;
}

inline void ParseArgs(int argc, char *argv[]) {
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("output", 'o', "the output directory", "output");
	parser->addArgument<int>("test", 't', "the test case index", 0);
	parser->addArgument<std::string>("mode", 'm', "the mode of advection");
	parser->addArgument<bool>("refmap", 'r', "enable the reference map", false);
	parser->addArgument<uint>("end", 'e', "the end frame (excluding)");
	parser->addArgument<int>("scale", 's', "the scale of grid");
	parser->addArgument<int>("number", 'n', "the number of PLS particles", 16);
	parser->parse(argc, argv);

	g_OutputDir   = std::any_cast<std::string>(parser->getValueByName("output"));
	g_TestCaseIdx = std::any_cast<int>(parser->getValueByName("test"));
	g_Mode        = std::any_cast<std::string>(parser->getValueByName("mode"));
	g_RefMap      = std::any_cast<bool>(parser->getValueByName("refmap"));
	g_EndFrame    = std::any_cast<uint>(parser->getValueByName("end"));
	g_Scale       = std::any_cast<int>(parser->getValueByName("scale"));
	g_NPPC        = std::any_cast<int>(parser->getValueByName("number"));

	if (!Linear() && !Quadratic() && !Cubic() && !BFECC() && !Hermite() && !Particle()) {
		std::cerr << "Error: [main] emcountered invalid mode." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

inline void Transfer(GridBasedData<2, double> &phiEx, const VectorField<2> &xiView, const ScalarField<2> &phi0View) {
	parallelForEach(phiEx.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phiEx.grid.position(cell);
		phiEx[cell] = phi0View(xiView(pos));
	});
}

inline void Transfer(GridBasedData<2, double> &phi, GridBasedData<2, Vector<2, double>> &psi, const VectorField<2> &xiView, const ScalarField<2> &phi0View) {
	parallelForEach(phi.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phi.grid.position(cell);
		phi[cell] = phi0View(xiView(pos));
		psi[cell] = xiView.gradient(pos) * phi0View.gradient(xiView(pos));
	});
}

inline void Transfer(GridBasedData<2, double> &phiEx, const ScalarField<2> &phiView) {
	parallelForEach(phiEx.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phiEx.grid.position(cell);
		phiEx[cell] = phiView(pos);
	});
}

inline double ComputeError(const GridBasedData<2, double> &phi, const TestField &testField, const int bandSize = 3) {
	double err = 0;
	int cnt = 0;
	double dx = phi.grid.spacing;
	forEach(phi.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phi.grid.position(cell);
		const auto ck = phi[cell];
		const auto gt = testField.getLevelSetGt(pos);
		if (bandSize < 0 || std::abs(gt) <= bandSize * dx) {
			err += std::abs(phi[cell] - gt);
			cnt++;
		}
	});
	return err / cnt;
}

inline double Solve(const bool log) {
	// Initialize data structures;
	auto sGrid  = StaggeredGrid<2>(0, c_Length / g_Scale, g_Scale * Vector2i::Ones());
	auto phi0   = GridBasedData<2, double>(sGrid.cellGrid);
	auto psi0   = GridBasedData<2, Vector2d>(sGrid.cellGrid);
	auto phi    = GridBasedData<2, double>(sGrid.cellGrid);
	auto psi    = GridBasedData<2, Vector2d>(sGrid.cellGrid);
	auto xi     = GridBasedData<2, Vector2d>(sGrid.cellGrid);
	auto pi     = GridBasedData<2, Matrix2d>(sGrid.cellGrid);
	auto gridEx = Refine::get(sGrid.cellGrid, 4);
	auto phiEx  = GridBasedData<2, double>(gridEx);
	auto gridFd = Grid<2>(c_Length / 800, 800 * Vector2i::Ones(), sGrid.cellGrid.origin);
	auto phiFd  = GridBasedData<2, double>(gridFd);
	// Initialize the velocity field.
	std::unique_ptr<TestField> velocityPtr;
	if (g_TestCaseIdx == 0) {
		velocityPtr = std::make_unique<Rotation>();
	} else if (g_TestCaseIdx == 1) {
		velocityPtr = std::make_unique<Twist>();
	} else if (g_TestCaseIdx == 2) {
		velocityPtr = std::make_unique<Bell>();
	}
	// Initialize the level set.
	parallelForEach(phi.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phi.grid.position(cell);
		auto [tphi, tpsi] = velocityPtr->getLevelSetInit(pos);
		phi0[cell] = tphi, psi0[cell] = tpsi;
		phi[cell] = tphi, psi[cell] = tpsi;
		xi[cell] = pos, pi[cell] = Matrix2d::Identity();
	});
	parallelForEach(phiEx.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phiEx.grid.position(cell);
		phiEx[cell] = velocityPtr->getLevelSetInit(pos).first;
	});
	// Initialize the initial level-set view.
	auto phi0ViewLinear    = ScalarFieldView<LinearIntrpl<2>>(phi0);
	auto phi0ViewQuadratic = ScalarFieldView<QuadraticIntrpl<2>>(phi0);
	auto phi0ViewCubic     = ScalarFieldView<CubicIntrpl<2>>(phi0);
	auto phi0ViewHermite   = ScalarFieldView<HermiteIntrpl<2, double>>(phi0, psi0);
	// Initialize particles.
	std::vector<Vector2d> pPos;
	std::vector<double>   pRad;
	std::vector<int>      pSign;
	if (Particle()) {
		PLS::resample(phi, pPos, pRad, pSign, g_NPPC);
	}
	// Prepare logger.
	PrepareWritter(sGrid.radius());
	auto volume0 = WriteFrame(phiEx, 0);
	auto volume  = volume0;
	if (log) {
		g_VolFout << std::format("{:.6e}", volume0) << std::endl;
		std::cout << std::format("Initial volume: {:.6e}", volume0) << std::endl;
	}

	const auto beginTime = std::chrono::steady_clock::now();
	for (uint i = 1; i < g_EndFrame; i++) {
		g_Time = i * c_Dt;
		if (log) std::cout << std::format("[{:>4}] Advect by {} scheme, {} reference map... ", i, g_Mode, g_RefMap ? "with" : "without") << std::flush;

		if (g_RefMap) {
			if (Linear()) {
				Advection::solve<SemiLagrangian<2, 3, LinearIntrpl<2>>>(xi, *velocityPtr, c_Dt);
				Transfer(phiEx, VectorFieldView<LinearIntrpl<2>>(xi), phi0ViewHermite);
			} else if (Quadratic()) {
				Advection::solve<SemiLagrangian<2, 3, EnhancedLinearIntrpl<2, Vector2d>>>(xi, *velocityPtr, c_Dt);
				Transfer(phiEx, VectorFieldView<EnhancedLinearIntrpl<2, Vector2d>>(xi), phi0ViewHermite);
			} else if (Cubic()) {
				Advection::solve<SemiLagrangian<2, 3, CubicIntrpl<2>>>(xi, *velocityPtr, c_Dt);
				Transfer(phiEx, VectorFieldView<CubicIntrpl<2>>(xi), phi0ViewHermite);
			} else if (BFECC()) {
				Advection::solve<MacCormack<2, 3>>(xi, *velocityPtr, c_Dt);
				Transfer(phiEx, VectorFieldView<LinearIntrpl<2>>(xi), phi0ViewHermite);
			} else if (Hermite()) {
				Advection::solve<SemiLagrangian<2, 3, HermiteIntrpl<2, Vector2d>>>(xi, pi, *velocityPtr, c_Dt);

				// if (g_Mode == "hnewton" && i % 50 == 1) {
				if (g_Mode == "hnewton") {
					auto points = GridBasedData<2, Vector2d>(sGrid.cellGrid);
					Transfer(phi, psi, VectorFieldView<HermiteIntrpl<2, Vector2d>>(xi, pi), phi0ViewHermite);
					parallelForEach(sGrid.cellGrid, [&](const Vector2i &cell) {
						phi0[cell] = phi[cell];
					});
					Reinitialization::solve<FastMarching<2>>(phi0, -1);
					phi0ViewLinear.reset();
					parallelForEach(points.grid, [&](const Vector2i &cell) {
						Vector2d pos = points.grid.position(cell);
						points[cell] = pos - phi0ViewLinear(pos) * phi0ViewLinear.gradient(pos);
					});
					HybridNewton<2> reinizer(sGrid.cellGrid);
					reinizer.perform(phi, psi, points, 30, 0.9);
					parallelForEach(sGrid.cellGrid, [&](const Vector2i &cell) {
						phi0[cell] = phi[cell];
						psi0[cell] = psi[cell];
						xi[cell] = sGrid.cellCenter(cell);
						pi[cell] = Matrix2d::Identity();
					});
					phi0ViewHermite.reset();
				}

				if (g_Mode == "hfmm" && i % 50 == 1) {
					Transfer(phi, psi, VectorFieldView<HermiteIntrpl<2, Vector2d>>(xi, pi), phi0ViewHermite);
					parallelForEach(sGrid.cellGrid, [&](const Vector2i &cell) {
						phi0[cell] = phi[cell];
					});
					Reinitialization::solve<FastMarching<2>>(phi0, -1);
					Derivatives::computeFirst<2>(phi0, psi0);
					parallelForEach(sGrid.cellGrid, [&](const Vector2i &cell) {
						xi[cell] = sGrid.cellCenter(cell);
						pi[cell] = Matrix2d::Identity();
					});
					phi0ViewHermite.reset();
				}

				Transfer(phiEx, VectorFieldView<HermiteIntrpl<2, Vector2d>>(xi, pi), phi0ViewHermite);
			}
		} else {
			if (Linear()) {
				Advection::solve<SemiLagrangian<2, 3, LinearIntrpl<2>>>(phi, *velocityPtr, c_Dt);
				Transfer(phiEx, ScalarFieldView<LinearIntrpl<2>>(phi));
			} else if (Quadratic()) {
				Advection::solve<SemiLagrangian<2, 3, EnhancedLinearIntrpl<2, double>>>(phi, *velocityPtr, c_Dt);
				Transfer(phiEx, ScalarFieldView<EnhancedLinearIntrpl<2, double>>(phi));
			} else if (Cubic()) {
				Advection::solve<SemiLagrangian<2, 3, CubicIntrpl<2>>>(phi, *velocityPtr, c_Dt);
				Transfer(phiEx, ScalarFieldView<CubicIntrpl<2>>(phi));
			} else if (BFECC()) {
				Advection::solve<MacCormack<2, 3>>(phi, *velocityPtr, c_Dt);
				Transfer(phiEx, ScalarFieldView<LinearIntrpl<2>>(phi));
			} else if (Hermite()) {
				Advection::solve<SemiLagrangian<2, 3, HermiteIntrpl<2, double>>>(phi, psi, *velocityPtr, c_Dt);
				Transfer(phiEx, ScalarFieldView<HermiteIntrpl<2, double>>(phi, psi));
			} else if (Particle()) {
				Advection::solve<SemiLagrangian<2, 3, LinearIntrpl<2>>>(phi, *velocityPtr, c_Dt);
				PLS::advectParticles(pPos, *velocityPtr, c_Dt);
				PLS::correctLevelSet(phi, pPos, pRad, pSign);
				Reinitialization::solve<FastMarching<2>>(phi, -1);
				PLS::correctLevelSet(phi, pPos, pRad, pSign);
				PLS::reradius(phi, pPos, pRad, pSign);
				Transfer(phiEx, ScalarFieldView<LinearIntrpl<2>>(phi));
			}
		}

		volume = WriteFrame(phiEx, i);
		if (log) {
			g_VolFout << std::format("{:.6e}", volume) << std::endl;
			std::cout << std::format("Volume: {:.6e}", volume) << std::endl;
		}
	}
	const auto endTime = std::chrono::steady_clock::now();

	if (g_RefMap) {
		if (Linear()) {
			Transfer(phiFd, VectorFieldView<LinearIntrpl<2>>(xi), phi0ViewHermite);
		} else if (Quadratic()) {
			Transfer(phiFd, VectorFieldView<EnhancedLinearIntrpl<2, Vector2d>>(xi), phi0ViewHermite);
		} else if (Cubic()) {
			Transfer(phiFd, VectorFieldView<CubicIntrpl<2>>(xi), phi0ViewHermite);
		} else if (BFECC()) {
			Transfer(phiFd, VectorFieldView<LinearIntrpl<2>>(xi), phi0ViewHermite);
		} else if (Hermite()) {
			Transfer(phiFd, VectorFieldView<HermiteIntrpl<2, Vector2d>>(xi, pi), phi0ViewHermite);
		}
	} else {
		if (Linear()) {
			Transfer(phiFd, ScalarFieldView<LinearIntrpl<2>>(phi));
		} else if (Quadratic()) {
			Transfer(phiFd, ScalarFieldView<EnhancedLinearIntrpl<2, double>>(phi));
		} else if (Cubic()) {
			Transfer(phiFd, ScalarFieldView<CubicIntrpl<2>>(phi));
		} else if (BFECC()) {
			Transfer(phiFd, ScalarFieldView<LinearIntrpl<2>>(phi));
		} else if (Hermite()) {
			Transfer(phiFd, ScalarFieldView<HermiteIntrpl<2, double>>(phi, psi));
		} else if (Particle()) {
			Transfer(phiFd, ScalarFieldView<LinearIntrpl<2>>(phi));
		}
	}
	const auto err = ComputeError(phiFd, *velocityPtr);
	if (log) {
		std::cout << std::format("dx: {:.6e}", sGrid.spacing) << std::endl;
		std::cout << std::format("Rel. Volume Loss: {:.6e}", (volume0 - volume) / volume0) << std::endl;
		std::cout << std::format("Err. L1 norm: {:.6e}", err) << std::endl;
		std::cout << std::format("Time: {:>10.3f}s used", std::chrono::duration<double>(endTime - beginTime).count()) << std::endl;
	}
	return err;
}

int main(int argc, char *argv[]) {
	ParseArgs(argc, argv);
	if (g_Scale > 0) {
		Solve(true);
	} else {
		for (int s : { 40, 50, 75, 100, 125, 150 }) {
			g_Scale = s;
			std::cout << std::format("{:.6e}", Solve(false)) << std::endl;
		}
	}
	return 0;
}
