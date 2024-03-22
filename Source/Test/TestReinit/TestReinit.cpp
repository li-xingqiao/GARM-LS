#include "ArgsParser.h"
#include "Contourer.h"
#include "ImplicitSurface.h"
#include "Reinitialization.h"
#include "HybridNewton.h"
#include "LambdaCorr.h"
#include "StaggeredGrid.h"
#include "Yaml.h"

#include <omp.h>

#include <chrono>
#include <filesystem>
#include <numbers>

using namespace PhysX;

constexpr double length = 10;

inline auto BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("output", 'o', "the output directory", "output");
	parser->addArgument<int>("test", 't', "the test case index", 0);
	parser->addArgument<std::string>("mode", 'm', "the mode of reinitialization");
	parser->addArgument<uint>("end", 'e', "the end frame (excluding)", 201);
	parser->addArgument<double>("cfl", 'c', "the CFL number", .5);
	parser->addArgument<int>("scale", 's', "the scale of grid", 200);
	parser->addArgument<double>("number", 'n', "compute error if |phi| < n * dx", 3);
	return parser;
}

inline void prepareWriter(const std::string outputDir, const double radius)
{
	std::filesystem::create_directory(outputDir);
	{ // Write description.
		YAML::Node root;
		root["dimension"] = 2;
		root["radius"] = float(radius);
		{ // Description of level set.
			YAML::Node node;
			node["name"] = "levelSet";
			node["data_mode"] = "dynamic";
			node["primitive_type"] = "point_list";
			node["indexed"] = false;
			node["color_map"]["enabled"] = true;
			root["objects"].push_back(node);
		}
		{ // Description of contour.
			YAML::Node node;
			node["name"] = "contour";
			node["data_mode"] = "dynamic";
			node["primitive_type"] = "line_list";
			node["material"]["diffuse_albedo"] = Vector4f(0, 0, 1, 1);
			node["indexed"] = true;
			root["objects"].push_back(node);
		}
		std::ofstream fout(outputDir + "/description.yaml");
		fout << root;
	}
}

inline void writeFrame(const GridBasedData<2, double> &phi, const std::string outputDir, const int frame)
{
	{ // Write the last frame.
		std::ofstream fout(outputDir + "/end_frame.txt");
		fout << frame + 1 << std::endl;
	}
	std::filesystem::create_directory(std::format("{}/frames", outputDir));
	const std::string frameDir = std::format("{}/frames/{}", outputDir, frame);
	std::filesystem::create_directory(frameDir);
	{ // Write level set.
		std::ofstream fout(frameDir + "/levelSet.mesh", std::ios::binary);
		IO::writeValue(fout, uint(phi.grid.count()));
		forEach(phi.grid, [&](const Vector2i &cell) {
			const Vector2d pos = phi.grid.position(cell);
			IO::writeValue(fout, pos.template cast<float>().eval());
		});
		forEach(phi.grid, [&](const Vector2i &cell) {
			IO::writeValue(fout, float(phi[cell]));
		});
	}
	{ // Write contour.
		std::ofstream fout(frameDir + "/contour.mesh", std::ios::binary);
		Contourer<2, false>::contour(phi).write(fout);
	}
}

inline void buildCase0(GridBasedData<2, double> &phi, GridBasedData<2, Vector2d> &psi, GridBasedData<2, double> &sdf)
{
	const double r0 = .2313 * length;
	const double maxDist = length / std::numbers::sqrt2 - r0;
	const double factor = maxDist / (std::exp(maxDist) - 1); // to normalize
	parallelForEach(phi.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phi.grid.position(cell);
		const double r = pos.norm();
		phi[cell] = (std::exp(r - r0) - 1) * factor;
		psi[cell] = pos * std::exp(r - r0) / r * factor;
		sdf[cell] = r - r0;
	});
}

inline void buildCase1(GridBasedData<2, double> &phi, GridBasedData<2, Vector2d> &psi, GridBasedData<2, double> &sdf)
{
	const double r0 = length / 4;
	const double factor = 1 / (2 * (length / 2 + r0) * (length / 2 + r0) + .1 * r0 * r0);
	parallelForEach(phi.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phi.grid.position(cell);
		const double r = pos.norm();
		const double f = (pos - Vector2d::Ones() * r0).squaredNorm() + .1 * r0 * r0;
		const double g = r - r0;
		const Vector2d f_ = (pos - Vector2d::Ones() * r0) * 2;
		const Vector2d g_ = pos / r;
		phi[cell] = f * g * factor;
		psi[cell] = (f_ * g + f * g_) * factor;
		sdf[cell] = r - r0;
	});
}

inline void buildCase2(GridBasedData<2, double> &phi, GridBasedData<2, Vector2d> &psi, GridBasedData<2, double> &sdf)
{
	const Vector2d c0 = Vector2d(-1, 1) * length / 4;
	const Vector2d c1 = Vector2d(.5, -.5) * length / 4;
	const double r0 = length / 8;
	const double r1 = std::numbers::pi / 3 * length / 4;
	const double maxDist2 = (Vector2d::Ones() * length / 2 - c1).squaredNorm();
	const double factor = (std::sqrt(maxDist2) - r1) / (std::exp(maxDist2) - std::exp(r1 * r1)); // to normalize
	parallelForEach(phi.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phi.grid.position(cell);
		const double d20 = (pos - c0).squaredNorm();
		const double d21 = (pos - c1).squaredNorm();
		const double f0 = std::exp(d20) - std::exp(r0 * r0);
		const double f1 = std::exp(d21) - std::exp(r1 * r1);
		if (f0 < f1) {
			phi[cell] = f0 * factor;
			psi[cell] = (pos - c0) * 2 * std::exp(d20) * factor;
			sdf[cell] = std::sqrt(d20) - r0;
		}
		else {
			phi[cell] = f1 * factor;
			psi[cell] = (pos - c1) * 2 * std::exp(d21) * factor;
			sdf[cell] = std::sqrt(d21) - r1;
		}
	});
}

inline void buildCase3(GridBasedData<2, double> &phi, GridBasedData<2, Vector2d> &psi, GridBasedData<2, double> &sdf)
{
	const double r0 = .2313 * length;
	const double maxDist = length / std::numbers::sqrt2 - r0;
	const double factor = maxDist / (std::exp(maxDist) - 1); // to normalize
	parallelForEach(phi.grid, [&](const Vector2i &cell) {
		const Vector2d pos = phi.grid.position(cell);
		const double r = pos.norm();
		// phi[cell] = (std::exp(r - r0) - 1) * factor;
		phi[cell] = (r - r0) * 2.;
		psi[cell] = pos * std::exp(r - r0) / r * factor;
		sdf[cell] = r - r0;
	});
}

inline double computeError(const GridBasedData<2, double> &phi, const GridBasedData<2, double> &sdf, const double n)
{
	double error = 0;
	int cnt = 0;
	forEach(phi.grid, [&](const Vector2i &cell) {
		if (n < 0 || abs(phi[cell]) <= n * phi.grid.spacing)
			cnt++, error += std::abs(phi[cell] - sdf[cell]);
	});
	return error / cnt;
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
	const auto end = std::any_cast<uint>(parser->getValueByName("end"));
	const auto cfl = std::any_cast<double>(parser->getValueByName("cfl"));
	const auto scale = std::any_cast<int>(parser->getValueByName("scale"));
	const auto number = std::any_cast<double>(parser->getValueByName("number"));

	const Vector2i resolution = scale * Vector2i::Ones();
	StaggeredGrid<2> sGrid(0, length / scale, resolution);
	GridBasedData<2, double> phi(sGrid.cellGrid);
	GridBasedData<2, Vector2d> psi(sGrid.cellGrid);
	GridBasedData<2, double> sdf(sGrid.cellGrid);

	switch (test) {
	case 0:
		buildCase0(phi, psi, sdf);
		break;
	case 1:
		buildCase1(phi, psi, sdf);
		break;
	case 2:
		buildCase2(phi, psi, sdf);
		break;
	case 3:
		buildCase3(phi, psi, sdf);
		break;
	default:
		std::cerr << "Error: [main] encountered invalid test case." << std::endl;
		std::exit(-1);
	}

	prepareWriter(output, sGrid.radius());
	writeFrame(phi, output, 0);
	std::cout << std::format("Initial error: {:.6e}", computeError(phi, sdf, number)) << std::endl;

	const auto beginTime = std::chrono::steady_clock::now();

	if (mode == "fast-marching") {
		Reinitialization::solve<FastMarching<2>>(phi, -1);
		writeFrame(phi, output, 1);
		std::cout << std::format("Error: {:.6e}", computeError(phi, sdf, number)) << std::endl;
	}
	else if (mode == "upwind" || mode == "eno" || mode == "weno") {
		for (uint i = 1; i < end; i++) {
			std::cout << std::format("[{:>4}] Reinitialize iteratively by {} scheme... ", i, mode) << std::flush;

			if (mode == "upwind") Reinitialization::solve<1>(phi, 1, cfl);
			else if (mode == "eno") Reinitialization::solve<3>(phi, 1, cfl);
			else Reinitialization::solve<5>(phi, 1, cfl);

			writeFrame(phi, output, i);
			std::cout << std::format("Error: {:.6e}", computeError(phi, sdf, number)) << std::endl;
		}
	}
	else if (mode == "gal-newton") {
		std::cout << std::format("[{:>4}] Reinitialize interface cells by {} scheme... ", 1, mode) << std::flush;
		HybridNewton<2> galNewton(sGrid.cellGrid);
		galNewton.perform(phi, psi, 0, 0);
		writeFrame(phi, output, 1);
		std::cout << std::format("Error: {:.6e}", computeError(phi, sdf, number)) << std::endl;

		for (uint i = 2; i < end; i++) {
			std::cout << std::format("[{:>4}] Reinitialize iteratively by {} scheme... ", i, mode) << std::flush;
			galNewton.advect(phi, psi, 1, cfl);
			writeFrame(phi, output, i);
			std::cout << std::format("Error: {:.6e}", computeError(phi, sdf, number)) << std::endl;
		}
	}
	else if (mode == "lambda-corr") {
		LambdaCorr<2, 5> lambdaCorr(sGrid.cellGrid);
		for (uint i = 1; i < end; i++) {
			std::cout << std::format("[{:>4}] Reinitialize iteratively by {} scheme... ", i, mode) << std::flush;

			lambdaCorr.perform(phi, 1, cfl);
			writeFrame(phi, output, i);
			std::cout << std::format("Error: {:.6e}", computeError(phi, sdf, number)) << std::endl;
		}
	}
	else {
		std::cerr << "Error: [main] encountered invalid mode." << std::endl;
		std::exit(-1);
	}

	const auto endTime = std::chrono::steady_clock::now();
	std::cout << std::format("Time: {:>10.3f}s used", std::chrono::duration<double>(endTime - beginTime).count()) << std::endl;

	return 0;
}
