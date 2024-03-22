#include "ArgsParser.h"
#include "LevelSetLiquidBuilder.h"
#include "Simulator.h"

#include <omp.h>

using namespace PhysX;

inline auto BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("output", 'o', "the output directory", "output");
	parser->addArgument<int>("dim", 'd', "the dimension of the simulation", 2);
	parser->addArgument<int>("test", 't', "the test case index", 0);
	parser->addArgument<uint>("begin", 'b', "the begin frame (including)", 0);
	parser->addArgument<uint>("end", 'e', "the end frame (excluding)", 201);
	parser->addArgument<uint>("rate", 'r', "the frame rate (frames per second)", 25);
	parser->addArgument<double>("cfl", 'c', "the CFL number", 0.8);
	parser->addArgument<int>("scale", 's', "the scale of grid", -1);
	parser->addArgument<std::string>("param", 'p', "other parameters", "");
	return parser;
}

int main(int argc, char *argv[])
{
#ifdef _OPENMP
	omp_set_num_threads(std::max(omp_get_num_procs() / 3, 1));
#endif

	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

	const auto output = std::any_cast<std::string>(parser->getValueByName("output"));
	const auto dim = std::any_cast<int>(parser->getValueByName("dim"));
	const auto test = std::any_cast<int>(parser->getValueByName("test"));
	const auto begin = std::any_cast<uint>(parser->getValueByName("begin"));
	const auto end = std::any_cast<uint>(parser->getValueByName("end"));
	const auto rate = std::any_cast<uint>(parser->getValueByName("rate"));
	const auto cfl = std::any_cast<double>(parser->getValueByName("cfl"));
	const auto scale = std::any_cast<int>(parser->getValueByName("scale"));
	const auto parameters = std::any_cast<std::string>(parser->getValueByName("param"));

	std::unique_ptr<Simulation> liquid;
	if (dim == 2)
		liquid = LevelSetLiquidBuilder<2>::build(scale, test, parameters);
	else if (dim == 3)
		liquid = LevelSetLiquidBuilder<3>::build(scale, test, parameters);
	else {
		std::cerr << "Error: [main] encountered invalid dimension." << std::endl;
		std::exit(-1);
	}
	auto simulator = std::make_unique<Simulator>(output, begin, end, rate, cfl, liquid.get());
	simulator->Simulate();

	return 0;
}
