#pragma once

#include "ImplicitSurface.h"
#include "LevelSetLiquid.h"

#include <format>

namespace PhysX {

template <int Dim>
class LevelSetLiquidBuilder
{
	DECLARE_DIM_TYPES(Dim)

public:

	static std::unique_ptr<LevelSetLiquid<Dim>> build(const int scale, const int option, std::string parameters = "");

protected:

	static std::unique_ptr<LevelSetLiquid<Dim>> buildCase0(int scale, std::string parameters = "");
	static std::unique_ptr<LevelSetLiquid<Dim>> buildCase1(int scale);

	static void reportError(const std::string &msg);
};

}
