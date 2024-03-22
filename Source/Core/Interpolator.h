#pragma once

#include "EnhancedIntrpl.h"
#include "HermiteIntrpl.h"
#include "LagrangeIntrpl.h"

#include <memory>

namespace PhysX {

template <typename Type, int Dim>
concept BasicIntrpl = requires (const Grid<Dim> &grid, const Vector<Dim, double> &pos) { Type::wtPoints(grid, pos); };

template <typename Intrpl, int Dim, typename Type>
std::unique_ptr<Intrpl> makeIntrpl(const GridBasedData<Dim, Type> &gbd)
{
	if constexpr (BasicIntrpl<Intrpl, Dim>) return std::make_unique<Intrpl>();
	else return std::make_unique<Intrpl>(gbd);
}

template <typename Intrpl, int Dim, typename Type>
std::unique_ptr<Intrpl> makeIntrpl(const GridBasedData<Dim, Type> &gbd, const GridBasedData<Dim, Upgrade<Type, Dim>> &drv)
{
	if constexpr (BasicIntrpl<Intrpl, Dim>) return std::make_unique<Intrpl>();
	else return std::make_unique<Intrpl>(gbd, drv);
}

}
