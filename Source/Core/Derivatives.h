#pragma once

#include "GridBasedData.h"

namespace PhysX {

template <int Dim, typename Type>
inline auto makeDrvArray(const GridBasedData<Dim, Type> &gbd)
{
	if constexpr (Dim == 2)
		return std::array { GridBasedData<Dim, Type>(gbd.grid), GridBasedData<Dim, Type>(gbd.grid) };
	else
		return std::array { GridBasedData<Dim, Type>(gbd.grid), GridBasedData<Dim, Type>(gbd.grid), GridBasedData<Dim, Type>(gbd.grid) };
}

template <int Dim, typename Type>
inline auto makeMixedDrvArray(const GridBasedData<Dim, Type> &gbd)
{
	if constexpr (Dim == 2)
		return std::array { GridBasedData<Dim, Type>(gbd.grid) };
	else
		return std::array { GridBasedData<Dim, Type>(gbd.grid), GridBasedData<Dim, Type>(gbd.grid), GridBasedData<Dim, Type>(gbd.grid), GridBasedData<Dim, Type>(gbd.grid) };
}

template <typename DrvType, typename Type>
void setC(DrvType &val, const int axis, const Type &vc)
{
	if constexpr (std::is_scalar_v<Type>)
		val[axis] = vc;
	else
		val.col(axis) = vc;
}

template <typename Type>
Downgrade<Type> getC(const Type &val, const int axis)
{
	if constexpr (std::is_scalar_v<Downgrade<Type>>)
		return val[axis];
	else
		return val.col(axis);
}

}

namespace PhysX::Derivatives {

template <typename Type> constexpr Type firstCentral(const std::array<Type, 3> &v, const double invDx) { return (v[2] - v[0]) / 2 * invDx; }
template <typename Type> constexpr Type firstCentral(const std::array<Type, 5> &v, const double invDx) { return (v[0] - 8 * v[1] + 8 * v[3] - v[4]) / 12 * invDx; }
template <typename Type> constexpr Type secondCentral(const std::array<Type, 3> &v, const double invDx) { return (v[0] - 2 * v[1] + v[2]) * invDx * invDx; }

template <int Order, int Dim, typename Type>
inline Type getFirst(const GridBasedData<Dim, Type> &gbd, const Vector<Dim, int> &coord, const int axis)
{
	DECLARE_DIM_TYPES(Dim)
	static_assert(Order == 2 || Order == 4, "Order must be 2 or 4.");
	return firstCentral(gbd.template stencil<Order>(coord - VectorDi::Unit(axis) * (Order >> 1), axis), gbd.grid.invSpacing);
}

template <int Order, int Dim, typename Type>
inline Type getSecond(const GridBasedData<Dim, Type> &gbd, const Vector<Dim, int> &coord, const int axis1, const int axis2 = -1)
{
	DECLARE_DIM_TYPES(Dim)
	static_assert(Order == 2, "Order must be 2.");
	if (axis2 < 0)
		return secondCentral(gbd.template stencil<Order>(coord - VectorDi::Unit(axis1) * (Order >> 1), axis1), gbd.grid.invSpacing);
	else
		return (getFirst<2, Dim, Type>(gbd, coord + VectorDi::Unit(axis2), axis1) - getFirst<2, Dim, Type>(gbd, coord - VectorDi::Unit(axis2), axis1)) / 2 * gbd.grid.invSpacing;
}

template <int Order, int Dim, typename Type>
inline Upgrade<Type, Dim> getFirst(const GridBasedData<Dim, Type> &gbd, const Vector<Dim, int> &coord)
{
	Upgrade<Type, Dim> drv;
	for (int axis = 0; axis < Dim; axis++)
		setC(drv, axis, getFirst<Order>(gbd, coord, axis));
	return drv;
}

template <int Order, int Dim, typename Type>
inline Upgrade<Upgrade<Type, Dim>, Dim> getSecond(const GridBasedData<Dim, Type> &gbd, const Vector<Dim, int> &coord)
{
	Upgrade<Upgrade<Type, Dim>, Dim> drv;
	for (int axis1 = 0; axis1 < Dim; axis1++)
		for (int axis2 = 0; axis2 < Dim; axis2++)
			drv(axis1, axis2) = getSecond<Order>(gbd, coord, axis1, axis2);
	return drv;
}

template <int Order, int Dim, typename Type>
inline void computeFirst(const GridBasedData<Dim, Type> &gbd, GridBasedData<Dim, Type> &drv, const int axis)
{
	DECLARE_DIM_TYPES(Dim)
	parallelForEach(gbd.grid, [&](const VectorDi &coord) {
		drv[coord] = getFirst<Order>(gbd, coord, axis);
	});
}

template <int Order, int Dim, typename Type>
inline void computeFirst(const GridBasedData<Dim, Type> &gbd, GridBasedData<Dim, Upgrade<Type, Dim>> &drv)
{
	DECLARE_DIM_TYPES(Dim)
	parallelForEach(gbd.grid, [&](const VectorDi &coord) {
		drv[coord] = getFirst<Order>(gbd, coord);
	});
}

template <int Order, int Dim, typename Type>
inline void computeSecond(const GridBasedData<Dim, Type> &gbd, GridBasedData<Dim, Type> &drv, const int axis1, const int axis2 = -1)
{
	DECLARE_DIM_TYPES(Dim)
	parallelForEach(gbd.grid, [&](const VectorDi &coord) {
		drv[coord] = getSecond<Order>(gbd, coord, axis1, axis2);
	});
}

constexpr std::pair<double, double> upwind1(const std::array<double, 3> &v, const double invDx) { return std::make_pair((v[1] - v[0]) * invDx, (v[2] - v[1]) * invDx); }

constexpr std::pair<double, double> eno3(const std::array<double, 7> &v, const double invDx)
{
	std::pair<double, double> drv;
	const std::array<double, 6> d1 = { v[1] - v[0], v[2] - v[1], v[3] - v[2], v[4] - v[3], v[5] - v[4], v[6] - v[5] };
	const std::array<double, 5> d2 = { d1[1] - d1[0], d1[2] - d1[1], d1[3] - d1[2], d1[4] - d1[3], d1[5] - d1[4] };
	for (int i = 0; i < 2; i++) {
		const int j = i - (std::abs(d2[i + 1]) < std::abs(d2[i + 2]));
		const std::array<double, 2> d3 = { d2[j + 2] - d2[j + 1], d2[j + 3] - d2[j + 2] };
		const double dQ1 = d1[i + 2];
		const double dQ2 = d2[j + 2] * (1 - 2 * i) / 2;
		const double dQ3 = d3[std::abs(d3[1]) < std::abs(d3[0])] * (3 * j * j - 1) / 6;
		(i ? drv.second : drv.first) = (dQ1 + dQ2 + dQ3) * invDx;
	}
	return drv;
}

constexpr std::pair<double, double> weno5(const std::array<double, 7> &v, const double invDx)
{
	const auto square = [](const double x) { return x * x; };

	std::pair<double, double> drv;
	for (int i = 0; i < 2; i++) {
		std::array<double, 5> dev = i ?
			std::array { v[6] - v[5], v[5] - v[4], v[4] - v[3], v[3] - v[2], v[2] - v[1] } :
			std::array { v[1] - v[0], v[2] - v[1], v[3] - v[2], v[4] - v[3], v[5] - v[4] };
		double eps = 0;
		for (int j = 0; j < 5; j++) eps = std::max(eps, dev[j] * dev[j]);
		eps = eps * 1e-6 + 1e-99;
		const double phix1 = dev[0] / 3 - dev[1] * 7 / 6 + dev[2] * 11 / 6;
		const double phix2 = -dev[1] / 6 + dev[2] * 5 / 6 + dev[3] / 3;
		const double phix3 = dev[2] / 3 + dev[3] * 5 / 6 - dev[4] / 6;
		const double s1 = square(dev[0] - 2 * dev[1] + dev[2]) * 13 / 12 + square(dev[0] - 4 * dev[1] + 3 * dev[2]) / 4;
		const double s2 = square(dev[1] - 2 * dev[2] + dev[3]) * 13 / 12 + square(dev[1] - dev[3]) / 4;
		const double s3 = square(dev[2] - 2 * dev[3] + dev[4]) * 13 / 12 + square(3 * dev[2] - 4 * dev[3] + dev[4]) / 4;
		const double alpha1 = .1 / square(s1 + eps);
		const double alpha2 = .6 / square(s2 + eps);
		const double alpha3 = .3 / square(s3 + eps);
		const double sum = alpha1 + alpha2 + alpha3;
		(i ? drv.second : drv.first) = (alpha1 * phix1 + alpha2 * phix2 + alpha3 * phix3) / sum * invDx;
	}
	return drv;
}

}
