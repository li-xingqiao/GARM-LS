#pragma once

#include "Types.h"

#include <array>
#include <functional>

namespace PhysX {

template <int Dim>
class Grid
{
	DECLARE_DIM_TYPES(Dim)

public:

	const double spacing;
	const double invSpacing;
	const VectorDi size;
	const VectorDd origin;

public:

	Grid(const double spacing_, const VectorDi &size_, const VectorDd &origin_) :
		spacing(spacing_),
		invSpacing(1 / spacing),
		size(size_),
		origin(origin_)
	{ }

	bool operator==(const Grid &other) const { return spacing == other.spacing && size == other.size && origin == other.origin; }

	bool isInside(const VectorDi &coord, const int offset = 0) const { return (coord.array() >= offset).all() && (coord.array() < size.array() - offset).all(); }
	bool isValid(const VectorDi &coord) const { return isInside(coord, 0); }

	size_t count() const { return size.template cast<size_t>().prod(); }
	VectorDd position(const VectorDi &coord) const { return origin + coord.template cast<double>() * spacing; }
	VectorDi clamp(const VectorDi &coord) const { return coord.cwiseMax(0).cwiseMin(size - VectorDi::Ones()); }

	size_t index(const VectorDi &coord) const
	{
		if constexpr (Dim == 2) return coord.x() + size_t(size.x()) * coord.y();
		else return coord.x() + size.x() * (coord.y() + size_t(size.y()) * coord.z());
	}

	VectorDi coordinate(const size_t index) const
	{
		if constexpr (Dim == 2) return VectorDi(index % size.x(), index / size.x());
		else return VectorDi(int(index % size.x()), int(index / size.x() % size.y()), int(index / size.x() / size.y()));
	}

	VectorDi getLinearLower(const VectorDd &pos) const { return ((pos - origin) * invSpacing).array().floor().template cast<int>().matrix(); }
	VectorDi getQuadraticLower(const VectorDd &pos) const { return ((pos - origin) * invSpacing - VectorDd::Ones() / 2).array().floor().template cast<int>().matrix(); }
	VectorDi getCubicLower(const VectorDd &pos) const { return ((pos - origin) * invSpacing).array().floor().template cast<int>().matrix() - VectorDi::Ones(); }
	VectorDd getLowerFrac(const VectorDd &pos, const VectorDi &lower) const { return (pos - origin - lower.template cast<double>() * spacing) * invSpacing; }

	static constexpr int numberOfNeighbors() { return Dim << 1; }
	static constexpr int neighborAxis(const int ord) { return ord >> 1; }
	static constexpr int neighborSide(const int ord) { return ord & 1 ? 1 : -1; }
	static VectorDi neighbor(const VectorDi &coord, const int ord) { return coord + VectorDi::Unit(neighborAxis(ord)) * neighborSide(ord); }
};

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<void(const Vector<Dim, int> &)>>
inline void forEach(const Grid<Dim> &grid, const Func func)
{
	DECLARE_DIM_TYPES(Dim)

	if constexpr (Dim == 2) {
		for (int j = 0; j < grid.size.y(); j++)
			for (int i = 0; i < grid.size.x(); i++)
				func(VectorDi(i, j));
	}
	else {
		for (int k = 0; k < grid.size.z(); k++)
			for (int j = 0; j < grid.size.y(); j++)
				for (int i = 0; i < grid.size.x(); i++)
					func(VectorDi(i, j, k));
	}
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<void(const Vector<Dim, int> &)>>
inline void parallelForEach(const Grid<Dim> &grid, const Func func)
{
	DECLARE_DIM_TYPES(Dim)

	if constexpr (Dim == 2) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int j = 0; j < grid.size.y(); j++)
			for (int i = 0; i < grid.size.x(); i++)
				func(VectorDi(i, j));
	}
	else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int k = 0; k < grid.size.z(); k++)
			for (int j = 0; j < grid.size.y(); j++)
				for (int i = 0; i < grid.size.x(); i++)
					func(VectorDi(i, j, k));
	}
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<void(const int, const Vector<Dim, int> &)>>
inline void forEach(const std::array<Grid<Dim>, Dim> &grids, const Func func)
{
	for (int axis = 0; axis < Dim; axis++)
		forEach(grids[axis], std::bind(func, axis, std::placeholders::_1));
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<void(const int, const Vector<Dim, int> &)>>
inline void parallelForEach(const std::array<Grid<Dim>, Dim> &grids, const Func func)
{
	for (int axis = 0; axis < Dim; axis++)
		parallelForEach<Dim>(grids[axis], std::bind(func, axis, std::placeholders::_1));
}

}
