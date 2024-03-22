#pragma once

#include "Grid.h"
#include "IO.h"

#include <algorithm>
#include <numeric>
#include <vector>

namespace PhysX {

template <int Dim, typename Type>
class GridBasedData
{
	DECLARE_DIM_TYPES(Dim)

public:

	const Grid<Dim> &grid;

protected:

	std::vector<Type> _vec;

public:

	GridBasedData(const Grid<Dim> &grid_, const Type &value = Zero<Type>()) : grid(grid_), _vec(grid.count(), value) { }

	GridBasedData &operator=(const GridBasedData &other)
	{
		assert(grid == other.grid);

		_vec = other._vec;
		return *this;
	}

	Type *data() { return _vec.data(); }
	const Type *data() const { return _vec.data(); }

	size_t index(const VectorDi &coord) const { return grid.index(coord); }
	VectorDi coordinate(const size_t index) const { return grid.coordinate(index); }

	Type &operator[](const size_t index) { return _vec[index]; }
	const Type &operator[](const size_t index) const { return _vec[index]; }
	Type &operator[](const VectorDi &coord) { return _vec[grid.index(coord)]; }
	const Type &operator[](const VectorDi &coord) const { return _vec[grid.index(coord)]; }
	const Type &at(const VectorDi &coord) const { return _vec[grid.index(grid.clamp(coord))]; }

	void setConstant(const Type &value) { std::fill(_vec.begin(), _vec.end(), value); }
	void setZero() { setConstant(Zero<Type>()); }

	template <typename AccType = Type>
	AccType sum() const { return std::accumulate(_vec.begin(), _vec.end(), Zero<AccType>()); }

	double normMax() const
	{
		if constexpr (requires (Type a) { { a.squaredNorm() } -> std::convertible_to<double>; }) {
			double squaredNormMax = 0;
			for (const auto &val : _vec) {
				squaredNormMax = std::max(squaredNormMax, val.squaredNorm());
			}
			return std::sqrt(squaredNormMax);
		}
		else {
			auto minmax = std::minmax_element(_vec.begin(), _vec.end());
			return std::max(std::abs(*minmax.first), std::abs(*minmax.second));
		}
	}

	template <int Order>
	auto stencil(const VectorDi &lower, const int axis) const
	{
		std::array<Type, Order + 1> vals;
		for (int i = 0; i <= Order; i++)
			vals[i] = at(lower + VectorDi::Unit(axis) * i);
		return vals;
	}

	template <int Order>
	auto stencil(const VectorDi &lower) const
	{
		if constexpr (Dim == 2) {
			std::array<Type, (Order + 1) * (Order + 1)> vals;
			for (int j = 0; j <= Order; j++)
				for (int i = 0; i <= Order; i++)
					vals[j * (Order + 1) + i] = at(lower + VectorDi(i, j));
			return vals;
		}
		else {
			std::array<Type, (Order + 1) * (Order + 1) * (Order + 1)> vals;
			for (int k = 0; k <= Order; k++)
				for (int j = 0; j <= Order; j++)
					for (int i = 0; i <= Order; i++)
						vals[(k * (Order + 1) + j) * (Order + 1) + i] = at(lower + VectorDi(i, j, k));
			return vals;
		}
	}

	auto asVectorXd() { return Eigen::Map<VectorXd, Eigen::Aligned>(reinterpret_cast<double *>(_vec.data()), _vec.size() * (sizeof(Type) / sizeof(double))); }
	auto asVectorXd() const { return Eigen::Map<const VectorXd, Eigen::Aligned>(reinterpret_cast<const double *>(_vec.data()), _vec.size() * (sizeof(Type) / sizeof(double))); }

	void load(std::istream &in) { IO::readArray(in, _vec.data(), _vec.size()); }
	void save(std::ostream &out) const { IO::writeArray(out, _vec.data(), _vec.size()); }
};

}
