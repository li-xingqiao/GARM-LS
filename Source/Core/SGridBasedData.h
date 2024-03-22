#pragma once

#include "GridBasedData.h"

namespace PhysX {

template <int Dim, typename Type>
class SGridBasedData
{
	DECLARE_DIM_TYPES(Dim)
	using VectorDT = Vector<Dim, Type>;

public:

	const std::array<Grid<Dim>, Dim> &grids;

protected:

	std::array<GridBasedData<Dim, Type>, Dim> _gbds;

public:

	SGridBasedData(const std::array<Grid<Dim>, Dim> &grids_, const VectorDT &value = Zero<VectorDT>()) :
		grids(grids_),
		_gbds(makeGbds(grids, value))
	{ }

	SGridBasedData &operator=(const SGridBasedData &other)
	{
		for (int axis = 0; axis < Dim; axis++)
			_gbds[axis] = other._gbds[axis];
		return *this;
	}

	GridBasedData<Dim, Type> &operator[](const int axis) { return _gbds[axis]; }
	const GridBasedData<Dim, Type> &operator[](const int axis) const { return _gbds[axis]; }

	void setConstant(const VectorDT &value) { for (int axis = 0; axis < Dim; axis++) _gbds[axis].setConstant(value[axis]); }
	void setZero() { setConstant(Zero<VectorDT>()); }

	template <typename AccType>
	AccType sum() const
	{
		if constexpr (Dim == 2) return _gbds[0].template sum<AccType>() + _gbds[1].template sum<AccType>();
		else return _gbds[0].template sum<AccType>() + _gbds[1].template sum<AccType>() + _gbds[2].template sum<AccType>();
	}

	double normMax() const
	{
		if constexpr (Dim == 2) return std::max(_gbds[0].normMax(), _gbds[1].normMax());
		else return std::max({ _gbds[0].normMax(), _gbds[1].normMax(), _gbds[2].normMax() });
	}

	void load(std::istream &in) { for (int axis = 0; axis < Dim; axis++) _gbds[axis].load(in); }
	void save(std::ostream &out) const { for (int axis = 0; axis < Dim; axis++) _gbds[axis].save(out); }

private:

	static auto makeGbds(const std::array<Grid<Dim>, Dim> &grids_, const Vector<Dim, Type> &value = Zero<Vector<Dim, Type>>())
	{
		if constexpr (Dim == 2)
			return std::array {
				GridBasedData<Dim, Type>(grids_[0], value[0]),
				GridBasedData<Dim, Type>(grids_[1], value[1])
			};
		else
			return std::array {
				GridBasedData<Dim, Type>(grids_[0], value[0]),
				GridBasedData<Dim, Type>(grids_[1], value[1]),
				GridBasedData<Dim, Type>(grids_[2], value[2])
			};
	}
};

}
