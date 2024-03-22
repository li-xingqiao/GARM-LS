#pragma once

#include "ScalarField.h"

namespace PhysX {

template <int Dim, typename Field, typename Intrpl>
class RefMapWrapper : public ScalarField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Field &_refField;
	const GridBasedData<Dim, VectorDd> &_refMap;
	const GridBasedData<Dim, MatrixDd> &_refMapGrad;
	std::unique_ptr<Intrpl> _intrpl;

public:

	RefMapWrapper(const Field &refField, const GridBasedData<Dim, VectorDd> &refMap, const GridBasedData<Dim, MatrixDd> &refMapGrad, const bool init = false) : _refField(refField), _refMap(refMap), _refMapGrad(refMapGrad) { if (init) reset(); }

	void reset() { _intrpl = makeIntrpl<Intrpl>(_refMap, _refMapGrad); }

	virtual double operator()(const VectorDd &pos) const override
	{
		const VectorDd orig = _intrpl->interpolate(_refMap, _refMapGrad, pos);
		return _refField(orig);
	}

	virtual VectorDd gradient(const VectorDd &pos) const override
	{
		const VectorDd orig = _intrpl->interpolate(_refMap, _refMapGrad, pos);
		const MatrixDd grad = _intrpl->interpolateGradient(_refMap, _refMapGrad, pos);
		return grad * _refField.gradient(orig);
		// const double dx = _refMap.grid.spacing;
		// const double invDx = _refMap.grid.invSpacing;
		// VectorDd res;
		// for (int axis = 0; axis < Dim; ++axis) {
		// 	res[axis] = this->operator()(pos + 0.2 * dx * VectorDd::Unit(axis)) - this->operator()(pos - 0.2 * dx * VectorDd::Unit(axis));
		// 	res[axis] *= 0.5 * invDx * 5;
		// }
		// return res;
	}
};

template <int Dim, typename Field, typename Intrpl>
class RefMapOnlyWrapper : public ScalarField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Field &_refField;
	const GridBasedData<Dim, VectorDd> &_refMap;
	std::unique_ptr<Intrpl> _intrpl;

public:

	RefMapOnlyWrapper(const Field &refField, const GridBasedData<Dim, VectorDd> &refMap, const bool init = false) : _refField(refField), _refMap(refMap) { if (init) reset(); }

	void reset() { _intrpl = makeIntrpl<Intrpl>(_refMap); }

	virtual double operator()(const VectorDd &pos) const override
	{
		const VectorDd orig = _intrpl->interpolate(_refMap, pos);
		return _refField(orig);
	}

	virtual VectorDd gradient(const VectorDd &pos) const override
	{
		const VectorDd orig = _intrpl->interpolate(_refMap, pos);
		MatrixDd grad = MatrixDd::Zero();
		for (const auto [coord, weight] : LinearIntrpl<Dim>::wtPoints(_refMap.grid, pos))
			grad += Derivatives::getFirst<2>(_refMap, coord) * weight;
		return grad * _refField.gradient(orig);
	}
};

template <typename Intrpl, typename Field, int Dim>
inline auto RefMapFieldView(const Field &rf, const GridBasedData<Dim, Vector<Dim, double>> &rm, const bool init = true) { return RefMapOnlyWrapper<Dim, Field, Intrpl>(rf, rm, init); }

template <typename Intrpl, typename Field, int Dim>
inline auto RefMapFieldView(const Field &rf, const GridBasedData<Dim, Vector<Dim, double>> &rm, const GridBasedData<Dim, Matrix<Dim, double>> &rmg, const bool init = true) { return RefMapWrapper<Dim, Field, Intrpl>(rf, rm, rmg, init); }

} // namespace PhysX
