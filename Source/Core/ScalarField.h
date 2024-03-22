#pragma once

#include "Interpolator.h"

namespace PhysX {

template <int Dim>
class ScalarField
{
	DECLARE_DIM_TYPES(Dim)

public:

	virtual ~ScalarField() = default;

	virtual double operator()(const VectorDd &pos) const = 0;
	virtual VectorDd gradient(const VectorDd &pos) const = 0;
};

template <int Dim, typename Data, typename Intrpl> class ScalarFieldWrapper;
template <int Dim, typename Data, typename Intrpl> class DrvScalarFieldWrapper;

template <int Dim, typename Intrpl>
class ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, Intrpl> : public ScalarField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const GridBasedData<Dim, double> &_gbd;
	std::unique_ptr<Intrpl> _intrpl;

public:

	ScalarFieldWrapper(const GridBasedData<Dim, double> &gbd, const bool init = false) : _gbd(gbd) { if (init) reset(); }

	void reset() { _intrpl = makeIntrpl<Intrpl>(_gbd); }

	virtual double operator()(const VectorDd &pos) const override { return _intrpl->interpolate(_gbd, pos); }

	virtual VectorDd gradient(const VectorDd &pos) const
	{
		VectorDd grad = VectorDd::Zero();
		for (const auto [coord, weight] : LinearIntrpl<Dim>::wtPoints(_gbd.grid, pos))
			grad += Derivatives::getFirst<2>(_gbd, coord) * weight;
		return grad;
	}
};

template <int Dim, typename Intrpl>
class DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, Intrpl> : public ScalarField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const GridBasedData<Dim, double> &_gbd;
	const GridBasedData<Dim, VectorDd> &_drv;
	std::unique_ptr<Intrpl> _intrpl;

public:

	DrvScalarFieldWrapper(const GridBasedData<Dim, double> &gbd, const GridBasedData<Dim, VectorDd> &drv, const bool init = false) : _gbd(gbd), _drv(drv) { if (init) reset(); }

	void reset() { _intrpl = makeIntrpl<Intrpl>(_gbd, _drv); }

	virtual double operator()(const VectorDd &pos) const override { return _intrpl->interpolate(_gbd, _drv, pos); }
	virtual VectorDd gradient(const VectorDd &pos) const override { return _intrpl->interpolateGradient(_gbd, _drv, pos); }
};

template <typename Intrpl, int Dim>
inline auto ScalarFieldView(const GridBasedData<Dim, double> &gbd, const bool init = true) { return ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, Intrpl>(gbd, init); }

template <typename Intrpl, int Dim>
inline auto ScalarFieldView(const GridBasedData<Dim, double> &gbd, const GridBasedData<Dim, Vector<Dim, double>> &drv, const bool init = true) { return DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, Intrpl>(gbd, drv, init); }

}
