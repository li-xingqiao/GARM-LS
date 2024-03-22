#pragma once

#include "Interpolator.h"
#include "SGridBasedData.h"

namespace PhysX {

template <int Dim>
class VectorField
{
	DECLARE_DIM_TYPES(Dim)

public:

	virtual ~VectorField() = default;

	virtual VectorDd operator()(const VectorDd &pos) const = 0;
	virtual MatrixDd gradient(const VectorDd &pos) const = 0;
};

template <int Dim, typename Data, typename Intrpl> class VectorFieldWrapper;
template <int Dim, typename Data, typename Intrpl> class DrvVectorFieldWrapper;

template <int Dim, typename Intrpl>
class VectorFieldWrapper<Dim, GridBasedData<Dim, Vector<Dim, double>>, Intrpl> : public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const GridBasedData<Dim, VectorDd> & _gbd;
	std::unique_ptr<Intrpl> _intrpl;

public:

	VectorFieldWrapper(const GridBasedData<Dim, VectorDd> &gbd, const bool init = false) : _gbd(gbd) { if (init) reset(); }

	void reset() { _intrpl = makeIntrpl<Intrpl>(_gbd); }

	virtual VectorDd operator()(const VectorDd &pos) const override { return _intrpl->interpolate(_gbd, pos); }

	virtual MatrixDd gradient(const VectorDd &pos) const override
	{
		MatrixDd grad = MatrixDd::Zero();
		for (const auto [coord, weight] : LinearIntrpl<Dim>::wtPoints(_gbd.grid, pos))
			grad += Derivatives::getFirst<2>(_gbd, coord) * weight;
		return grad;
	}
};

template <int Dim, typename Intrpl>
class VectorFieldWrapper<Dim, SGridBasedData<Dim, double>, Intrpl> : public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const SGridBasedData<Dim, double> &_sgbd;
	std::array<std::unique_ptr<Intrpl>, Dim> _intrpls;

public:

	VectorFieldWrapper(const SGridBasedData<Dim, double> &sgbd, const bool init = false) : _sgbd(sgbd) { if (init) reset(); }

	void reset() { for (int axis = 0; axis < Dim; axis++) _intrpls[axis] = makeIntrpl<Intrpl>(_sgbd[axis]); }

	virtual VectorDd operator()(const VectorDd &pos) const override
	{
		if constexpr (Dim == 2)
			return VectorDd(_intrpls[0]->interpolate(_sgbd[0], pos), _intrpls[1]->interpolate(_sgbd[1], pos));
		else
			return VectorDd(_intrpls[0]->interpolate(_sgbd[0], pos), _intrpls[1]->interpolate(_sgbd[1], pos), _intrpls[2]->interpolate(_sgbd[2], pos));
	}

	virtual MatrixDd gradient(const VectorDd &pos) const override
	{
		MatrixDd grad = MatrixDd::Zero();
		for (int axis = 0; axis < Dim; axis++)
			for (const auto [coord, weight] : LinearIntrpl<Dim>::wtPoints(_sgbd[axis].grid, pos))
				grad.row(axis) += Derivatives::getFirst<2>(_sgbd[axis], coord).transpose() * weight;
		return grad;
	}
};

template <int Dim, typename Intrpl>
class DrvVectorFieldWrapper<Dim, GridBasedData<Dim, Vector<Dim, double>>, Intrpl> : public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const GridBasedData<Dim, VectorDd> &_gbd;
	const GridBasedData<Dim, MatrixDd> &_drv;
	std::unique_ptr<Intrpl> _intrpl;

public:

	DrvVectorFieldWrapper(const GridBasedData<Dim, VectorDd> &gbd, const GridBasedData<Dim, MatrixDd> &drv, const bool init = false) : _gbd(gbd), _drv(drv) { if (init) reset(); }

	void reset() { _intrpl = makeIntrpl<Intrpl>(_gbd, _drv); }

	virtual VectorDd operator()(const VectorDd &pos) const override { return _intrpl->interpolate(_gbd, _drv, pos); }
	virtual MatrixDd gradient(const VectorDd &pos) const override { return _intrpl->interpolateGradient(_gbd, _drv, pos); }
};

template <int Dim, typename Intrpl>
class DrvVectorFieldWrapper<Dim, SGridBasedData<Dim, double>, Intrpl> : public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const SGridBasedData<Dim, double> &_sgbd;
	const SGridBasedData<Dim, VectorDd> &_drv;
	std::array<std::unique_ptr<Intrpl>, Dim> _intrpls;

public:

	DrvVectorFieldWrapper(const SGridBasedData<Dim, double> &sgbd, const SGridBasedData<Dim, VectorDd> &drv, const bool init = false) : _sgbd(sgbd), _drv(drv) { if (init) reset(); }

	void reset() { for (int axis = 0; axis < Dim; axis++) _intrpls[axis] = makeIntrpl<Intrpl>(_sgbd[axis], _drv[axis]); }

	virtual VectorDd operator()(const VectorDd &pos) const override
	{
		VectorDd val;
		for (int axis = 0; axis < Dim; axis++)
			val[axis] = _intrpls[axis]->interpolate(_sgbd[axis], _drv[axis], pos);
		return val;
	}

	virtual MatrixDd gradient(const VectorDd &pos) const override
	{
		MatrixDd grad;
		for (int axis = 0; axis < Dim; axis++)
			grad.row(axis) = _intrpls[axis]->interpolateGradient(_sgbd[axis], _drv[axis], pos).transpose();
		return grad;
	}
};

template <typename Intrpl, int Dim>
inline auto VectorFieldView(const GridBasedData<Dim, Vector<Dim, double>> &gbd, const bool init = true) { return VectorFieldWrapper<Dim, GridBasedData<Dim, Vector<Dim, double>>, Intrpl>(gbd, init); }

template <typename Intrpl, int Dim>
inline auto VectorFieldView(const SGridBasedData<Dim, double> &sgbd, const bool init = true) { return VectorFieldWrapper<Dim, SGridBasedData<Dim, double>, Intrpl>(sgbd, init); }

template <typename Intrpl, int Dim>
inline auto VectorFieldView(const GridBasedData<Dim, Vector<Dim, double>> &gbd, const GridBasedData<Dim, Matrix<Dim, double>> &drv, const bool init = true) { return DrvVectorFieldWrapper<Dim, GridBasedData<Dim, Vector<Dim, double>>, Intrpl>(gbd, drv, init); }

template <typename Intrpl, int Dim>
inline auto VectorFieldView(const SGridBasedData<Dim, double> &sgbd, const SGridBasedData<Dim, Vector<Dim, double>> &drv, const bool init = true) { return DrvVectorFieldWrapper<Dim, SGridBasedData<Dim, double>, Intrpl>(sgbd, drv, init); }
}
