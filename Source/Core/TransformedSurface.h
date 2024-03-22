#pragma once

#include <Eigen/Geometry>

#include "Surface.h"

namespace PhysX {

using Eigen::AngleAxisf;
using Eigen::AngleAxisd;

template <int Dim>
class MovedSurface : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	VectorDd _displacement;
	const Surface<Dim> &_surface;

public:

	MovedSurface(const Surface<Dim> &surface, const VectorDd &displacement) : _surface(surface), _displacement(displacement) { }

	virtual VectorDd closestPosition(const VectorDd &pos) const override { return _surface.closestPosition(pos - _displacement) + _displacement; }
	virtual VectorDd closestNormal(const VectorDd &pos) const override { return _surface.closestNormal(pos - _displacement); }
	virtual double signedDistance(const VectorDd &pos) const override { return _surface.signedDistance(pos - _displacement); }
};

template <int Dim> class RotatedSurface;

template <>
class RotatedSurface<3> : public Surface<3>
{
protected:

	Vector3d _center;
	AngleAxisd _rotaa;
	const Surface<3> &_surface;

public:

	RotatedSurface(const Surface<3> &surface, const Vector3d &center, const AngleAxisd &rotaa) : _surface(surface), _center(center), _rotaa(rotaa) { }

	Vector3d getRotatedPosition(const Vector3d &pos) const { return _center + _rotaa * (pos - _center); }
	Vector3d getInverseRotatedPosition(const Vector3d &pos) const { return _center + _rotaa.inverse() * (pos - _center); }

	virtual Vector3d closestPosition(const Vector3d &pos) const override { return getInverseRotatedPosition(_surface.closestPosition(getRotatedPosition(pos))); }
	virtual Vector3d closestNormal(const Vector3d &pos) const override { return _rotaa.inverse() * _surface.closestNormal(getRotatedPosition(pos)); }
	virtual double signedDistance(const Vector3d &pos) const override { return _surface.signedDistance(getRotatedPosition(pos)); }
};

template <int Dim>
class ScaledSurface : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	VectorDd _center;
	double _scalef;
	const Surface<Dim> &_surface;

public:

	ScaledSurface(const Surface<Dim> &surface, const VectorDd &center, const double &scalef) : _surface(surface), _center(center), _scalef(scalef) { }

	VectorDd getScaledPosition(const VectorDd &pos) const { return _center + 1 / _scalef * (pos - _center); }
	VectorDd getInverseScaledPosition(const VectorDd &pos) const { return _center + _scalef * (pos - _center); }

	virtual VectorDd closestPosition(const VectorDd &pos) const override { return getInverseScaledPosition(_surface.closestPosition(getScaledPosition(pos))); }
	virtual VectorDd closestNormal(const VectorDd &pos) const override { return _surface.closestNormal(getScaledPosition(pos)); }
	virtual double signedDistance(const VectorDd &pos) const override { return _surface.signedDistance(getScaledPosition(pos)); }
};

template <int Dim>
class ReversedSurface : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Surface<Dim> &_surface;

public:

	ReversedSurface(const Surface<Dim> &surface) : _surface(surface) { }

	virtual VectorDd closestPosition(const VectorDd &pos) const override { return _surface.closestPosition(pos); }
	virtual VectorDd closestNormal(const VectorDd &pos) const override { return -_surface.closestNormal(pos); }
	virtual double signedDistance(const VectorDd &pos) const override { return -_surface.signedDistance(pos); }
};

} // namespace PhysX
