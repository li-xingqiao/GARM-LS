#pragma once

#include "Surface.h"

namespace PhysX {

template <int Dim>
class ImplicitSphere : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDd _center;
	const double _radius;

public:

	ImplicitSphere(const VectorDd &center, const double radius) : _center(center), _radius(radius) { }

	virtual VectorDd closestPosition(const VectorDd &pos) const override { return _center + closestNormal(pos) * _radius; }
	virtual VectorDd closestNormal(const VectorDd &pos) const override { return (pos - _center).normalized(); }
	virtual double signedDistance(const VectorDd &pos) const override { return (pos - _center).norm() - _radius; }
};

template <int Dim>
class ImplicitBox : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDd _center;
	const VectorDd _halfLengths;

public:

	ImplicitBox(const VectorDd &minCorner, const VectorDd &lengths) : _center(minCorner + lengths / 2), _halfLengths(lengths / 2) { }

	virtual VectorDd closestNormal(const VectorDd &pos) const override
	{
		const VectorDd phi = (pos - _center).cwiseAbs() - _halfLengths;
		VectorDd normal;
		if ((phi.array() <= 0).all()) {
			int axis;
			phi.maxCoeff(&axis);
			normal = VectorDd::Unit(axis);
		}
		else normal = phi.cwiseMax(0);
		return normal.cwiseProduct((pos - _center).cwiseSign()).normalized();
	}

	virtual double signedDistance(const VectorDd &pos) const override
	{
		const VectorDd phi = (pos - _center).cwiseAbs() - _halfLengths;
		if ((phi.array() <= 0).all()) return phi.maxCoeff();
		else return phi.cwiseMax(0).norm();
	}
};

template <int Dim>
class ImplicitCylinder : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDd _center;
	const double _radius;
	const double _halfHeight;

public:

	ImplicitCylinder(const VectorDd &baseCenter, const double radius, const double height) : _center(baseCenter + VectorDd::Unit(1) * height * .5), _radius(radius), _halfHeight(height * .5) { }

	virtual VectorDd closestNormal(const VectorDd &pos) const override
	{
		const double phiRadius = ((pos - pos[1] * VectorDd::Unit(1)) - (_center - _center[1] * VectorDd::Unit(1))).norm() - _radius;
		const double phiHeight = std::abs(pos[1] - _center[1]) - _halfHeight;
		if (phiRadius <= 0 && phiHeight <= 0) {
			if (phiRadius > phiHeight)
				return ((pos - pos[1] * VectorDd::Unit(1)) - (_center - _center[1] * VectorDd::Unit(1))).normalized();
			else
				return VectorDd::Unit(1) * (pos[1] > _center[1] ? +1 : -1);
		}
		else if (phiRadius <= 0)
			return VectorDd::Unit(1) * (pos[1] > _center[1] ? +1 : -1);
		else if (phiHeight <= 0)
			return ((pos - pos[1] * VectorDd::Unit(1)) - (_center - _center[1] * VectorDd::Unit(1))).normalized();
		else {
			const VectorDd deltaR = ((pos - pos[1] * VectorDd::Unit(1)) - (_center - _center[1] * VectorDd::Unit(1))).normalized() * phiRadius;
			const VectorDd deltaH = phiHeight * VectorDd::Unit(1) * (pos[1] > _center[1] ? +1 : -1);
			return (deltaR + deltaH).normalized();
		}
	}

	virtual double signedDistance(const VectorDd &pos) const override
	{
		const double phiRadius = ((pos - pos[1] * VectorDd::Unit(1)) - (_center - _center[1] * VectorDd::Unit(1))).norm() - _radius;
		const double phiHeight = std::abs(pos[1] - _center[1]) - _halfHeight;
		if (phiRadius <= 0 && phiHeight <= 0)
			return std::max(phiRadius, phiHeight);
		else {
			const auto square = [](const double x) { return x * x; };
			return std::sqrt(square(std::max(phiRadius, .0)) + square(std::max(phiHeight, .0)));
		}
	}
};

template <int Dim>
class ImplicitPlane : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDd _position;
	const VectorDd _normal;

public:

	ImplicitPlane(const VectorDd &position, const VectorDd &direction) : _position(position), _normal(direction.normalized()) { }

	virtual VectorDd closestNormal(const VectorDd &pos) const override { return _normal; }
	virtual double signedDistance(const VectorDd &pos) const override { return (pos - _position).dot(_normal); }
};

template <int Dim>
class ImplicitEllipsoid : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDd _center;
	const VectorDd _semiAxels;

public:

	ImplicitEllipsoid(const VectorDd &center, const VectorDd &semiAxels) : _center(center), _semiAxels(semiAxels) { }

	virtual VectorDd closestPosition(const VectorDd &pos) const override { return pos / pos.cwiseQuotient(_semiAxels).norm(); } // not accurate solution
	virtual VectorDd closestNormal(const VectorDd &pos) const override { return (pos - closestPosition(pos)).normalized() * (isInside(pos) ? -1 : 1); }
	virtual double distance(const VectorDd &pos) const override { return (pos - closestPosition(pos)).norm(); }
	virtual double signedDistance(const VectorDd &pos) const override { return distance(pos) * (isInside(pos) ? -1 : 1); }
	virtual bool isInside(const VectorDd &pos) const override { return pos.cwiseQuotient(_semiAxels).squaredNorm() <= 1; }
};

template <int Dim>
class ImplicitSemiPlane : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	VectorDd _center;
	double _radius;

public:

	ImplicitSemiPlane(const VectorDd &center, const double radius) : _center(center), _radius(radius) { }

	virtual VectorDd closestNormal(const VectorDd &pos) const override
	{
		VectorDd diff = pos - _center;
		if (diff(0) > 0) {
			return diff(1) > 0 ? VectorDd::Unit(1).eval() : (-1. * VectorDd::Unit(1)).eval();
		}
		else {
			if constexpr (Dim == 3) {
				diff = diff - diff(2) * VectorDd::Unit(2);
			}
			return diff.normalized();
		}
	}

	virtual double signedDistance(const VectorDd &pos) const override
	{
		VectorDd diff = pos - _center;
		if (diff(0) > 0) {
			return std::abs(diff(1)) - _radius;
		}
		else {
			if constexpr (Dim == 3) {
				diff = diff - diff(2) * VectorDd::Unit(2);
			}
			return diff.norm() - _radius;
		}
	}
};

template <int Dim>
class ImplicitCone : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDd _vertex;
	const double _halfAngle;

public:

	ImplicitCone(const VectorDd &vertex, const double hAngle) : _vertex(vertex), _halfAngle(hAngle) { }

	virtual VectorDd closestPosition(const VectorDd &pos) const override {
		VectorDd diff = pos - _vertex;
		VectorDd generatrix = diff;
		generatrix[1] = 0;
		double radius = generatrix.norm();
		generatrix[1] = -radius / std::tan(_halfAngle);
		if (generatrix.dot(diff) < 0) {
			return _vertex;
		} else {
			return _vertex + generatrix * generatrix.dot(diff) / generatrix.squaredNorm();
		}
	}
	virtual VectorDd closestNormal(const VectorDd &pos) const override {
		VectorDd diff = pos - _vertex;
		VectorDd generatrix = diff;
		generatrix[1] = 0;
		double radius = generatrix.norm();
		generatrix[1] = -radius / std::tan(_halfAngle);
		if (generatrix.dot(diff) < 0) {
			return diff.normalized();
		} else {
			generatrix[1] = radius * std::tan(_halfAngle); // generatrix become normal!
			return generatrix.normalized();
		}
	}
	virtual double signedDistance(const VectorDd &pos) const override {
		VectorDd diff = pos - _vertex;
		VectorDd generatrix = diff;
		generatrix[1] = 0;
		double radius = generatrix.norm();
		generatrix[1] = -radius / std::tan(_halfAngle);
		if (generatrix.dot(diff) < 0) {
			return diff.norm();
		} else {
			generatrix[1] = radius * std::tan(_halfAngle); // generatrix become normal!
			return diff.dot(generatrix.normalized());
		}
	}
};

}
