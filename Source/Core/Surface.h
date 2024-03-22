#pragma once

#include "Types.h"

namespace PhysX {

template <int Dim>
class Surface
{
	DECLARE_DIM_TYPES(Dim)

public:

	virtual ~Surface() = default;

	virtual VectorDd closestPosition(const VectorDd &pos) const { return pos - signedDistance(pos) * closestNormal(pos); }
	virtual VectorDd closestNormal(const VectorDd &pos) const = 0;
	virtual double distance(const VectorDd &pos) const { return std::abs(signedDistance(pos)); }
	virtual double signedDistance(const VectorDd &pos) const = 0;
	virtual bool isInside(const VectorDd &pos) const { return signedDistance(pos) <= 0; }
};

}
