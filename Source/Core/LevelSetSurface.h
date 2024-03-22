#pragma once

#include "Surface.h"
#include "ScalarField.h"

namespace PhysX {

template <int Dim, typename Intrpl = CubicIntrpl<Dim>>
class LevelSetSurfaceMap : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:
	ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, Intrpl> _levelSetView;

public:
	LevelSetSurfaceMap(const GridBasedData<Dim, double> &levelSet) : _levelSetView(levelSet) { _levelSetView.reset(); }

	virtual VectorDd closestPosition(const VectorDd &pos) const override { return pos - _levelSetView(pos) * _levelSetView.gradient(pos); }
	virtual VectorDd closestNormal(const VectorDd &pos) const override { return _levelSetView.gradient(pos); }
	virtual double distance(const VectorDd &pos) const override { return std::abs(_levelSetView(pos)); }
	virtual double signedDistance(const VectorDd &pos) const override { return _levelSetView(pos); }
	virtual bool isInside(const VectorDd &pos) const override { return _levelSetView(pos) < 0.0; }
};

template <int Dim, typename Intrpl = CubicIntrpl<Dim>>
class LevelSetSurface : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:
	ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, Intrpl> _levelSetView;
	GridBasedData<Dim, double> _levelSetData;

public:
	LevelSetSurface(const GridBasedData<Dim, double> &levelSet) : _levelSetData(levelSet.grid), _levelSetView(_levelSetData) { _levelSetData = levelSet; _levelSetView.reset(); }

	virtual VectorDd closestPosition(const VectorDd &pos) const override { return pos - _levelSetView(pos) * _levelSetView.gradient(pos); }
	virtual VectorDd closestNormal(const VectorDd &pos) const override { return _levelSetView.gradient(pos); }
	virtual double distance(const VectorDd &pos) const override { return std::abs(_levelSetView(pos)); }
	virtual double signedDistance(const VectorDd &pos) const override { return _levelSetView(pos); }
	virtual bool isInside(const VectorDd &pos) const override { return _levelSetView(pos) < 0.0; }
};

template <int Dim, typename S>
void unionSurfaceLevelSet(GridBasedData<Dim, double> &levelSet, const S &surface)
{
	parallelForEach(levelSet.grid, [&](const Vector<Dim, int> &cell) {
		const Vector<Dim, double> pos = levelSet.grid.position(cell);
		levelSet[cell] = std::min(levelSet[cell], surface.signedDistance(pos));
	});
}

template <int Dim, typename S>
void intersectSurfaceLevelSet(GridBasedData<Dim, double> &levelSet, const S &surface)
{
	parallelForEach(levelSet.grid, [&](const Vector<Dim, int> &cell) {
		const Vector<Dim, double> pos = levelSet.grid.position(cell);
		levelSet[cell] = std::max(levelSet[cell], surface.signedDistance(pos));
	});
}

}
