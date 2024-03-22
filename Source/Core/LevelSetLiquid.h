#pragma once

#include "Boundary.h"
// #include "GridProjector.h"
#include "Projector.h"
#include "ScalarField.h"
#include "Simulation.h"

#include <random>
#include <numbers>

namespace PhysX {

template <int Dim> class LevelSetLiquidBuilder;

template <int Dim>
class LevelSetLiquid : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class LevelSetLiquidBuilder<Dim>;

protected:

	const StaggeredGrid<Dim> _sGrid;

	SGridBasedData<Dim, double> _velocity;
	GridBasedData<Dim, double> _levelSet;

#ifdef _ENABLE_PARTICLE_
    struct Particle {
		VectorDd pos;
		double rad;
		int sign;
	};
	std::vector<Particle> _particles;
#endif

#if defined(_DISABLE_RM_) && !defined(_DISABLE_GA_)
#error GA cannot be enabled when RM is disabled now.
#endif

#ifdef _ENABLE_NEWTON_
	GridBasedData<Dim, VectorDd> _levelSetGrad;
#endif

#ifdef _USE_LERP_
	ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, LinearIntrpl<Dim>> _levelSetView;
#else
	ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, CubicIntrpl<Dim>> _levelSetView;
#endif

#ifndef _DISABLE_RM_
	GridBasedData<Dim, VectorDd> _refMap;
	GridBasedData<Dim, MatrixDd> _refMapGrad; // used for check restart
	GridBasedData<Dim, double> _tmpLevelSet;  // used for rendering
	GridBasedData<Dim, uchar> _advRegion;     // record advection region
	GridBasedData<Dim, double> _oriLevelSet;  // used for extrapolation
	GridBasedData<Dim, double> _refLevelSet;
#ifndef _DISABLE_GA_
	GridBasedData<Dim, VectorDd> _refLevelSetGrad;
	DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>> _refLevelSetView;
#else
	ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, CubicIntrpl<Dim>> _refLevelSetView;
#endif
#endif

	GridBasedData<Dim, double> _velocityDiv;

	bool _enableGravity = true;
	bool _enableRotation = false;
	bool _enableSurfaceTension = false;
	bool _enableBounce = false;
	bool _enableRaining = false;
	bool _enableMovingBoundary = false;
	bool _enableContactAngle = false;
	bool _enableColorOutput = false;
	bool _enablePrettifiedOutput = true;
	bool _enableResume = true;
	bool _enableRefine = false;
	double _elapsedTime = 0.0;   // time elapsed
	double _rainingRate = 10.0;  // droplets per second
	double _rainingCount = 0.5;  // initial droplets
	double _rainingSpeed = 10.0; // speed of droplets
	double _density = 1e3;
	double _surfaceTensionCoeff = 7.28e-2;
	double _surfaceTensionDeltaMin = 2.0;
	double _surfaceTensionDeltaMax = 8.0;
	double _surfaceTensionDeltaCritTime = 1.0;
	double _boundaryExtrapolationRadius = 6.0;
	double _airExtrapolationRadius = 15.0;
	double _contactAngle = std::numbers::pi / 2;

	std::function<VectorDd(double)> _boundaryTransformation;
	std::function<VectorDd(double,VectorDd)> _boundaryVelocity;

	int _reinitCount = 0;

	std::vector<int> _restarted;
	std::vector<double> _volume;
	struct TimeInfo {
		double tadvect = 0.0;
		double tproject = 0.0;
		double trestart = 0.0;
		double textrap = 0.0;
	};
	std::vector<TimeInfo> _profile;
	GridBasedData<Dim, uchar> _reinitmark;
	GridBasedData<Dim, double> _refColor;
	ScalarFieldWrapper<Dim, GridBasedData<Dim, double>, LinearIntrpl<Dim>> _refColorView;
	GridBasedData<Dim, double> _color;

	Boundary<Dim> _boundary;

	// GridProjector<Dim> _projector;
	Projector<Dim> _projector;

	std::unique_ptr<Surface<Dim>> _source;
	GridBasedData<Dim, double> _initialBoundary;

public:

	LevelSetLiquid(const StaggeredGrid<Dim> &sGrid);

	virtual ~LevelSetLiquid() = default;

	virtual double getTimeStep(const uint frameRate, const double stepRate) const override { return stepRate * _sGrid.spacing / _velocity.normMax(); }

	virtual int dimension() const override { return Dim; }
	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	virtual void initialize() override;
	virtual void advance(const double dt) override;

	void unionLevelSet(const Surface<Dim> &surface);
	void intersectLevelSet(const Surface<Dim> &surface);

protected:

	void advectFields(const double dt);
	void applyBodyForces(const double dt);
	void applyFluidSources(const double dt);
	void transformBoundary();
	void enforceBoundaryLevelSet(GridBasedData<Dim, double> &sdf, const double tolerance = 0);
	void enforceBoundaryLevelSet(GridBasedData<Dim, double> &sdf, GridBasedData<Dim, Vector<Dim, double>> &drv, const double tolerance = 0);
	void projectVelocity(const double dt);
	void applySurfaceTension(const double dt);
	void trickLittleDrop(const double dt);

	void extrapolateFieldsToAir(bool onlyv = false);
	void extrapolateFieldsToCollider();
	void enforceBoundaryConditions();

	void reinitializeLevelSet();
	void cutBoundaryLevelSet(GridBasedData<Dim, double> &sdf, GridBasedData<Dim, Vector<Dim, double>> &drv, const double tolerance = 0);

	double calculateVolume() const;
	double calculateVelocityDiv();

#ifndef _DISABLE_RM_
	void generateCurrentFields();
	bool checkReferenceMap(const double theta, const double bandStep, const double ratio) const;
	void restart();
#endif

#ifdef _ENABLE_PARTICLE_
	void advectParticle(double dt);
	void correctLevelSet();
    void resampleParticles();
	void reradiusParticles();
#endif

	bool checkDrop(const VectorDi &cell) const
	{
		if (_levelSet[cell] > 0) return false;

		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
			const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
			const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
			if (_boundary.fraction[axis][face] || _levelSet[nbCell] <= 0)
				return false;
		}
		return true;
	}
};

}
