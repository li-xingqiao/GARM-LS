#pragma once

#include "Yaml.h"

#include <fstream>
#include <string>

namespace PhysX {

class Simulation
{
protected:

	double _time = 0;

public:

	Simulation() = default;
	Simulation(const Simulation &rhs) = delete;
	Simulation &operator=(const Simulation &rhs) = delete;
	virtual ~Simulation() = default;

	void setTime(const double time) { _time = time; }
	virtual double getTimeStep(const uint frameRate, const double stepRate) const { return 1.0 / frameRate / stepRate; }

	virtual int dimension() const = 0;
	virtual void writeDescription(YAML::Node &root) const = 0;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const = 0;
	virtual void saveFrame(const std::string &frameDir) const = 0;
	virtual void loadFrame(const std::string &frameDir) = 0;

	virtual void initialize() { }
	virtual void advance(const double dt) = 0;
};

}
