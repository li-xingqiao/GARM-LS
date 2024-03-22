#pragma once

#include <chrono>

namespace PhysX {

class StepTimer final
{
public:

	using duration = std::chrono::duration<double>;
	using steady_clock = std::chrono::steady_clock;
	using time_point = std::chrono::time_point<std::chrono::steady_clock, duration>;

protected:

	duration _deltaTime = duration::zero();
	duration _pausedTime = duration::zero();

	time_point _baseTime;
	time_point _stopTime;
	time_point _prevTime;
	time_point _currTime;

	bool _stopped = false;

public:

	StepTimer() = default;
	virtual ~StepTimer() = default;

	auto totalTime() const { return (_stopped ? _stopTime : _currTime) - _baseTime - _pausedTime; }
	auto deltaTime() const { return _deltaTime; }

	void reset()
	{
		_baseTime = steady_clock::now();
		_prevTime = _baseTime;
		_pausedTime = duration::zero();
		_stopped = false;
	}

	void start()
	{
		if (_stopped) {
			_prevTime = steady_clock::now();
			_pausedTime += _prevTime - _stopTime;
			_stopped = false;
		}
	}

	void stop()
	{
		if (!_stopped) {
			_stopTime = steady_clock::now();
			_stopped = true;
		}
	}

	void tick()
	{
		if (_stopped) {
			_deltaTime = duration::zero();
			return;
		}
		_currTime = steady_clock::now();
		_deltaTime = _currTime - _prevTime;
		_prevTime = _currTime;
		if (_deltaTime < duration::zero())
			_deltaTime = duration::zero();
	}
};

}
