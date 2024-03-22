#pragma once

#include "Camera.h"

#include <algorithm>

#include <cmath>

namespace PhysX {

class OrthoCamera : public Camera
{
protected:

	const Vector3f _savedPos;
	const float _savedYScale;
	const float _zScale;

public:

	OrthoCamera(const Vector3f &pos, const float yScale, const float zScale) :
		_savedPos(pos),
		_savedYScale(yScale),
		_zScale(zScale)
	{
		_pos = pos, _yScale = yScale;
		_projDirty = true;
	}

	virtual ~OrthoCamera() = default;

	void update()
	{
		if (_projDirty) {
			_xScale = _yScale / _aspect;
			_projView << _xScale, 0, 0, -_pos.x() * _xScale,
				0, _yScale, 0, -_pos.y() * _yScale,
				0, 0, -_zScale , _pos.z() *_zScale,
				0, 0, 0, 1;
		}
		_projDirty = false;
	}

	void reset()
	{
		_pos = _savedPos, _yScale = _savedYScale;
		_projDirty = true;
	}

	virtual void translate(const float dx, const float dy) override
	{
		static constexpr float kTranslateRatio = 0.0025f;
		_pos += kTranslateRatio * Vector3f(-dx, dy, 0) / _yScale;
		_projDirty = true;
	}

	void scale(const float dy)
	{
		static constexpr float kScaleRatio = 0.05f;
		_yScale *= std::exp(kScaleRatio * dy);
		_projDirty = true;
	}
};

}
