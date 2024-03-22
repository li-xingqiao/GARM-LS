#pragma once

#include "Types.h"

namespace PhysX {

class Camera
{
protected:

	float _aspect = 1.0f;
	float _yScale = 1.0f;
	float _xScale = 1.0f;

	Vector3f _pos = Vector3f::Zero();
	Matrix4f _projView = Matrix4f::Identity();

	bool _projDirty = true;

public:

	Camera() = default;
	virtual ~Camera() = default;

	virtual void update() = 0;
	virtual void reset() = 0;

	const Vector3f &pos() const { return _pos; }
	const Matrix4f &projView() const { return _projView; }

	void setAspectRatio(const float aspect)
	{
		_aspect = aspect;
		_projDirty = true;
	}

	virtual void rotate(const float dx, const float dy) { }
	virtual void translate(const float dx, const float dy) { }
	virtual void scale(const float dy) { }
};
}
