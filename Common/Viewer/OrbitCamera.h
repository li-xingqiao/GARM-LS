#pragma once

#include "Camera.h"

#include <algorithm>
#include <numbers>

#include <cmath>

namespace PhysX {

class OrbitCamera : public Camera
{
protected:

	float _fovy;
	float _zNear;
	float _zFar;

	const float _savedRadius;
	const float _savedPhi;
	const float _savedTheta;
	const Vector3f _savedTarget;

	float _radius;
	float _phi;
	float _theta;
	Vector3f _target;

	Vector3f _front = Vector3f::Zero();
	Vector3f _up = Vector3f::Zero();
	Vector3f _right = Vector3f::Zero();

	Matrix4f _proj = Matrix4f::Identity();
	Matrix4f _view = Matrix4f::Identity();

	bool _viewDirty = true;

public:

	OrbitCamera(
		const float fovy,
		const float zNear,
		const float zFar,
		const float radius,
		const float phi,
		const float theta,
		const Vector3f &target)
		:
		_savedRadius(radius),
		_savedPhi(phi),
		_savedTheta(theta),
		_savedTarget(target)
	{
		setPerspective(fovy, zNear, zFar);
		setSpherical(radius, phi, theta, target);
	}

	virtual ~OrbitCamera() = default;

	virtual void update() override
	{
		if (_projDirty)
			updateProjMatrix();
		if (_viewDirty)
			updateViewMatrix();
		if (_projDirty || _viewDirty)
			_projView = _proj * _view;
		_projDirty = false;
		_viewDirty = false;
	}

	virtual void reset() override
	{
		setSpherical(_savedRadius, _savedPhi, _savedTheta, _savedTarget);
	}

	virtual void rotate(const float dx, const float dy) override
	{
		static constexpr float kRotateRatio = 0.25f * float(kPi) / 180.0f;
		_phi += kRotateRatio * dx;
		_theta += kRotateRatio * dy;
		_phi = std::fmod(_phi, 2.0f * float(kPi));
		if (_phi < 0.0f) _phi += 2.0f * float(kPi);
		_theta = std::clamp(_theta, 0.1f, float(kPi) - 0.1f);
		lookAt(sphericalToCartesian(_radius, _phi, _theta), _target);
	}

	virtual void translate(const float dx, const float dy) override
	{
		static constexpr float kTranslateRatio = 0.001f;
		_target += kTranslateRatio * _radius * (dy * _up - dx * _right);
		lookAt(sphericalToCartesian(_radius, _phi, _theta), _target);
	}

	virtual void scale(const float dy) override
	{
		static constexpr float kScaleRatio = 0.05f;
		_radius /= std::exp(kScaleRatio * dy);
		_radius = std::clamp(_radius, _zNear, _zFar);
		lookAt(sphericalToCartesian(_radius, _phi, _theta), _target);
	}

	static Vector3f sphericalToCartesian(const float radius, const float phi, const float theta)
	{
		return Vector3f(std::sin(theta) * std::sin(phi), std::cos(theta), std::sin(theta) * std::cos(phi)) * radius;
	}

protected:

	void updateProjMatrix()
	{
		_yScale = 1.0f / std::tan(_fovy / 2);
		_xScale = _yScale / _aspect;

		_proj << _xScale, 0, 0, 0,
			0, _yScale, 0, 0,
			0, 0, -(_zFar + _zNear) / (_zFar - _zNear), -2 * _zNear * _zFar / (_zFar - _zNear),
			0, 0, -1, 0;
	}

	void updateViewMatrix()
	{
		_view << _right.x(), _right.y(), _right.z(), -_right.dot(_pos),
			_up.x(), _up.y(), _up.z(), -_up.dot(_pos),
			-_front.x(), -_front.y(), -_front.z(), _front.dot(_pos),
			0, 0, 0, 1;
	}

	void setPerspective(const float fovy, const float zNear, const float zFar)
	{
		_fovy = fovy;
		_zNear = zNear;
		_zFar = zFar;
		_projDirty = true;
	}

	void setSpherical(const float radius, const float phi, const float theta, const Vector3f &target)
	{
		_radius = radius;
		_phi = phi;
		_theta = theta;
		lookAt(sphericalToCartesian(radius, phi, theta), target);
	}

	void lookAt(const Vector3f &pos, const Vector3f &target)
	{
		_pos = pos;
		_target = target;

		_front = (_target - _pos).normalized();
		_right = _front.cross(Vector3f::Unit(1)).normalized();
		_up = _right.cross(_front).normalized();

		_viewDirty = true;
	}
};

}
