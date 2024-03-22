#pragma once

#include "Constants.h"
#include "Types.h"

#include <vector>

namespace PhysX {

class GeometryGenerator
{
public:

	struct Vertex
	{
		Vector3f pos;
		Vector3f normal;
	};

	struct Data
	{
		std::vector<Vertex> vertices;
		std::vector<uint> indices32;
	};

public:

	static Data createBox(const Vector3f &lengths);
	static Data createBox(const float width, const float height, const float depth) { return createBox(Vector3f(width, height, depth)); }
	static Data createUVSphere(const float radius, const uint sliceCnt, const uint stackCnt);
	static Data createCircle(const float radius, const uint edgeCnt);
};

}
