#include "GeometryGenerator.h"

#include <iostream>
#include <numbers>

#include <cmath>
#include <cstdlib>

namespace PhysX {

GeometryGenerator::Data GeometryGenerator::createBox(const Vector3f &lengths)
{
	Data data;
	Vector3f halfLengths = lengths * 0.5f;

	Vector3i coeff;
	for (int axis = 0; axis < 3; axis++) {
		for (int dir = 0; dir < 2; dir++) {
			uint base = uint(data.vertices.size());
			coeff[axis] = dir * 2 - 1;
			for (int i = 0; i < 4; i++) {
				coeff[(axis + 1) % 3] = (i & 1) * 2 - 1;
				coeff[(axis + 2) % 3] = (i & 2) - 1;
				Vector3f pos = halfLengths.cwiseProduct(coeff.template cast<float>());
				Vector3f normal = Vector3f::Unit(axis) * coeff[axis];
				data.vertices.push_back({ pos, normal });
			}
			if (dir) {
				data.indices32.insert(data.indices32.end(), { base + 0, base + 1, base + 3 });
				data.indices32.insert(data.indices32.end(), { base + 3, base + 2, base + 0 });
			}
			else {
				data.indices32.insert(data.indices32.end(), { base + 3, base + 1, base + 0 });
				data.indices32.insert(data.indices32.end(), { base + 0, base + 2, base + 3 });
			}
		}
	}
	return data;
}

GeometryGenerator::Data GeometryGenerator::createUVSphere(const float radius, const uint sliceCnt, const uint stackCnt)
{
	assert(radius > 0.0f && sliceCnt >= 3 && stackCnt >= 2);

	Data data;
	const float deltaPhi = float(kPi) / stackCnt; // latitude, in [-pi / 2, +pi / 2]
	const float deltaLambda = float(kPi) * 2.0f / sliceCnt; // longitude, in [-pi, pi]

	data.vertices.push_back({ Vector3f::Unit(1) * (-radius), Vector3f::Unit(1) * (-1) });
	for (uint i = 1; i <= stackCnt - 1; i++) {
		const float phi = i * deltaPhi - float(kPi) * 0.5f;
		for (uint j = 0; j < sliceCnt; j++) {
			const float lambda = j * deltaLambda - float(kPi);
			Vertex vert;
			vert.normal = Vector3f(std::cos(phi) * std::sin(lambda), std::sin(phi), std::cos(phi) * std::cos(lambda));
			vert.pos = vert.normal * radius;
			data.vertices.push_back(vert);
		}
	}
	data.vertices.push_back({ Vector3f::Unit(1) * radius, Vector3f::Unit(1) });

	for (uint i = 1; i <= sliceCnt; i++)
		data.indices32.insert(data.indices32.end(), { 0, 1 + i % sliceCnt, i });
	for (uint i = 0; i < stackCnt - 2; i++)
		for (uint j = 0; j < sliceCnt; j++) {
			const uint index0 = 1 + i * sliceCnt + j;
			const uint index1 = 1 + i * sliceCnt + (j + 1) % sliceCnt;
			const uint index2 = 1 + (i + 1) * sliceCnt + j;
			const uint index3 = 1 + (i + 1) * sliceCnt + (j + 1) % sliceCnt;
			data.indices32.insert(data.indices32.end(), { index2, index0, index1 });
			data.indices32.insert(data.indices32.end(), { index1, index3, index2 });
		}
	const uint topIndex = uint(data.vertices.size()) - 1;
	for (uint i = 1; i <= sliceCnt; i++)
		data.indices32.insert(data.indices32.end(), { topIndex, topIndex - i % sliceCnt - 1, topIndex - i });

	return data;
}

GeometryGenerator::Data GeometryGenerator::createCircle(const float radius, const uint edgeCnt)
{
	Data data;
	const float deltaAlpha = float(kPi) * 2.0f / edgeCnt;
	for (uint i = 0; i < edgeCnt; i++) {
		const float alpha = i * deltaAlpha;
		Vertex vert;
		vert.pos = Vector3f(std::cos(alpha), std::sin(alpha), 0) * radius;
		vert.normal = Vector3f::Unit(2);
		data.vertices.push_back(vert);
	}
	for (uint i = 1; i < edgeCnt; i++) {
		data.indices32.insert(data.indices32.end(), { 0, i, (i + 1) % edgeCnt });
	}
	return data;
}

}
