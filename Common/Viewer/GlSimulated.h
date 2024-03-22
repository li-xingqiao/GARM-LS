#pragma once

#include "GlRenderItem.h"
#include "Yaml.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace PhysX {

class GlSimulated : public GlRenderItem
{
public:

	static inline std::unordered_map<std::string, GLenum> strToMode = {
		{ "point_list", GL_POINTS },
		{ "line_list", GL_LINES },
		{ "triangle_list", GL_TRIANGLES },
		{ "line_strip", GL_LINE_STRIP }
	};

protected:

	GLuint _vbo;
	GLuint _ebo;

	Vector4f _diffuseAlbedo;
	Vector3f _fresnelR0;
	float _roughness;
	bool _enableColorMap;

	uint _currentFrame = 0;

	std::vector<uint> _vtxFrameOffset;
	std::vector<uint> _idxFrameOffset;

public:

	GlSimulated(GlProgram *const program, const std::string &outputDir, const uint endFrame, const int dim, const YAML::Node &node);

	GlSimulated(const GlSimulated &rhs) = delete;
	GlSimulated &operator=(const GlSimulated &rhs) = delete;

	virtual ~GlSimulated()
	{
		glDeleteBuffers(1, &_ebo);
		glDeleteBuffers(1, &_vbo);
	}

	virtual void beginDraw() override;
	void setCurrentFrame(const uint frame) { _currentFrame = frame; }
	bool isTransparent() const { return _diffuseAlbedo.w() < 1.0f; }

protected:

	void loadMesh(
		const std::string &fileName,
		const int dim,
		std::vector<Vector3f> *const positions,
		std::vector<Vector3f> *const normals,
		std::vector<float> *const heats,
		std::vector<uint> *const indices);

	void reportError(const std::string &msg) const;
};

}
