#pragma once

#include "GlRenderItem.h"

namespace PhysX {

class GlAxes : public GlRenderItem
{
protected:

	GLuint _vbo;

public:

	GlAxes(GlProgram *const program, const int dim, const float length) : GlRenderItem(program)
	{
		float vertices[] = {
			0.0f,	0.0f,	0.0f,	1.0f, 0.0f, 0.0f, 1.0f,
			length,	0.0f,	0.0f,	1.0f, 0.0f, 0.0f, 1.0f,
			0.0f,	0.0f,	0.0f,	0.0f, 1.0f, 0.0f, 1.0f,
			0.0f,	length,	0.0f,	0.0f, 1.0f, 0.0f, 1.0f,
			0.0f,	0.0f,	0.0f,	0.0f, 0.0f, 1.0f, 1.0f,
			0.0f,	0.0f,	length,	0.0f, 0.0f, 1.0f, 1.0f
		};
		glCreateBuffers(1, &_vbo);
		glNamedBufferStorage(_vbo, 14 * sizeof(float) * dim, vertices, 0);

		glVertexArrayVertexBuffer(_vao, 0, _vbo, 0, 7 * sizeof(float));
		glVertexArrayAttribBinding(_vao, 0, 0);
		glVertexArrayAttribBinding(_vao, 1, 0);
		glVertexArrayAttribFormat(_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
		glVertexArrayAttribFormat(_vao, 1, 4, GL_FLOAT, GL_FALSE, 3 * sizeof(float));
		glEnableVertexArrayAttrib(_vao, 0);
		glEnableVertexArrayAttrib(_vao, 1);

		_mode = GL_LINES;
		_count = 2 * dim;
	}

	GlAxes(const GlAxes &rhs) = delete;
	GlAxes &operator=(const GlAxes &rhs) = delete;

	virtual ~GlAxes() { glDeleteBuffers(1, &_vbo); }
};

}
