#pragma once

#include "GlProgram.h"

namespace PhysX {

class GlRenderItem
{
protected:

	GLenum _mode = GL_TRIANGLES;
	GLint _first = 0;
	GLsizei _count = 0;
	GLenum _type = GL_UNSIGNED_INT;
	const void *_indices = 0;
	GLsizei _instanceCount = 1;
	GLint _baseVertex = 0;
	GLuint _baseInstance = 0;

	bool _indexed = false;
	bool _visible = true;

	GlProgram *const _program;

	GLuint _vao;

public:

	GlRenderItem(GlProgram *const program) : _program(program) { glCreateVertexArrays(1, &_vao); }

	GlRenderItem(const GlRenderItem &rhs) = delete;
	GlRenderItem &operator=(const GlRenderItem &rhs) = delete;

	virtual ~GlRenderItem() { glDeleteVertexArrays(1, &_vao); }

	virtual void beginDraw() { }
	virtual void endDraw() { }

	void draw();

	bool isVisible() const { return _visible; }
	void flipVisibility() { _visible = !_visible; }
};

}
