#include "GlRenderItem.h"

namespace PhysX {

void GlRenderItem::draw()
{
	beginDraw();

	_program->use();
	glBindVertexArray(_vao);
	if (_indexed)
		glDrawElementsInstancedBaseVertexBaseInstance(_mode, _count, _type, _indices, _instanceCount, _baseVertex, _baseInstance);
	else
		glDrawArraysInstancedBaseInstance(_mode, _first, _count, _instanceCount, _baseInstance);
	glBindVertexArray(0);

	endDraw();
}

}
