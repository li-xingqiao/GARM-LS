#include "GlText.h"

#include <iostream>

#include <cstdlib>

namespace PhysX {

GlText::GlText(GlProgram *const program) :
	GlRenderItem(program)
{
	glCreateTextures(GL_TEXTURE_2D, 1, &_texture);
	glTextureParameteri(_texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTextureParameteri(_texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTextureParameteri(_texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTextureParameteri(_texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTextureStorage2D(_texture, 5, GL_R8, _bmConsolas.scaleW, _bmConsolas.scaleH);
	glTextureSubImage2D(_texture, 0, 0, 0, _bmConsolas.scaleW, _bmConsolas.scaleH, GL_RED, GL_UNSIGNED_BYTE, _bmConsolas.data);
	glGenerateTextureMipmap(_texture);

	glCreateBuffers(1, &_vbo);
	glNamedBufferStorage(_vbo, sizeof(Vertex) * _kTextBufferSize, nullptr, GL_DYNAMIC_STORAGE_BIT);

	glVertexArrayVertexBuffer(_vao, 0, _vbo, 0, sizeof(Vertex));
	glVertexArrayAttribFormat(_vao, 0, 4, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribFormat(_vao, 1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float));
	glVertexArrayAttribFormat(_vao, 2, 4, GL_FLOAT, GL_FALSE, 8 * sizeof(float));
	glVertexArrayAttribBinding(_vao, 0, 0);
	glVertexArrayAttribBinding(_vao, 1, 0);
	glVertexArrayAttribBinding(_vao, 2, 0);
	glVertexArrayBindingDivisor(_vao, 0, 1);
	glEnableVertexArrayAttrib(_vao, 0);
	glEnableVertexArrayAttrib(_vao, 1);
	glEnableVertexArrayAttrib(_vao, 2);

	_mode = GL_TRIANGLE_STRIP;
	_count = 4;
}

void GlText::set(const std::string &text, const Vector2f &pos, const Vector2f &scale, const Vector4f &color, const Alignment alignment)
{
	const float headX = pos.x() * 2.0f - 1.0f;
	const float headY = (1.0f - pos.y()) * 2.0f - 1.0f;
	float x = headX;
	float y = headY;

	for (size_t lineBegin = 0, lineEnd; lineBegin < text.size(); lineBegin = lineEnd + 1) {
		lineEnd = lineBegin;
		while (lineEnd < text.size() && text[lineEnd] != '\n') lineEnd++;
		for (size_t i = 0; i < lineEnd - lineBegin; i++) {
			char c = text[alignment == Alignment::Left ? lineBegin + i : lineEnd - i - 1];
			const BitmapFont::FontCharN *fc = _bmConsolas.charsN + c;

			if (alignment == Alignment::Right)
				x -= fc->xAdvance * scale.x();

			_vertices[_instanceCount++] = {
				Vector4f(x + fc->xOffset * scale.x(), y - fc->yOffset * scale.y(), fc->width * scale.x(), fc->height * scale.y()),
				Vector4f(fc->x, fc->y, fc->texWidth, fc->texHeight),
				color
			};
			if (_instanceCount >= _kTextBufferSize) {
				std::cerr << "Error: [GlText] ran out of its text buffer." << std::endl;
				std::exit(-1);
			}

			if (alignment == Alignment::Left)
				x += fc->xAdvance * scale.x();
		}
		if (lineEnd < text.size()) {
			x = headX;
			y -= _bmConsolas.lineHeightN * scale.y();
		}
	}
}

void GlText::beginDraw()
{
	glNamedBufferSubData(_vbo, 0, sizeof(Vertex) * _instanceCount, _vertices);
	glBindTextureUnit(0, _texture);
}

void GlText::endDraw()
{
	glBindTextureUnit(0, 0);
}

}
