#pragma once

#include "BitmapFont.h"
#include "GlRenderItem.h"

namespace PhysX {

class GlText : public GlRenderItem
{
public:

	enum class Alignment : uchar { Left, Right };

protected:

	static constexpr size_t _kTextBufferSize = 1024;

	struct Vertex
	{
		Vector4f pos;
		Vector4f texCoord;
		Vector4f color;
	};

	Vertex _vertices[_kTextBufferSize];

	GLuint _texture;
	GLuint _vbo;

	BitmapFont _bmConsolas = {
		64,
		51,
		512,
		512,
		1,
		128,
		reinterpret_cast<const uchar *>(BitmapFont::kBitmapConsolas[0]),
		BitmapFont::kBitmapConsolasChars
	};

public:

	GlText(GlProgram *const program);

	GlText(const GlText &rhs) = delete;
	GlText &operator=(const GlText &rhs) = delete;

	virtual ~GlText() { glDeleteBuffers(1, &_vbo); }

	void reset() { _instanceCount = 0; }
	void set(
		const std::string &text,
		const Vector2f &pos,
		const Vector2f &scale,
		const Vector4f &color,
		const Alignment alignment = Alignment::Left);
	void resize(const int width, const int height) { _bmConsolas.normalize(width, height); }

	virtual void beginDraw() override;
	virtual void endDraw() override;
};

}
