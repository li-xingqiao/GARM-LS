#pragma once

namespace PhysX {

struct BitmapFont final
{
	struct FontChar
	{
		int x;
		int y;
		int width;
		int height;
		int xOffset;
		int yOffset;
		int xAdvance;
	};

	struct FontCharN
	{
		float x;
		float y;
		float width;
		float texWidth;
		float height;
		float texHeight;
		float xOffset;
		float yOffset;
		float xAdvance;
	};

#include "BitmapConsolas.inc"

	int lineHeight;
	int base;
	int scaleW;
	int scaleH;
	int nChnls;
	int count;
	const unsigned char *data;
	const FontChar *chars;

	float lineHeightN;
	float baseN;
	FontCharN charsN[128];

	void normalize(const int clientWidth, const int clientHeight);
};

};
