#include "BitmapFont.h"

#include "Types.h"

namespace PhysX {

void BitmapFont::normalize(const int clientWidth, const int clientHeight)
{
	lineHeightN = float(lineHeight) / clientHeight;
	baseN = float(base) / clientHeight;

	for (int i = 0; i < count; i++) {
		charsN[i] = {
			float(chars[i].x) / scaleW,
			float(chars[i].y) / scaleH,
			float(chars[i].width) / clientWidth,
			float(chars[i].width) / scaleW,
			float(chars[i].height) / clientHeight,
			float(chars[i].height) / scaleH,
			float(chars[i].xOffset) / clientWidth,
			float(chars[i].yOffset) / clientHeight,
			float(chars[i].xAdvance) / clientWidth
		};
	}
}

}
