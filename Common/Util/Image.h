#pragma once

#include "Types.h"

#include <memory>
#include <string>

namespace PhysX {

class Image
{
protected:

	uint _width = 0;
	uint _height = 0;
	uint _numChnls = 0;
	size_t _size = 0;

	std::unique_ptr<unsigned char[]> _buffer;

public:

	Image() = default;
	Image(const uint width, const uint height, const unsigned char *buffer, const bool transparent = false) { load(width, height, buffer, transparent); }
	Image(const uint width, const uint height, unsigned char *&&buffer, const bool transparent = false) { load(width, height, buffer, transparent); }

	uint width() const { return _width; }
	uint height() const { return _height; }
	size_t size() const { return _size; }
	bool hasAlpha() const { return _numChnls == 4; }

	size_t index(const uint x, const uint y) const { return size_t(y) * _width + x; }

	void load(const uint width, const uint height, const unsigned char *buffer, const bool transparent = false);
	void load(const uint width, const uint height, unsigned char *&&buffer, const bool transparent = false);

	std::string write(const std::string &fileName) const;
	std::string writeBMP(const std::string &fileName) const;

	static std::string toLower(std::string str);

protected:

	void reportError(const std::string &msg) const;
};

}
