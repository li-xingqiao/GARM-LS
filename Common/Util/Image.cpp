#include "Image.h"

#include "IO.h"

#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>

#include <cctype>
#include <cstring>

namespace PhysX {

void Image::load(const uint width, const uint height, const unsigned char *buffer, const bool transparent)
{
	_width = width, _height = height, _numChnls = 3 + transparent;
	_size = size_t(_numChnls) * _width * _height;
	_buffer = std::make_unique<unsigned char[]>(_size);
	std::memcpy(_buffer.get(), buffer, _size);
}

void Image::load(const uint width, const uint height, unsigned char *&&buffer, const bool transparent)
{
	_width = width, _height = height, _numChnls = 3 + transparent;
	_size = size_t(_numChnls) * _width * _height;
	_buffer.reset(buffer);
}

std::string Image::write(const std::string &fileName) const
{
	const std::string ext = std::filesystem::path(fileName).extension().string();
	if (toLower(ext) == ".bmp")
		return writeBMP(fileName);
	else if (ext == "")
		reportError("unspecified file format");
	else
		reportError("unsupported format " + ext);
	return "";
}

std::string Image::writeBMP(const std::string &fileName) const
{
	auto path = std::filesystem::path(fileName);
	if (toLower(path.extension().string()) != ".bmp") path += ".bmp";
	std::filesystem::create_directories(path.parent_path());
	std::ofstream fout(path.string(), std::ios::binary);
	uint bmpSize = _height * (_numChnls == 4 ? 4 * _width : (3 * _width + _width % 4));
	// Write bitmap file header.
	IO::writeValue(fout, ushort(19778)); // bfType
	IO::writeValue(fout, uint(54 + bmpSize)); // bfSize
	IO::writeValue(fout, ushort(0)); // bfReserved1
	IO::writeValue(fout, ushort(0)); // bfReserved2
	IO::writeValue(fout, uint(54)); // bfOffBits
	// Write bitmap info header.
	IO::writeValue(fout, uint(40)); // biSize
	IO::writeValue(fout, int(_width)); // biWidth
	IO::writeValue(fout, int(_height)); // biHeight
	IO::writeValue(fout, ushort(1)); // biPlanes
	IO::writeValue(fout, ushort(_numChnls * 8)); // biBitCount
	IO::writeValue(fout, uint(0)); // biCompression
	IO::writeValue(fout, uint(bmpSize)); // biSizeImage
	IO::writeValue(fout, int(0)); // biXPelsPerMeter
	IO::writeValue(fout, int(0)); // biYPelsPerMeter
	IO::writeValue(fout, uint(0)); // biClrUsed
	IO::writeValue(fout, uint(0)); // biClrImportant
	// Write bitmap pixels.
	for (uint y = 0; y < _height; y++) {
		for (uint x = 0; x < _width; x++) {
			IO::writeValue(fout, _buffer[index(x, y) * _numChnls + 2]);
			IO::writeValue(fout, _buffer[index(x, y) * _numChnls + 1]);
			IO::writeValue(fout, _buffer[index(x, y) * _numChnls + 0]);
			if (_numChnls == 4)
				IO::writeValue(fout, _buffer[index(x, y) * _numChnls + 3]);
		}
		if (_numChnls != 4)
			for (uint x = 0; x < _width % 4; x++)
				IO::writeValue(fout, uchar(0));
	}
	return path.lexically_normal().string();
}

std::string Image::toLower(std::string str)
{
	for (auto& ch : str) ch = std::tolower(ch);
	return str;
}

void Image::reportError(const std::string &msg) const
{
	std::cerr << std::format("Error: [Image] encountered {}.", msg) << std::endl;
	std::exit(-1);
}

}
