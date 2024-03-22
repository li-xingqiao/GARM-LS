#pragma once

#include "Types.h"

#include <fstream>
#include <iostream>

namespace PhysX::IO {

template <typename Type>
inline void read(std::istream &in, Type *const data, const size_t cnt) { in.read(reinterpret_cast<char *>(data), cnt); }

template <typename Type>
inline void readValue(std::istream &in, Type &val) { read(in, &val, sizeof(Type)); }

template <typename Type>
inline void readArray(std::istream &in, Type *const data, const size_t cnt) { if (cnt > 0) read(in, data, sizeof(Type) * cnt); }

template <typename Type>
inline void write(std::ostream &out, const Type *const data, const size_t cnt) { out.write(reinterpret_cast<const char *>(data), cnt); }

template <typename Type>
inline void writeValue(std::ostream &out, const Type &val) { write(out, &val, sizeof(Type)); }

template <typename Type>
inline void writeArray(std::ostream &out, const Type *const data, const size_t cnt) { if (cnt > 0) write(out, data, sizeof(Type) * cnt); }

}
