#pragma once

#ifdef _WITH_OPENVDB

#include <openvdb/openvdb.h>

#include "GridBasedData.h"

namespace PhysX {

void convertSdfToVDB(const GridBasedData<3, double> &sdf, openvdb::DoubleGrid &pgrid);

void saveVDB(const GridBasedData<3, double> &sdf, const std::string &filename);

} // namespace PhysX

#endif
