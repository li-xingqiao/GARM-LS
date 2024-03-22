#include "SdfToVDB.h"

#ifdef _WITH_OPENVDB

namespace PhysX {

void convertSdfToVDB(const GridBasedData<3, double> &sdf, openvdb::DoubleGrid &pgrid)
{
	const Vector3i size = sdf.grid.size;
	openvdb::Coord ijk;
	openvdb::DoubleGrid::Accessor accessor = pgrid.getAccessor();
	for (ijk[0] = 0; ijk[0] < size[0]; ++ijk[0])
		for (ijk[1] = 0; ijk[1] < size[1]; ++ijk[1])
			for (ijk[2] = 0; ijk[2] < size[2]; ++ijk[2]) {
				Vector3i coord(ijk[0], ijk[1], ijk[2]);
				accessor.setValue(ijk, sdf[coord]);
			}
}

void saveVDB(const GridBasedData<3, double> &sdf, const std::string &filename)
{
	openvdb::initialize();
	openvdb::DoubleGrid::Ptr grid = openvdb::DoubleGrid::create(1e5);
	convertSdfToVDB(sdf, *grid);

	openvdb::io::File file(filename);
	openvdb::GridPtrVec grids;
	grids.push_back(grid);
	file.write(grids);
	file.close();
}

} // namespace PhysX

#endif
