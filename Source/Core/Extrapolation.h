#pragma once

#include "ScalarField.h"
#include "SGridBasedData.h"

namespace PhysX::Extrapolation {

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const int, const Vector<Dim, int> &)>>
void solve(SGridBasedData<Dim, double> &sgbd, const ScalarField<Dim> &sdf, const double maxSteps, const Func isLiquidFace)
{
	DECLARE_DIM_TYPES(Dim)

	const double bandWidth = maxSteps * sgbd.grids[0].spacing;
	auto newSgbd = sgbd;

	parallelForEach(sgbd.grids, [&](const int axis, const VectorDi &face) {
		if (isLiquidFace(axis, face)) return;
		int cnt = 0;
		double sum = 0;
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi &nbFace = Grid<Dim>::neighbor(face, i);
			if (newSgbd[axis].grid.isValid(nbFace) && isLiquidFace(axis, nbFace))
				sum += sgbd[axis][nbFace], cnt++;
			if (cnt > 0) newSgbd[axis][face] = sum / cnt;
		}
	});

	parallelForEach(sgbd.grids, [&](const int axis, const VectorDi &face) {
		if (isLiquidFace(axis, face)) return;
		const VectorDd pos = sgbd.grids[axis].position(face);
		const double phi = sdf(pos);
		if (bandWidth < 0 || phi < bandWidth) {
			const VectorDd normal = sdf.gradient(pos).normalized();
			sgbd[axis][face] = LinearIntrpl<Dim>::interpolate(newSgbd[axis], pos - phi * normal);
		}
		else
		 	sgbd[axis][face] = 0; // for velocity only
	});
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const int, const Vector<Dim, int> &)>>
void solve(SGridBasedData<Dim, double> &sgbd, const double maxSteps, const Func isLiquidFace)
{
	DECLARE_DIM_TYPES(Dim)

	// extrapolate without levelset
	SGridBasedData<Dim, double> newsgbd = sgbd;
	SGridBasedData<Dim, uchar> visited0(sgbd.grids);
	SGridBasedData<Dim, uchar> visited(sgbd.grids);

	parallelForEach(sgbd.grids, [&](const int axis, const VectorDi &face) {
		if (isLiquidFace(axis, face)) {
			visited0[axis][face] = 1;
		}
		else {
			visited0[axis][face] = 0;
		}
	});

	for (int iter = 0; iter < maxSteps; iter++) {
		parallelForEach(sgbd.grids, [&](const int axis, const VectorDi &face) {
			if (!visited0[axis][face]) {
				int cnt = 0;
				double sum = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbFace = Grid<Dim>::neighbor(face, i);
					if (newsgbd.grids[axis].isValid(nbFace) && visited0[axis][nbFace])
						sum += newsgbd[axis][nbFace], cnt++;
				}
				if (cnt > 0) {
					newsgbd[axis][face] = sum / cnt;
					visited[axis][face] = true;
				} else {
					newsgbd[axis][face] = 0; // for velocity only
				}
			}
			else visited[axis][face] = true;
		});
		std::swap(visited0, visited);
	}
	sgbd = newsgbd;
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const Vector<Dim, int> &)>>
void solve(GridBasedData<Dim, double> &gbd, const double maxSteps, const Func isFixed)
{
	DECLARE_DIM_TYPES(Dim)

	// extrapolate without levelset
	GridBasedData<Dim, double> newgbd = gbd;
	GridBasedData<Dim, uchar> visited0(gbd.grid);
	GridBasedData<Dim, uchar> visited(gbd.grid);

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell) && gbd.grid.isValid(cell)) {
			visited0[cell] = 1;
		}
		else {
			visited0[cell] = 0;
		}
	});

	for (int iter = 0; iter < maxSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited0[cell]) {
				int cnt = 0;
				double sum = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited0[nbCell])
						sum += newgbd[nbCell], cnt++;
				}
				if (cnt > 0) {
					newgbd[cell] = sum / cnt;
					visited[cell] = true;
				}
			}
			else visited[cell] = true;
		});
		std::swap(visited0, visited);
	}
	gbd = newgbd;
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const Vector<Dim, int> &)>>
void solve(GridBasedData<Dim, uchar> &gbd, const double maxSteps, const Func isFixed)
{
	DECLARE_DIM_TYPES(Dim)

	// extrapolate without levelset
	GridBasedData<Dim, uchar> newgbd = gbd;
	GridBasedData<Dim, uchar> visited0(gbd.grid);
	GridBasedData<Dim, uchar> visited(gbd.grid);

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell) && gbd.grid.isValid(cell)) {
			visited0[cell] = 1;
		}
		else {
			visited0[cell] = 0;
		}
	});

	for (int iter = 0; iter < maxSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited0[cell]) {
				int cnt = 0;
				double sum = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited0[nbCell])
						sum += newgbd[nbCell], cnt++;
				}
				if (cnt > 0) {
					newgbd[cell] = 1;
					visited[cell] = true;
				}
			}
			else visited[cell] = true;
		});
		std::swap(visited0, visited);
	}
	
	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (!visited0[cell]) {
			newgbd[cell] = 0;
		}
	});
	gbd = newgbd;
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const Vector<Dim, int> &)>>
void solve(GridBasedData<Dim, Vector<Dim, double>> &gbd, const double maxSteps, const Func isFixed)
{
	DECLARE_DIM_TYPES(Dim)

	// extrapolate without levelset
	GridBasedData<Dim, Vector<Dim, double>> newgbd = gbd;
	GridBasedData<Dim, uchar> visited0(gbd.grid);
	GridBasedData<Dim, uchar> visited(gbd.grid);

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell) && gbd.grid.isValid(cell)) {
			visited0[cell] = 1;
		}
		else {
			visited0[cell] = 0;
		}
	});

	for (int iter = 0; iter < maxSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited0[cell]) {
				int cnt = 0;
				VectorDd sum = VectorDd::Zero();
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited0[nbCell])
						sum += newgbd[nbCell], cnt++;
				}
				if (cnt > 0) {
					newgbd[cell] = sum / cnt;
					visited[cell] = true;
				} else {
					newgbd[cell] = gbd.grid.position(cell); // for reference map only
				}
			}
			else visited[cell] = true;
		});
		std::swap(visited0, visited);
	}
	gbd = newgbd;
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const Vector<Dim, int> &)>>
void solve(GridBasedData<Dim, double> &gbd, const double maxSteps, const double fixSteps, const Func isFixed)
{
	DECLARE_DIM_TYPES(Dim)

	// extrapolate without levelset
	GridBasedData<Dim, double> newgbd = gbd;
	GridBasedData<Dim, uchar> visited0(gbd.grid);
	GridBasedData<Dim, uchar> visited(gbd.grid);

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell) && gbd.grid.isValid(cell)) {
			visited[cell] = 1;
		}
		else {
			visited[cell] = 0;
		}
	});

	for (int iter = 0; iter < fixSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited[cell]) {
				int cnt = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited[nbCell]) {
						visited0[cell] = true;
					}
				}
			}
			else visited0[cell] = true;
		});
		std::swap(visited0, visited);
	}

	for (int iter = 0; iter < maxSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited0[cell]) {
				int cnt = 0;
				double sum = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited0[nbCell])
						sum += newgbd[nbCell], cnt++;
				}
				if (cnt > 0) {
					newgbd[cell] = sum / cnt;
					visited[cell] = true;
				}
			}
			else visited[cell] = true;
		});
		std::swap(visited0, visited);
	}
	gbd = newgbd;
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const Vector<Dim, int> &)>>
void solve(GridBasedData<Dim, Vector<Dim, double>> &gbd, const double maxSteps, const double fixSteps, const Func isFixed)
{
	DECLARE_DIM_TYPES(Dim)

	// extrapolate without levelset
	GridBasedData<Dim, Vector<Dim, double>> newgbd = gbd;
	GridBasedData<Dim, uchar> visited0(gbd.grid);
	GridBasedData<Dim, uchar> visited(gbd.grid);

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell) && gbd.grid.isValid(cell)) {
			visited[cell] = 1;
		}
		else {
			visited[cell] = 0;
		}
	});

	for (int iter = 0; iter < fixSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited[cell]) {
				int cnt = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited[nbCell]) {
						visited0[cell] = true;
					}
				}
			}
			else visited0[cell] = true;
		});
		std::swap(visited0, visited);
	}

	for (int iter = 0; iter < maxSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited0[cell]) {
				int cnt = 0;
				VectorDd sum = VectorDd::Zero();
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited0[nbCell])
						sum += newgbd[nbCell], cnt++;
				}
				if (cnt > 0) {
					newgbd[cell] = sum / cnt;
					visited[cell] = true;
				}
			}
			else visited[cell] = true;
		});
		std::swap(visited0, visited);
	}
	gbd = newgbd;
}

template <int Dim, typename Field, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const Vector<Dim, int> &)>>
void solve(GridBasedData<Dim, Vector<Dim, double>> &gbd, GridBasedData<Dim, Matrix<Dim, double>> &drv, const Field &ref, const double maxSteps, const Func isFixed)
{
	DECLARE_DIM_TYPES(Dim)

	// extrapolate without levelset
	GridBasedData<Dim, Vector<Dim, double>> newgbd = gbd;
	GridBasedData<Dim, Matrix<Dim, double>> newdrv = drv;
	GridBasedData<Dim, uchar> visited0(gbd.grid);
	GridBasedData<Dim, uchar> visited(gbd.grid);

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell) && gbd.grid.isValid(cell)) {
			visited0[cell] = 1;
		}
		else {
			visited0[cell] = 0;
		}
	});

	for (int iter = 0; iter < maxSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited0[cell]) {
				int cnt = 0;
				VectorDd sum = VectorDd::Zero();
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited0[nbCell]) {
						cnt++;
						sum += newgbd[nbCell];
					}
				}
				if (cnt > 0) {
					newgbd[cell] = sum / cnt;
					newdrv[cell] = MatrixDd::Identity();
					visited[cell] = true;
				}
			}
			else visited[cell] = true;
		});
		std::swap(visited0, visited);
	}
	gbd = newgbd;
	drv = newdrv;
}

template <int Dim, typename Field, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const Vector<Dim, int> &)>>
void solve(GridBasedData<Dim, Vector<Dim, double>> &gbd, GridBasedData<Dim, Matrix<Dim, double>> &drv, const Field &ref, const double maxSteps, const double fixSteps, const Func isFixed)
{
	DECLARE_DIM_TYPES(Dim)

	// extrapolate without levelset
	GridBasedData<Dim, Vector<Dim, double>> newgbd = gbd;
	GridBasedData<Dim, Matrix<Dim, double>> newdrv = drv;
	GridBasedData<Dim, uchar> visited0(gbd.grid);
	GridBasedData<Dim, uchar> visited(gbd.grid);

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell) && gbd.grid.isValid(cell)) {
			visited[cell] = 1;
		}
		else {
			visited[cell] = 0;
		}
	});

	for (int iter = 0; iter < fixSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited[cell]) {
				int cnt = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited[nbCell]) {
						visited0[cell] = true;
					}
				}
			}
			else visited0[cell] = true;
		});
		std::swap(visited0, visited);
	}

	for (int iter = 0; iter < maxSteps; iter++) {
		parallelForEach(gbd.grid, [&](const VectorDi &cell) {
			if (!visited0[cell]) {
				int cnt = 0;
				VectorDd sum = VectorDd::Zero();
				VectorDd refsum = VectorDd::Zero();
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
					if (newgbd.grid.isValid(nbCell) && visited0[nbCell]) {
						cnt++;
						sum += newgbd[nbCell];
						refsum += newdrv[nbCell].transpose() * ref.gradient(newgbd[nbCell]);
					}
				}
				if (cnt > 0) {
					newgbd[cell] = sum / cnt;
					VectorDd refgrad = ref.gradient(newgbd[cell]);
					newdrv[cell] = refgrad * refsum.normalized().transpose() / refgrad.squaredNorm();
					visited[cell] = true;
				}
			}
			else visited[cell] = true;
		});
		std::swap(visited0, visited);
	}
	gbd = newgbd;
	drv = newdrv;
}

template <int Dim, typename Func>
	requires std::is_convertible_v<Func, std::function<bool(const Vector<Dim, int> &)>>
void solve(GridBasedData<Dim, Vector<Dim, double>> &gbd, const ScalarField<Dim> &sdf, const double maxSteps, const double fixSteps, const Func isFixed)
{
	DECLARE_DIM_TYPES(Dim)

	const double bandWidth = maxSteps * gbd.grid.spacing;
	const double fixWidth = fixSteps * gbd.grid.spacing;
	auto newgbd = gbd;

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell)) return;
		int cnt = 0;
		VectorDd sum = VectorDd::Zero();
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
			if (newgbd.grid.isValid(nbCell) && isFixed(nbCell))
				sum += newgbd[nbCell], cnt++;
			if (cnt > 0) newgbd[cell] = sum / cnt;
		}
	});

	parallelForEach(gbd.grid, [&](const VectorDi &cell) {
		if (isFixed(cell)) return;
		const VectorDd pos = gbd.grid.position(cell);
		const double phi = sdf(pos) - fixWidth;
		if (bandWidth < 0 || phi < bandWidth) {
			const VectorDd normal = sdf.gradient(pos).normalized();
			gbd[cell] = LinearIntrpl<Dim>::interpolate(newgbd, pos - phi * normal);
		}
	});
}

}
