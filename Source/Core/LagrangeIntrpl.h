#pragma once

#include "GridBasedData.h"

#include <array>

namespace PhysX {

template <int Dim>
class LinearIntrpl
{
	DECLARE_DIM_TYPES(Dim)

	using WtPoint = std::pair<VectorDi, double>;

public:

	static std::array<VectorDi, Dim == 2 ? 4 : 8> points(const Grid<Dim> &grid, const VectorDd &pos)
	{
		const VectorDi lower = grid.getLinearLower(pos);
		if constexpr (Dim == 2) {
			return {
				lower + Vector2i(0, 0), lower + Vector2i(1, 0),
				lower + Vector2i(0, 1), lower + Vector2i(1, 1)
			};
		}
		else {
			return {
				lower + Vector3i(0, 0, 0), lower + Vector3i(1, 0, 0),
				lower + Vector3i(0, 1, 0), lower + Vector3i(1, 1, 0),
				lower + Vector3i(0, 0, 1), lower + Vector3i(1, 0, 1),
				lower + Vector3i(0, 1, 1), lower + Vector3i(1, 1, 1)
			};
		}
	}

	static auto wtPoints(const Grid<Dim> &grid, const VectorDd &pos)
	{
		const VectorDi lower = grid.getLinearLower(pos);
		const VectorDd frac = grid.getLowerFrac(pos, lower);
		const std::array<VectorDd, 2> w = {
			VectorDd::Ones() - frac,
			frac
		};

		if constexpr (Dim == 2) {
			return std::array {
				WtPoint(lower + Vector2i(0, 0), w[0][0] * w[0][1]),
				WtPoint(lower + Vector2i(1, 0), w[1][0] * w[0][1]),
				WtPoint(lower + Vector2i(0, 1), w[0][0] * w[1][1]),
				WtPoint(lower + Vector2i(1, 1), w[1][0] * w[1][1])
			};
		}
		else {
			return std::array {
				WtPoint(lower + Vector3i(0, 0, 0), w[0][0] * w[0][1] * w[0][2]),
				WtPoint(lower + Vector3i(1, 0, 0), w[1][0] * w[0][1] * w[0][2]),
				WtPoint(lower + Vector3i(0, 1, 0), w[0][0] * w[1][1] * w[0][2]),
				WtPoint(lower + Vector3i(1, 1, 0), w[1][0] * w[1][1] * w[0][2]),
				WtPoint(lower + Vector3i(0, 0, 1), w[0][0] * w[0][1] * w[1][2]),
				WtPoint(lower + Vector3i(1, 0, 1), w[1][0] * w[0][1] * w[1][2]),
				WtPoint(lower + Vector3i(0, 1, 1), w[0][0] * w[1][1] * w[1][2]),
				WtPoint(lower + Vector3i(1, 1, 1), w[1][0] * w[1][1] * w[1][2])
			};
		}
	}

	template <typename Type>
	static Type interpolate(const GridBasedData<Dim, Type> &gbd, const VectorDd &pos)
	{
		Type val = Zero<Type>();
		for (const auto [coord, weight] : wtPoints(gbd.grid, pos))
			val += gbd.at(coord) * weight;
		return val;
	}

	template <typename Type>
	static Type interpolate(const GridBasedData<Dim, Type> &gbd, const GridBasedData<Dim, Upgrade<Type, Dim>> &drv, const VectorDd &pos) { return interpolate(gbd, pos); }

	template <typename Type>
	static Upgrade<Type, Dim> interpolateGradient(const GridBasedData<Dim, Type> &gbd, const GridBasedData<Dim, Upgrade<Type, Dim>> &drv, const VectorDd &pos) { return interpolate(drv, pos); }
};

template <int Dim>
class QuadraticIntrpl
{
	DECLARE_DIM_TYPES(Dim)

	using WtPoint = std::pair<VectorDi, double>;

public:

	static auto wtPoints(const Grid<Dim> &grid, const VectorDd &pos)
	{
		const VectorDi lower = grid.getCubicLower(pos);
		const VectorDd frac = grid.getLowerFrac(pos, lower);
		const std::array<VectorDd, 3> w = {
			((frac.array() - 1) * (frac.array() - 2) / 2).matrix(),
			(frac.array() * (frac.array() - 2) / -1).matrix(),
			(frac.array() * (frac.array() - 1) / 2).matrix()
		};

		std::array<WtPoint, Dim == 2 ? 9 : 27> wtPts;
		if constexpr (Dim == 2) {
			for (int j = 0; j < 3; j++)
				for (int i = 0; i < 3; i++)
					wtPts[j * 3 + i] = WtPoint(lower + Vector2i(i, j), w[i][0] * w[j][1]);
		}
		else {
			for (int k = 0; k < 3; k++)
				for (int j = 0; j < 3; j++)
					for (int i = 0; i < 3; i++)
						wtPts[k * 9 + j * 3 + i] = WtPoint(lower + Vector3i(i, j, k), w[i][0] * w[j][1] * w[k][2]);
		}
		return wtPts;
	}

	template <typename Type>
	static Type interpolate(const GridBasedData<Dim, Type> &gbd, const VectorDd &pos)
	{
		Type val = Zero<Type>();
		for (const auto [coord, weight] : wtPoints(gbd.grid, pos))
			val += gbd.at(coord) * weight;
		return val;
	}
};

template <int Dim>
class CubicIntrpl
{
	DECLARE_DIM_TYPES(Dim)

	using WtPoint = std::pair<VectorDi, double>;

public:

	static auto wtPoints(const Grid<Dim> &grid, const VectorDd &pos)
	{
		const VectorDi lower = grid.getCubicLower(pos);
		const VectorDd frac = grid.getLowerFrac(pos, lower);
		const std::array<VectorDd, 4> w = {
			((frac.array() - 1) * (frac.array() - 2) * (frac.array() - 3) / -6).matrix(),
			(frac.array() * (frac.array() - 2) * (frac.array() - 3) / 2).matrix(),
			(frac.array() * (frac.array() - 1) * (frac.array() - 3) / -2).matrix(),
			(frac.array() * (frac.array() - 1) * (frac.array() - 2) / 6).matrix()
		};

		std::array<WtPoint, Dim == 2 ? 16 : 64> wtPts;
		if constexpr (Dim == 2) {
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					wtPts[j << 2 | i] = WtPoint(lower + Vector2i(i, j), w[i][0] * w[j][1]);
		}
		else {
			for (int k = 0; k < 4; k++)
				for (int j = 0; j < 4; j++)
					for (int i = 0; i < 4; i++)
						wtPts[k << 4 | j << 2 | i] = WtPoint(lower + Vector3i(i, j, k), w[i][0] * w[j][1] * w[k][2]);
		}
		return wtPts;
	}

	template <typename Type>
	static Type interpolate(const GridBasedData<Dim, Type> &gbd, const VectorDd &pos)
	{
		Type val = Zero<Type>();
		for (const auto [coord, weight] : wtPoints(gbd.grid, pos))
			val += gbd.at(coord) * weight;
		return val;
	}
};

}
