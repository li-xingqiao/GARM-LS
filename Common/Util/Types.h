#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <concepts>
#include <type_traits>

#include <cassert>
#include <cstddef>

namespace PhysX {

namespace detail {
template <typename Type, int Dim> struct UpgradeImpl { using DrvType = Eigen::Matrix<Type, Dim, 1>; };
template <typename Type, int Rows, int Dim> struct UpgradeImpl<Eigen::Matrix<Type, Rows, 1>, Dim> { using DrvType = Eigen::Matrix<Type, Rows, Dim>; };
template <typename Type> struct DowngradeImpl;
template <typename Type, int Dim> struct DowngradeImpl<Eigen::Matrix<Type, Dim, 1>> { using IntType = Type; };
template <typename Type, int Rows, int Dim> requires (Dim > 1) struct DowngradeImpl<Eigen::Matrix<Type, Rows, Dim>> { using IntType = Eigen::Matrix<Type, Rows, 1>; };
}
template <typename Type, int Dim> using Upgrade = typename detail::UpgradeImpl<Type, Dim>::DrvType;
template <typename Type> using Downgrade = typename detail::DowngradeImpl<Type>::IntType;

#define DECLARE_EIGEN_VECTOR_TYPES(type, t)							\
using	Array2##t			=	Eigen::Array2##t;					\
using	Array3##t			=	Eigen::Array3##t;					\
using	Array4##t			=	Eigen::Array4##t;					\
using	ArrayX##t			=	Eigen::ArrayX##t;					\
using	Vector2##t			=	Eigen::Vector2##t;					\
using	Vector3##t			=	Eigen::Vector3##t;					\
using	Vector4##t			=	Eigen::Vector4##t;					\
using	VectorX##t			=	Eigen::VectorX##t;

#define DECLARE_EIGEN_MATRIX_TYPES(type, t)							\
using	Matrix2##t			=	Eigen::Matrix2##t;					\
using	Matrix3##t			=	Eigen::Matrix3##t;					\
using	Matrix4##t			=	Eigen::Matrix4##t;					\
using	MatrixX##t			=	Eigen::MatrixX##t;					\
using	Triplet##t			=	Eigen::Triplet<type>;				\
using	SparseMatrix##t		=	Eigen::SparseMatrix<type>;

DECLARE_EIGEN_VECTOR_TYPES(int, i)
DECLARE_EIGEN_VECTOR_TYPES(float, f)
DECLARE_EIGEN_VECTOR_TYPES(double, d)

DECLARE_EIGEN_MATRIX_TYPES(float, f)
DECLARE_EIGEN_MATRIX_TYPES(double, d)

#undef DECLARE_EIGEN_VECTOR_TYPES
#undef DECLARE_EIGEN_MATRIX_TYPES

template <int Dim, typename Scalar> using Array = Eigen::Array<Scalar, Dim, 1, Eigen::ColMajor, Dim, 1>;
template <int Dim, typename Scalar> using Vector = Eigen::Matrix<Scalar, Dim, 1, Eigen::ColMajor, Dim, 1>;
template <int Dim, typename Scalar> using Matrix = Eigen::Matrix<Scalar, Dim, Dim, Eigen::ColMajor, Dim, Dim>;

#define DECLARE_DIM_TYPES(Dim)										\
static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");	\
using	ArrayDd				=	Array<Dim, double>;					\
using	ArrayDf				=	Array<Dim, float>;					\
using	ArrayDi				=	Array<Dim, int>;					\
using	VectorDd			=	Vector<Dim, double>;				\
using	VectorDf			=	Vector<Dim, float>;					\
using	VectorDi			=	Vector<Dim, int>;					\
using	MatrixDd			=	Matrix<Dim, double>;

using	uchar				=	unsigned char;
using	ushort				=	unsigned short;
using	uint				=	unsigned int;
using	llong				=	long long;
using	ullong				=	unsigned long long;

using	std::size_t;

template <typename Type>
inline Type Zero()
{
	if constexpr (requires { { Type::Zero() } -> std::convertible_to<Type>; })
		return Type::Zero();
	else
		return 0;
}

template <typename Type>
inline auto norm(const Type a)
{
	if constexpr (requires { a.squaredNorm(); })
		return a.squaredNorm();
	else if constexpr (requires { a.norm(); })
		return a.norm();
	else
		return std::abs(a);
}

template <typename Type>
inline Type cwiseMin(const Type a, const Type b)
{
	if constexpr (requires { { a.cwiseMin(b) } -> std::convertible_to<Type>; })
		return a.cwiseMin(b);
	else
		return a < b ? a : b;
}

template <typename Type>
inline Type cwiseMax(const Type a, const Type b)
{
	if constexpr (requires { { a.cwiseMax(b) } -> std::convertible_to<Type>; })
		return a.cwiseMax(b);
	else
		return a > b ? a : b;
}

template <typename Type>
inline bool inRange(const Type a, const Type l, const Type r)
{
	if constexpr (requires (Type a, Type b) { { (a.array() <= b.array()).all() } -> std::convertible_to<bool>; })
		return (l.array() <= a.array()).all() && (a.array() <= r.array()).all();
	else
		return l <= a && a <= r;
}

} // namespaxe PhysX
