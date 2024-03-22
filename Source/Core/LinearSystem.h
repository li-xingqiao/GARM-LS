#pragma once

#include "Types.h"

#include <format>
#include <iostream>

namespace PhysX {

template <typename MatrixType, typename Preconditioner = Eigen::IdentityPreconditioner>
using CG = Eigen::ConjugateGradient<MatrixType, Eigen::Lower | Eigen::Upper, Preconditioner>;

template <typename MatrixType>
using ICPCG = CG<MatrixType, Eigen::IncompleteCholesky<double>>;

template <typename MatrixType, typename Preconditioner = Eigen::IdentityPreconditioner>
using BiCGSTAB = Eigen::BiCGSTAB<MatrixType, Preconditioner>;

}

namespace PhysX::LinearSystem {

template <typename MatrixType, typename Solver = ICPCG<MatrixType>>
inline void solve(const MatrixType &A, Eigen::Ref<VectorXd, Eigen::Aligned> x, const Eigen::Ref<const VectorXd, Eigen::Aligned> &b, const int maxIterations = -1, const double tolerance = 1e-6)
{
	Solver solver(A);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: [IterativeSolver] failed to factorize matrix." << std::endl;
		std::exit(-1);
	}
	if (maxIterations >= 0) solver.setMaxIterations(maxIterations);
	solver.setTolerance(std::max(tolerance, std::numeric_limits<double>::epsilon() / b.norm()));
	x = solver.solveWithGuess(b, x);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: [IterativeSolver] failed to solve linear system." << std::endl;
		std::exit(-1);
	}
	std::cout << std::format("{:>6} iters", solver.iterations());
}

}
