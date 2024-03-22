#pragma once

#include "GridBasedData.h"
#include "ScalarField.h"
#include "VectorField.h"

// #define ENABLE_NLOPT
// #define ENABLE_TRACING
#define ENABLE_MODIFIED
// #define ENABLE_LAGRANGE
#ifdef ENABLE_NLOPT
#include <nlopt.hpp>
#endif

namespace PhysX {

// #define _LOG_NEWTON

template <int Dim>
class HybridNewton
{
	DECLARE_DIM_TYPES(Dim)

public:

	class Velocity : public VectorField<Dim>
	{
	protected:

		const DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>> &_phiView;

	public:

		Velocity(const DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>> &phiView) : _phiView(phiView) { }

		virtual VectorDd operator()(const VectorDd &pos) const override
		{
			const VectorDd psi = _phiView.gradient(pos);
			return psi.any() ? psi.normalized() : VectorDd::Unit(0);
		}

		virtual MatrixDd gradient(const VectorDd &pos) const override { return MatrixDd::Zero(); }
	};

public:

	GridBasedData<Dim, uchar> _mark;
	GridBasedData<Dim, uchar> _diverged;
	GridBasedData<Dim, VectorDi> _across;

	HybridNewton(const Grid<Dim> &grid) : _mark(grid), _diverged(grid), _across(grid) { }

	void perform(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const int maxSteps, const double cfl);
	void perform(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const GridBasedData<Dim, VectorDd> &points, const int maxSteps, const double cfl);

	template <typename Field = DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>>>
	void perform(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const Field &phiView, const int maxSteps, const double cfl)
	{
		initializeInterface(phi, psi, phiView);
#ifdef _ENABLE_FMM_
		sweep(phi, psi, maxSteps);
#else
		advect(phi, psi, maxSteps, cfl);
#endif
	}

	template <typename Field = DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>>>
	void perform(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const GridBasedData<Dim, VectorDd> &points, const Field &phiView, const int maxSteps, const double cfl)
	{
		initializeInterface(phi, psi, points, phiView);
#ifdef _ENABLE_FMM_
		sweep(phi, psi, maxSteps);
#else
		advect(phi, psi, maxSteps, cfl);
#endif
	}

	void advect(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const int maxSteps, const double cfl) const;

	void sweep(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const int maxSteps);

protected:

	void prepareInterface(const GridBasedData<Dim, double> &phi);

	void initializeInterface(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi);
	void initializeInterface(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const GridBasedData<Dim, VectorDd> &points);
	void extrapolate(const int maxSteps);

	template <typename Field = DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>>>
	void initializeInterface(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const Field &phiView)
	{
		prepareInterface(phi);
		// Initialize interface cells by modified Newton's method.
		auto newPhi = phi;
		auto newPsi = psi;
		auto cposes = psi;
#ifdef _LOG_NEWTON
		forEach(phi.grid, [&](const VectorDi &coord) {
#else
		parallelForEach(phi.grid, [&](const VectorDi &coord) {
#endif
			if (_mark[coord]) {
				const VectorDd pos = phi.grid.position(coord);
				VectorDd closestPos = getClosestPosition(coord, pos, phiView);
				if (_diverged[coord]) {
					closestPos = getClosestPosition(coord, pos, closestPos, phiView);
				}
				newPhi[coord] = (pos - closestPos).norm() * (phi[coord] < 0 ? -1 : 1);
				newPsi[coord] = phiView.gradient(closestPos).normalized();
				cposes[coord] = closestPos;
			}
		});
// 		for (int iter = 0; iter < 5; ++iter) {
// #ifdef _LOG_NEWTON
// 			forEach(phi.grid, [&](const VectorDi &coord) {
// #else
// 			parallelForEach(phi.grid, [&](const VectorDi &coord) {
// #endif
// 				if (_diverged[coord]) {
// 					VectorDd sum = VectorDd::Zero();
// 					int cnt = 0;
// 					for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
// 						const VectorDi &nbCoord = Grid<Dim>::neighbor(coord, i);
// 						if (phi.grid.isValid(nbCoord) && _mark[nbCoord] && !_diverged[nbCoord])
// 							sum += cposes[nbCoord], cnt++;
// 					}
// 					if (cnt > 0) {
// 						const VectorDd pos = phi.grid.position(coord);
// 						const VectorDd closestPos = getClosestPosition(coord, pos, sum / cnt, phiView);
// 						newPhi[coord] = (pos - closestPos).norm() * (phi[coord] < 0 ? -1 : 1);
// 						newPsi[coord] = phiView.gradient(closestPos).normalized();
// 						cposes[coord] = closestPos;
// 					}
// 				}
// 			});
// 		}
		phi = newPhi, psi = newPsi;
	}

	template <typename Field = DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>>>
	void initializeInterface(GridBasedData<Dim, double> &phi, GridBasedData<Dim, VectorDd> &psi, const GridBasedData<Dim, VectorDd> &points, const Field &phiView)
	{
		prepareInterface(phi);
		// Initialize interface cells by modified Newton's method.
		auto newPhi = phi;
		auto newPsi = psi;
		auto cposes = psi;
#ifdef _LOG_NEWTON
		forEach(phi.grid, [&](const VectorDi &coord) {
#else
		parallelForEach(phi.grid, [&](const VectorDi &coord) {
#endif
			if (_mark[coord]) {
				const VectorDd pos = phi.grid.position(coord);
				VectorDd closestPos = getClosestPosition(coord, pos, points[coord], phiView);
				// if (_diverged[coord]) {
				// 	closestPos = points[coord];
				// }
				newPhi[coord] = (pos - closestPos).norm() * (phi[coord] < 0 ? -1 : 1);
				newPsi[coord] = phiView.gradient(closestPos).normalized();
				cposes[coord] = closestPos;
			}
		});
// 		for (int iter = 0; iter < 5; ++iter) {
// #ifdef _LOG_NEWTON
// 			forEach(phi.grid, [&](const VectorDi &coord) {
// #else
// 			parallelForEach(phi.grid, [&](const VectorDi &coord) {
// #endif
// 				if (_diverged[coord]) {
// 					VectorDd sum = VectorDd::Zero();
// 					int cnt = 0;
// 					for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
// 						const VectorDi &nbCoord = Grid<Dim>::neighbor(coord, i);
// 						if (phi.grid.isValid(nbCoord) && _mark[nbCoord] && !_diverged[nbCoord])
// 							sum += cposes[nbCoord], cnt++;
// 					}
// 					if (cnt > 0) {
// 						const VectorDd pos = phi.grid.position(coord);
// 						const VectorDd closestPos = getClosestPosition(coord, pos, sum / cnt, phiView);
// 						newPhi[coord] = (pos - closestPos).norm() * (phi[coord] < 0 ? -1 : 1);
// 						newPsi[coord] = phiView.gradient(closestPos).normalized();
// 						cposes[coord] = closestPos;
// 					}
// 				}
// 			});
// 		}
		phi = newPhi, psi = newPsi;
	}

	template <typename Field = DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>>>
	inline VectorDd getClosestPosition(const VectorDi &coord, const VectorDd &startPos, const Field &phiView, const double epsl = 0.0)
	{
		return getClosestPosition(coord, startPos, startPos, phiView, epsl);
	}

	template <typename Field = DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>>>
	VectorDd getClosestPosition(const VectorDi &coord, const VectorDd &startPos, const VectorDd &initPos, const Field &phiView, const double epsl = 0.0)
	{
		// using VectorAd = Vector<Dim + 1, double>;
		// using MatrixAd = Matrix<Dim + 1, double>;

		const double tolerance = 1e-3 * _mark.grid.spacing * _mark.grid.spacing * (Dim == 2 ? 1 : _mark.grid.spacing);
		const double epsphi = epsl * _mark.grid.spacing;
		bool converged = false;

#ifdef _LOG_NEWTON
		Eigen::IOFormat cleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]");
		std::cout.setf(std::ios::fixed);
		std::cout.precision(5);
		std::cout << "====== new pos =======" << std::endl;
		std::cout << "initial pos=" << (initPos / _mark.grid.spacing).transpose().format(cleanFmt) << std::endl;
		std::cout << "start   pos=" << (startPos / _mark.grid.spacing).transpose().format(cleanFmt) << std::endl;
#endif

#ifdef ENABLE_NLOPT

		struct Input {
			const VectorDd &qpos;
			const Field &view;
			Input(const VectorDd &p, const Field &v) : qpos(p), view(v) {}
		};
		Input state(startPos, phiView);

		nlopt::vfunc objective = [](const std::vector<double> &x, std::vector<double> &grad, void *f_data) -> double {
			Input *pinput = reinterpret_cast<Input *>(f_data);
			VectorDd pos, psi;
			double phi, res = 0.0;

			for (int i = 0; i < Dim; ++i) pos(i) = x[i];
			phi = pinput->view(pos);
			psi = pinput->view.gradient(pos).normalized();
			res += x[Dim] * phi;
			for (int i = 0; i < Dim; ++i) {
				grad[i] = x[i] - pinput->qpos(i) + x[Dim] * psi(i);
				res += 0.5 * (pinput->qpos(i) - x[i]) * (pinput->qpos(i) - x[i]);
			}
			return res;
		};
		
		nlopt::opt opt(nlopt::LD_LBFGS, Dim + 1);
		opt.set_min_objective(objective, &state);
		opt.set_xtol_rel(tolerance);
		opt.set_maxeval(20);

		VectorDd pos, psi = phiView.gradient(initPos);
		std::vector<double> var;
		for (int i = 0; i < Dim; ++i) var.push_back(initPos(i));
		var.push_back(-(initPos - startPos).dot(psi) / psi.squaredNorm());

		double f;
		nlopt::result result = opt.optimize(var, f);
		for (int i = 0; i < Dim; ++i) pos(i) = var[i];
		if (result != nlopt::SUCCESS) {
			converged = false;
		}
		return pos;

#else

		double radius = .9 * _mark.grid.spacing;
		double phi, norm;
		VectorDd pos, psi;
#if defined(ENABLE_LAGRANGE)
		double sigma = 1.0;
		double eta = 1.0;
		double lambda;
		VectorDd target, grad, delta, tdelta;
		psi = phiView.gradient(initPos);
		target = initPos;
		pos = initPos;
		phi = phiView(pos);
		lambda = psi.dot(startPos - pos) / psi.squaredNorm();
#elif defined(ENABLE_MODIFIED)
		double K = 1.0;
		int lsi;
		double pnorm;
		VectorDd delta, delta1, delta2, delta3, tdelta;
		VectorDd pdelta, ppos;
		pos = initPos;
		pdelta = VectorDd::Zero();
#elif defined(ENABLE_TRACING)
		double cphi;
		VectorDd cpos;
		VectorDd delta, dir;
		pos = initPos;
#else
		VectorDd grad, delta, tdelta;
		MatrixDd hessian;
		pos = initPos;
#endif

#ifdef ENABLE_LAGRANGE
		constexpr double RHO = 1.1; // Augmented Lagrangian Method
		for (int si = 0, iter = 0; si < 40; si += iter) {
			// solve for current sigma
			for (iter = 0; iter < 10; iter++) {
				psi = phiView.gradient(pos);

				// line search
				auto loss = [&phiView, &startPos, &sigma, &lambda](const VectorDd &x) -> double {
					double phii = phiView(x);
					return 0.5 * (startPos - x).squaredNorm() + lambda * phii + 0.5 * sigma * phii * phii;
				};

				grad = pos - startPos + lambda * psi + sigma * phi * psi;
				// for (int axis = 0; axis < Dim; ++axis) {
				// 	grad(axis) = (loss(pos + phi * VectorDd::Unit(axis)) - loss(pos - phi * VectorDd::Unit(axis))) / (2 * phi);
				// }
				delta = -grad;

				constexpr double GAMMA = 0.5;
				tdelta = delta * GAMMA;
				int lsi = 0;
				while (loss(pos + tdelta) < loss(pos + delta) && lsi < 5) {
					delta = tdelta;
					tdelta *= GAMMA;
					++lsi;
				}
				if (loss(pos + delta) > loss(pos)) break;

#ifdef _LOG_NEWTON
				std::cout << "current[" << si + 1 << "]: phi=" << (phi / _mark.grid.spacing) << " \tsigma=" << sigma << " \tloss=" << loss(pos) << " \tlambda=" << lambda << " \tpos=" << (pos / _mark.grid.spacing).transpose().format(cleanFmt) << " \tpsi=" << psi.transpose().format(cleanFmt) << " \tgrad=" << (grad / _mark.grid.spacing).transpose().format(cleanFmt) << " \tdelta=" << (delta / _mark.grid.spacing).transpose().format(cleanFmt) << std::endl;
#endif

				pos += delta;
				phi = phiView(pos);
				if (grad.norm() < tolerance) {
					converged = true;
					break;
				}
			}
			if (converged && std::abs(phi) < tolerance) {
				break;
			} else {
				converged = false;
			}
			target = pos;
			lambda += phi; // update lambda
			sigma = sigma * RHO;   // update sigma
		}
#elif defined(ENABLE_TRACING)
#ifdef _LOG_NEWTON
		auto intersect = [&phiView, &tolerance, &cleanFmt, this] (const VectorDd &pos, const VectorDd &dir) -> VectorDd {
#else
		auto intersect = [&phiView, &tolerance, this] (const VectorDd &pos, const VectorDd &dir) -> VectorDd {
#endif
			// requires dir to be normalized
			constexpr double RHOI = 0.9;
			double dp;
			VectorDd res = pos;
			#ifdef _LOG_NEWTON
				std::cout << "intersect: pos=" << (pos / _mark.grid.spacing).transpose().format(cleanFmt) << " \tdir=" << dir.transpose().format(cleanFmt);
				std::cout << std::endl;
			#endif
			for (int i = 0; i < 10; ++i) {
				double val = phiView(res);
				double grad = phiView.gradient(res).dot(dir);
				// line search
				dp = 20 * tolerance * (grad * val > 0 ? -1 : 1);
				while (std::abs(phiView(res + dp * RHOI * dir)) > std::abs(val) - std::abs(0.1 * grad * dp)) // Armijo
					dp *= RHOI;
#ifdef _LOG_NEWTON
				std::cout << "current: phi=" << (val / _mark.grid.spacing)  << " \tpos=" << (res / _mark.grid.spacing).transpose().format(cleanFmt) << " \tpsi=" << phiView.gradient(res).transpose().format(cleanFmt) << " \tdist=" << (res.dot(dir) / _mark.grid.spacing) << " \tgrad=" << grad << " \tdp=" << (dp / _mark.grid.spacing) << " \tcostheta=" << phiView.gradient(res).normalized().dot(dir) << std::endl;
#endif
				res += dp * dir;
				if (std::abs(phiView(res)) < tolerance || std::abs(dp) < 0.1 * tolerance) {
#ifdef _LOG_NEWTON
				std::cout << "current: phi=" << (phiView(res) / _mark.grid.spacing)  << " \tpos=" << (res / _mark.grid.spacing).transpose().format(cleanFmt) << " \tpsi=" << phiView.gradient(res).transpose().format(cleanFmt) << " \tdist=" << (res.dot(dir) / _mark.grid.spacing) << " \tgrad=" << grad << " \tdp=" << (dp / _mark.grid.spacing) << " \tcostheta=" << phiView.gradient(res).normalized().dot(dir) << std::endl;
#endif
					return res;
				}
			}
			return res;
		};

		phi = phiView(pos);
		psi = phiView.gradient(pos);
		if ((initPos - startPos).norm() < 1e-3) {
			delta = -phi * psi / psi.squaredNorm();
			dir = phiView.gradient(pos + delta).normalized();
		} else {
			dir = (initPos - startPos).normalized();
		}
		cpos = intersect(startPos + (pos - startPos).dot(dir) * dir, dir);
		if (std::abs(phiView(cpos)) < tolerance && (cpos - startPos).norm() < cphi) {
			pos = cpos;
			cphi = (cpos - startPos).norm();
		} else {
			dir = (_mark.grid.position(_across[coord]) - startPos).normalized();
		}

		cphi = 2e7;
		for (int iter = 0; iter < 10; ++iter) {
			cpos = intersect(startPos + (cpos - startPos).dot(dir) * dir, dir);
			if (std::abs(phiView(cpos)) < tolerance && (cpos - startPos).norm() < cphi) {
				pos = cpos;
				cphi = (cpos - startPos).norm();
				phi = phiView(cpos);
				psi = phiView.gradient(cpos).normalized();
				if ((dir - psi).norm() < 1e-3) break;
				dir = (5 * dir - psi).normalized();
			} else {
#ifdef _LOG_NEWTON
				std::cout << "no intersection!" << std::endl;
#endif
				phi = phiView(cpos);
				psi = phiView.gradient(cpos).normalized();
				dir = (cpos - phi * psi / psi.squaredNorm() - startPos).normalized();
			}
		}
		if (std::abs(phiView(pos)) < tolerance) {
			converged = true;
		} else {
			pos = cpos;
		}

#else
		// for (;;) {
		for (int iter = 0; iter < 20; iter++) {
#ifdef ENABLE_MODIFIED
			phi = phiView(pos) - epsphi;
			psi = phiView.gradient(pos);

			delta1 = -phi * psi / psi.squaredNorm();
			delta2 = ((startPos - pos) - psi * (startPos - pos).dot(psi) / psi.squaredNorm());
			// if (delta2.norm() > K * delta1.norm())
			// 	delta2 = K * delta1.norm() * delta2.normalized();
			delta = delta1 + delta2;
			norm = delta.norm();
			if (norm > radius) {
				delta *= radius / norm;
				norm = radius;
			}
#ifdef _ENABLE_FMM_
			// go back and forth
			pnorm = (delta + pdelta).norm();
			if (pnorm < 0.1 * norm) {
				pdelta = (pos + ppos) / 2;
				ppos = pos;
				pos = pdelta;
				pdelta = VectorDd::Zero();
			}
			else {
				ppos = pos;
				pdelta = delta;
				pos += delta;
			}
#else
			pos += delta;
#endif

			// // line search
			// constexpr double GAMMA = 0.9;
			// delta = VectorDd::Zero();
			// delta1 = -phi * psi / psi.squaredNorm();
			// tdelta = delta1 * GAMMA;
			// lsi = 0;
			// while (std::abs(phiView(pos + delta + tdelta)) < std::abs(phiView(pos + delta + delta1)) && lsi < 5) {
			// 	delta1 = tdelta;
			// 	tdelta *= GAMMA;
			// 	++lsi;
			// }
			// delta += delta1;
			// delta2 = ((startPos - pos) - psi * (startPos - pos).dot(psi) / psi.squaredNorm());
			// auto loss = [&phiView, &startPos](const VectorDd &x) -> double {
			// 	VectorDd psii = phiView.gradient(x);
			// 	if constexpr (Dim == 3) {
			// 		return ((startPos - x).cross(psii)).norm();
			// 	} else {
			// 		return std::abs((startPos - x)(0) * psii(1) - (startPos - x)(1) * psii(0));
			// 	}
			// };
			// tdelta = delta2 * GAMMA;
			// lsi = 0;
			// while (loss(pos + delta + tdelta) < loss(pos + delta + delta2) && lsi < 5) {
			// 	delta2 = tdelta;
			// 	tdelta *= GAMMA;
			// 	++lsi;
			// }
			// delta += delta2;
			// if constexpr (Dim == 3) {
			// 	delta3 = delta1.cross(delta2).normalized();
			// 	delta3 *= (startPos - pos).norm() * delta3.dot(phiView.gradient(pos + delta).normalized());
			// 	tdelta = delta3 * GAMMA;
			// 	lsi = 0;
			// 	while (loss(pos + delta + tdelta) + std::abs(phiView(pos + delta + tdelta)) < loss(pos + delta + delta3) + std::abs(phiView(pos + delta + tdelta)) && lsi < 5) {
			// 		delta3 = tdelta;
			// 		tdelta *= GAMMA;
			// 		++lsi;
			// 	}
			// 	delta += delta3;
			// }
			// norm = delta.norm();
			// ppos = pos;
			// pos += delta;
#ifdef _LOG_NEWTON
			std::cout << "current: phi=" << (phi / _mark.grid.spacing) << " \tpos=" << ((pos - delta) / _mark.grid.spacing).transpose().format(cleanFmt) << " \tpsi=" << psi.transpose().format(cleanFmt) << " \tdelta1=" << (delta1 / _mark.grid.spacing).transpose().format(cleanFmt) << " \tdelta2=" << (delta2 / _mark.grid.spacing).transpose().format(cleanFmt);
			// std::cout << "current: phi=" << (phi / _mark.grid.spacing) << " \tloss="  << loss(ppos) << " \tpos=" << (ppos / _mark.grid.spacing).transpose().format(cleanFmt) << " \tpsi=" << psi.transpose().format(cleanFmt) << " \tdelta1=" << (delta1 / _mark.grid.spacing).transpose().format(cleanFmt) << " \tdelta2=" << (delta2 / _mark.grid.spacing).transpose().format(cleanFmt);
			// if constexpr (Dim == 3) {
			// 	std::cout << " \tdelta3=" << (delta3 / _mark.grid.spacing).transpose().format(cleanFmt);
			// }
			std::cout << std::endl;
#endif
			pos = pos - epsphi * phiView.gradient(pos);
			// radius *= double(iter + 9) / (iter + 10);
			if (norm < tolerance) {
				converged = true; break;
			}
#else
			constexpr double RHO = 2.0; // Penalty-Based Method

			phi = phiView(pos);
			psi = phiView.gradient(pos);

			grad = pos - startPos + RHO * phi * psi;

			hessian = MatrixDd::Identity() + RHO * psi * psi.transpose(); // + lambda * hessian of phi. evaluated as zero

			// delta = -hessian.householderQr().solve(grad);
			delta = -grad;
			norm = delta.norm();

			// line search
			auto loss = [&phiView, &startPos](const VectorDd &x) -> double {
				double phii = phiView(x);
				return 0.5 * (startPos - x).squaredNorm() + 0.5 * RHO * phii * phii;
			};

			constexpr double GAMMA = 0.9;
			decltype(delta) tdelta = GAMMA * delta;
			while (loss(pos + tdelta) < loss(pos + delta) && delta.norm() > 0.5 * tolerance)
				delta = tdelta, tdelta *= GAMMA;

#ifdef _LOG_NEWTON
			std::cout << "current: phi=" << (phi / _mark.grid.spacing) << " \tloss=" << loss(pos) << " \tpos=" << (pos / _mark.grid.spacing).transpose().format(cleanFmt) << " \tpsi=" << psi.transpose().format(cleanFmt) << " \tgrad=" << (grad / _mark.grid.spacing).transpose().format(cleanFmt) << " \tdelta=" << (delta / _mark.grid.spacing).transpose().format(cleanFmt) << std::endl;
#endif

			// if (norm > radius) {
			// 	delta = (radius / norm) * delta;
			// }
			pos += delta;
			if (norm < tolerance) {
				converged = true; break;
			}
#endif
		}
#endif

#endif

		if (!converged) {
#ifdef _LOG_NEWTON
			std::cout << "diverged! current: phi=" << (phi / _mark.grid.spacing) << " pos=(" << (pos / _mark.grid.spacing).transpose().format(cleanFmt) << ") psi=(" << psi.transpose().format(cleanFmt) << ") delta=(" << (delta / _mark.grid.spacing).transpose().format(cleanFmt) << ") " << std::endl;
#endif
// 			// // find from nearest
// 			VectorDd min, max, mid;
// 			min = _mark.grid.position(_across[coord]);
// 			max = startPos;
// 			mid = (min + max) / 2;
// 			while ((max - min).norm() > tolerance) {
// 				if (phiView(mid) * phiView(min) >= 0) {
// 					min = mid;
// 				} else {
// 					max = mid;
// 				}
// 				mid = (min + max) / 2;
// #ifdef _LOG_NEWTON
// 			std::cout << "binary search: phi=" << (phiView(mid) / _mark.grid.spacing) << " pos=(" << (mid / _mark.grid.spacing).transpose().format(cleanFmt) << ")" << std::endl;
// #endif
// 			}
// 			pos = mid;
			_diverged[coord] = true;
			// _mark[coord] = false;
		}
		return pos;
	}

	// template <typename Field = DrvScalarFieldWrapper<Dim, GridBasedData<Dim, double>, HermiteIntrpl<Dim, double>>>
	// VectorDd getClosestPosition(const VectorDi &coord, const VectorDd &startPos, const VectorDd &initPos, const Field &phiView, const double epsl = 0.0)
	// {
	// 	const double tolerance = 1e-5 * _mark.grid.spacing * _mark.grid.spacing * (Dim == 2 ? 1 : _mark.grid.spacing);
	// 	double radius = .9 * _mark.grid.spacing;
	// 	VectorDd pos = initPos;
	// 	bool converged = false;
	// 	double phi;
	// 	VectorDd psi, delta;
	// 	// for (;;) {
	// 	for (int iter = 0; iter < 20; iter++) {
	// 		phi = phiView(pos);
	// 		psi = phiView.gradient(pos);
	// 		delta = -phi * psi / psi.squaredNorm() + (startPos - pos) - psi * (startPos - pos).dot(psi) / psi.squaredNorm();
	// 		double norm = delta.norm();
	// 		if (norm > radius) {
	// 			delta *= radius / norm;
	// 			norm = radius;
	// 		}
	// 		// std::cout << "current: phi=" << phi << " pos=" << pos.transpose() << " psi=" << psi.transpose() << " delta=" << delta.transpose() << std::endl;
	// 		pos += delta;
	// 		if (norm < tolerance) {
	// 			converged = true; break;
	// 		}
	// 	}
	// 	if (!converged) {
	// 			// std::cout << "current: phi=" << phi << " psi=(" << psi.transpose() << ") pos=(" << pos.transpose() << ") delta=(" << delta.transpose() << ") " << delta.norm() << " " << radius << std::endl;
	// 		_mark[coord] = false;
	// 	}
	// 	return pos;
	// }

};

}
