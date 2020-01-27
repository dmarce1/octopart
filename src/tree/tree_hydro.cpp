#include <octopart/math.hpp>
#include <octopart/tree.hpp>
#include <octopart/options.hpp>
#include <octopart/profiler.hpp>


void tree::compute_gradients(fixed_real t) {
	static auto opts = options::get();
	if (leaf) {
		PROFILE();
		if (!opts.first_order_space) {
			for (int i = 0; i < nparts0; i++) {
				auto &pi = parts[i];
				if (pi.t == t || opts.global_time) {
					primitive_state max_ngb;
					primitive_state min_ngb;
					const auto piW = pi.W;
					max_ngb = min_ngb = piW;
					for (int dim = 0; dim < NDIM; dim++) {
						pi.dW[dim].zero();
					}
					for (const auto &pj : parts) {
						const auto r = abs(pi.x - pj.x);
						const auto &h = pi.h;
						if (r < h) {
							const auto pjW = pj.W;
							for (int dim = 0; dim < NDIM; dim++) {
								real psi_a_j = 0.0;
								for (int m = 0; m < NDIM; m++) {
									psi_a_j += pi.B[dim][m] * (pj.x[m] - pi.x[m]) * W(r, h) * pi.V;
								}
								pi.dW[dim] = pi.dW[dim] + (pjW - piW) * psi_a_j;
							}
						}
					}
					real max_dx = 0.0;
					for (const auto &pj : parts) {
						const auto r = abs(pi.x - pj.x);
						const auto &h = pi.h;
						if (r < h) {
							const auto pjW = pj.W;
							vect dx = pj.x - pi.x;
							vect xij = pi.x + dx * (pi.h) / (pi.h + pj.h);
							vect mid_dx = xij - pi.x;
							max_dx = max(max_dx, sqrt(mid_dx.dot(mid_dx)));
							max_ngb = max(max_ngb, pjW);
							min_ngb = min(min_ngb, pjW);
						}
					}
				//	const auto beta = max(1.0, min(2.0, 100.0 / pi.Nc));
					const auto beta = 0.5;
					real alpha;
					for (int k = 0; k < STATE_SIZE; k++) {
						const auto dmax_ngb = max_ngb[k] - piW[k];
						const auto dmin_ngb = piW[k] - min_ngb[k];
						real grad_abs = 0.0;
						for (int dim = 0; dim < NDIM; dim++) {
							grad_abs += pi.dW[dim][k] * pi.dW[dim][k];
						}
						grad_abs = sqrt(grad_abs);
						const real den = grad_abs * max_dx;
						if (den == 0.0) {
							alpha = 1.0;
						} else {
							alpha = min(1.0, beta * min(dmax_ngb / den, dmin_ngb / den));
						}
						for (int dim = 0; dim < NDIM; dim++) {
							pi.dW[dim][k] = pi.dW[dim][k] * alpha;
						}

					}
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_gradients_action>(children[ci], t);
		}
		hpx::wait_all(futs);
	}
	parts.resize(nparts0);
}

void tree::compute_conservative_update(fixed_real t, fixed_real dt) {
	static auto opts = options::get();
	if (leaf) {
		for (int i = 0; i < parts.size(); i++) {
			auto &pi = parts[i];
			for (int j = 0; j < parts.size(); j++) {
				auto &pj = parts[j];
				if (i >= nparts0 && j >= nparts0) {
					continue;
				}
				const auto iact = opts.global_time || pi.t + pi.dt == t + dt;
				const auto jact = opts.global_time || pj.t + pj.dt == t + dt;
				if (!iact && !jact) {
					continue;
				}
				if (pi.x > pj.x) {
					const auto r = abs(pj.x - pi.x);
					if (r < std::max(pi.h, pj.h)) {
						const auto xij = (pi.x * pj.h + pj.x * pi.h) / (pi.h + pj.h);
						const auto W_i = pi.W;
						const auto W_j = pj.W;
						auto WL = W_i;
						auto WR = W_j;
						const auto dx = pj.x - pi.x;
						const auto uij = (pi.vf * (pj.x - xij).dot(dx) + pj.vf * (xij - pi.x).dot(dx)) / (dx.dot(dx));
						if (!opts.first_order_time) {
							WL = WL.boost_to(uij);
							WR = WR.boost_to(uij);
							WL = WL + WL.dW_dt(pi.dW, pi.g) * double(fixed_real(0.5) * dt);
							WR = WR + WR.dW_dt(pj.dW, pj.g) * double(fixed_real(0.5) * dt);
							WL = WL.boost_to(-uij);
							WR = WR.boost_to(-uij);
						}
						if (!opts.first_order_space) {
							constexpr auto psi1 = 0.5;
							constexpr auto psi2 = 0.25;
							const auto dxi = xij - pi.x;
							const auto dxj = xij - pj.x;
							for (int dim = 0; dim < NDIM; dim++) {
								WL = WL + pi.dW[dim] * dxi[dim];
								WR = WR + pj.dW[dim] * dxj[dim];
							}
							const auto dW_abs = abs(W_j, W_i);
							const auto delta_1 = dW_abs * psi1;
							const auto delta_2 = dW_abs * psi2;
							const auto W_min = min(W_i, W_j);
							const auto W_max = max(W_i, W_j);
							const auto W_bar_i = W_i + (W_j - W_i) * abs(xij - pi.x) / abs(pi.x - pj.x);
							const auto W_bar_j = W_j + (W_i - W_j) * abs(xij - pj.x) / abs(pi.x - pj.x);
							for (int f = 0; f < STATE_SIZE; f++) {
								real W_m;
								real W_p;
								if ((W_min[f] - delta_1[f]) * W_min[f] >= 0.0) {
									W_m = W_min[f] - delta_1[f];
								} else {
									W_m = W_min[f] * abs(W_min[f]) / (abs(W_min[f]) + delta_1[f]);
								}
								if ((W_max[f] + delta_1[f]) * W_max[f] >= 0.0) {
									W_p = W_max[f] + delta_1[f];
								} else {
									W_p = W_max[f] * abs(W_max[f]) / (abs(W_max[f]) + delta_1[f]);
								}
								if (W_i[f] == W_j[f]) {
									WR[f] = WL[f] = W_i[f];
								} else if (W_i[f] < W_j[f]) {
									WL[f] = max(W_m, min(W_bar_i[f] + delta_2[f], WL[f]));
									WR[f] = min(W_p, max(W_bar_j[f] - delta_2[f], WR[f]));
								} else {
									WL[f] = min(W_p, max(W_bar_i[f] - delta_2[f], WL[f]));
									WR[f] = max(W_m, min(W_bar_j[f] + delta_2[f], WR[f]));
								}
							}
						}
						vect psi_a_ij;
						vect psi_a_ji;
						for (int n = 0; n < NDIM; n++) {
							psi_a_ij[n] = 0.0;
							psi_a_ji[n] = 0.0;
							for (int m = 0; m < NDIM; m++) {
								psi_a_ij[n] = psi_a_ij[n] + pj.B[n][m] * (pi.x[m] - pj.x[m]) * W(r, pj.h) * pj.V;
								psi_a_ji[n] = psi_a_ji[n] + pi.B[n][m] * (pj.x[m] - pi.x[m]) * W(r, pi.h) * pi.V;
							}
						}
						const auto da = psi_a_ji * pi.V - psi_a_ij * pj.V;
						const auto norm = da / abs(da);
						WL = WL.boost_to(uij);
						WR = WR.boost_to(uij);
						WL = WL.rotate_to(norm);
						WR = WR.rotate_to(norm);
						flux_state F;
						riemann_solver(F, WL, WR);
						F = F.rotate_from(norm);
						F = F.boost_from(uij);
						const auto this_dt = opts.global_time ? dt : min(pi.dt, pj.dt);
						if (i < nparts0) {
							pi.Q = pi.Q + F * abs(da) * -double(this_dt);
						}
						if (j < nparts0) {
							pj.Q = pj.Q + F * abs(da) * +double(this_dt);
						}
					}
				}
			}
		}
		parts.resize(nparts0);
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_conservative_update_action>(children[ci], t, dt);
		}
		hpx::wait_all(futs);
	}
}


void tree::con_to_prim(fixed_real t) {
	static auto opts = options::get();
	const auto use_grav = opts.gravity || opts.problem == "kepler";
	if (leaf) {
		PROFILE();
		parts.resize(nparts0);
		for (int i = 0; i < nparts0; i++) {
			auto &p = parts[i];
			if (p.t + p.dt == t || opts.global_time) {
				p.t += p.dt;
				p.con_to_prim();
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<con_to_prim_action>(children[ci], t);
		}
		hpx::wait_all(futs);
	}
}
