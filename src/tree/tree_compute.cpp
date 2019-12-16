#include <octopart/math.hpp>
#include <octopart/tree.hpp>
#include <octopart/options.hpp>
#include <octopart/profiler.hpp>

#if(NDIM == 1 )
constexpr real CV = 1.0;
constexpr int NNGB = 4;
#else
#if( NDIM == 2 )
constexpr real CV = M_PI;
constexpr int NNGB = 16;
#else
constexpr real CV = 4.0 * M_PI / 3.0;
constexpr int NNGB = 32;
#endif
#endif

void tree::compute_drift(real dt) {
	if (leaf) {
		std::vector<std::vector<particle>> send_parts(siblings.size());
		{
			PROFILE();
			int sz = parts.size();
			for (int i = 0; i < sz; i++) {
				auto &pi = parts[i];
				pi.x = pi.x + pi.u * dt;
				if (!in_range(pi.x, box)) {
					for (int j = 0; j < siblings.size(); j++) {
						if (in_range(pi.x, shift_range(siblings[j].box, siblings[j].pshift))) {
							send_parts[j].push_back(pi);
							break;
						}
					}
					sz--;
					parts[i] = parts[sz];
					i--;
					parts.resize(sz);
				}
			}
		}
		std::vector<hpx::future<void>> futs(siblings.size());
		for (int j = 0; j < siblings.size(); j++) {
			futs[j] = hpx::async<send_particles_action>(siblings[j].id, std::move(send_parts[j]), siblings[j].pshift);
		}
		hpx::wait_all(futs);
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_drift_action>(children[ci], dt);
		}
		hpx::wait_all(futs);
	}

}

void tree::compute_gradients() {
	static auto opts = options::get();
	if (leaf) {
		assert(nparts0 == parts.size());
		grad.resize(nparts0);
		grad_lim.resize(nparts0);
		range sbox = null_range();
		std::vector<hpx::future<std::vector<particle>>> futs(siblings.size());
		for (const auto &pi : parts) {
			for (int dim = 0; dim < NDIM; dim++) {
				sbox.min[dim] = min(sbox.min[dim], pi.x[dim] - pi.h);
				sbox.max[dim] = max(sbox.max[dim], pi.x[dim] + pi.h);
			}
		}
		for (int i = 0; i < siblings.size(); i++) {
			futs[i] = hpx::async<get_particles_action>(siblings[i].id, sbox, box, siblings[i].pshift);
		}
		for (int i = 0; i < siblings.size(); i++) {
			const auto these_parts = futs[i].get();
			std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
			parts.insert(parts.end(), these_parts.begin(), these_parts.end());
		}
		if (opts.reflecting) {
			for (int i = 0; i < 2 * NDIM; i++) {
				std::vector<particle> rparts;
				const auto sz = parts.size();
				const auto dim = i / 2;
				real axis = i % 2 ? root_box.max[dim] : root_box.min[dim];
				if (axis == (i % 2 ? box.max[dim] : box.min[dim])) {
					const auto rsbox = reflect_range(sbox, dim, axis);
					const auto rbox = reflect_range(box, dim, axis);
					for (int j = 0; j < sz; j++) {
						auto pj = parts[j];
						if (in_range(pj.x, rsbox) || ranges_intersect(range_around(pj.x, pj.h), rbox)) {
							pj.x[dim] = 2.0 * axis - pj.x[dim];
							pj.u[dim] = -pj.u[dim];
							for (int n = 0; n < NDIM; n++) {
								pj.B[n][dim] = -pj.B[n][dim];
								pj.B[dim][n] = -pj.B[dim][n];
							}
							rparts.push_back(pj);
						}
					}
				}
				std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
				parts.insert(parts.end(), rparts.begin(), rparts.end());
			}
		}
		{
			PROFILE();
			if (!opts.first_order_space) {
				for (int i = 0; i < nparts0; i++) {
					const auto &pi = parts[i];
					primitive_state max_ngb;
					primitive_state min_ngb;
					const auto piV = pi.to_prim();
					for (int j = 0; j < STATE_SIZE; j++) {
						max_ngb[j] = min_ngb[j] = piV[j];
						for (int dim = 0; dim < NDIM; dim++) {
							grad[i][dim][j] = 0.0;
						}
					}
					for (const auto &pj : parts) {
						const auto r = abs(pi.x - pj.x);
						const auto &h = pi.h;
						if (r < h) {
							const auto pjV = pj.to_prim();
							for (int dim = 0; dim < NDIM; dim++) {
								real psi_a_j = 0.0;
								for (int m = 0; m < NDIM; m++) {
									psi_a_j += pi.B[dim][m] * (pj.x[m] - pi.x[m]) * W(r, h) * pi.V;
								}
								grad[i][dim] = grad[i][dim] + (pjV - piV) * psi_a_j;
							}
						}
					}
					grad_lim[i] = grad[i];
					real max_dx = 0.0;
					for (const auto &pj : parts) {
						const auto r = abs(pi.x - pj.x);
						const auto &h = pi.h;
						if (r < h) {
							const auto pjV = pj.to_prim();
							vect dx = pj.x - pi.x;
							vect xij = pi.x + dx * (pi.h) / (pi.h + pj.h);
							vect mid_dx = xij - pi.x;
							max_dx = max(max_dx, sqrt(mid_dx.dot(mid_dx)));
							max_ngb = max(max_ngb, pjV);
							min_ngb = min(min_ngb, pjV);
						}
					}
					const auto beta = max(1.0, min(2.0, 100.0 / Ncond[i]));
					real alpha;
					for (int k = 0; k < STATE_SIZE; k++) {
						const auto dmax_ngb = max_ngb[k] - piV[k];
						const auto dmin_ngb = piV[k] - min_ngb[k];
						real grad_abs = 0.0;
						for (int dim = 0; dim < NDIM; dim++) {
							grad_abs += grad_lim[i][dim][k] * grad_lim[i][dim][k];
						}
						grad_abs = sqrt(grad_abs);
						const real den = grad_abs * max_dx;
						if (den == 0.0) {
							alpha = 1.0;
						} else {
							alpha = min(1.0, beta * min(dmax_ngb / den, dmin_ngb / den));
						}
						for (int dim = 0; dim < NDIM; dim++) {
							grad_lim[i][dim][k] = grad_lim[i][dim][k] * alpha;
						}

					}
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_gradients_action>(children[ci]);
		}
		hpx::wait_all(futs);
	}
}

void tree::compute_time_derivatives(real dt) {
	static auto opts = options::get();
	if (leaf) {
		if (!opts.first_order_space) {
			range sbox = null_range();
			std::vector<hpx::future<std::vector<gradient>>> futs(siblings.size());
			for (int i = 0; i < nparts0; i++) {
				const auto &pi = parts[i];
				for (int dim = 0; dim < NDIM; dim++) {
					sbox.min[dim] = min(sbox.min[dim], pi.x[dim] - pi.h);
					sbox.max[dim] = max(sbox.max[dim], pi.x[dim] + pi.h);
				}
			}
			for (int i = 0; i < siblings.size(); i++) {
				futs[i] = hpx::async<get_gradients_action>(siblings[i].id, sbox, box, siblings[i].pshift);
			}
			for (int i = 0; i < siblings.size(); i++) {
				const auto these_grads = futs[i].get();
				std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
				for (int j = 0; j < these_grads.size(); j += 2) {
					grad.push_back(these_grads[j]);
					grad_lim.push_back(these_grads[j + 1]);
				}
			}
			if (opts.reflecting) {
				for (int i = 0; i < 2 * NDIM; i++) {
					std::vector<gradient> rgrad;
					std::vector<gradient> rgrad_lim;
					const auto sz = grad.size();
					const auto dim = i / 2;
					real axis = i % 2 ? root_box.max[dim] : root_box.min[dim];
					if (axis == (i % 2 ? box.max[dim] : box.min[dim])) {
						const auto rsbox = reflect_range(sbox, dim, axis);
						const auto rbox = reflect_range(box, dim, axis);
						for (int j = 0; j < sz; j++) {
							const auto &pj = parts[j];
							if (in_range(pj.x, rsbox) || ranges_intersect(range_around(pj.x, pj.h), rbox)) {
								auto g = grad[j];
								auto gl = grad_lim[j];
								g[dim] = -g[dim];
								gl[dim] = -gl[dim];
								g[dim].vel()[dim] = -g[dim].vel()[dim];
								gl[dim].vel()[dim] = -gl[dim].vel()[dim];
								rgrad.push_back(g);
								rgrad_lim.push_back(gl);
							}
						}
					}
					std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
					grad.insert(grad.end(), rgrad.begin(), rgrad.end());
					grad_lim.insert(grad_lim.end(), rgrad_lim.begin(), rgrad_lim.end());
				}
			}
		}
		PROFILE();
		constexpr auto psi1 = 0.5;
		constexpr auto psi2 = 0.25;
		dudt.resize(nparts0);
		for (auto &du : dudt) {
			for (int i = 0; i < STATE_SIZE; i++) {
				du[i] = 0.0;
			}
		}
		for (int i = 0; i < parts.size(); i++) {
			const auto &pi = parts[i];
			for (int j = 0; j < parts.size(); j++) {
				const auto &pj = parts[j];
				if (i >= nparts0 && j >= nparts0) {
					continue;
				}
				if (pi.x > pj.x) {
					const auto r = abs(pj.x - pi.x);
					if (r < std::max(pi.h, pj.h)) {
						const auto xij = (pi.x * pj.h + pj.x * pi.h) / (pi.h + pj.h);
						const auto V_i = pi.to_prim();
						const auto V_j = pj.to_prim();
						auto VL = V_i;
						auto VR = V_j;
						if (!opts.first_order_space) {
							const auto dxi = xij - pi.x;
							const auto dxj = xij - pj.x;
							for (int dim = 0; dim < NDIM; dim++) {
								VL = VL + grad_lim[i][dim] * dxi[dim];
								VR = VR + grad_lim[j][dim] * dxj[dim];
							}
							const auto dV_abs = abs(V_j, V_i);
							const auto delta_1 = dV_abs * psi1;
							const auto delta_2 = dV_abs * psi2;
							const auto V_min = min(V_i, V_j);
							const auto V_max = max(V_i, V_j);
							const auto V_bar_i = V_i + (V_j - V_i) * abs(xij - pi.x) / abs(pi.x - pj.x);
							const auto V_bar_j = V_j + (V_i - V_j) * abs(xij - pj.x) / abs(pi.x - pj.x);
							for (int f = 0; f < STATE_SIZE; f++) {
								real V_m;
								real V_p;
								if ((V_min[f] - delta_1[f]) * V_min[f] >= 0.0) {
									V_m = V_min[f] - delta_1[f];
								} else {
									V_m = V_min[f] * abs(V_min[f]) / (abs(V_min[f]) + delta_1[f]);
								}
								if ((V_max[f] + delta_1[f]) * V_max[f] >= 0.0) {
									V_p = V_max[f] + delta_1[f];
								} else {
									V_p = V_max[f] * abs(V_max[f]) / (abs(V_max[f]) + delta_1[f]);
								}
								if (V_i[f] == V_j[f]) {
									VR[f] = VL[f] = V_i[f];
								} else if (V_i[f] < V_j[f]) {
									VL[f] = max(V_m, min(V_bar_i[f] + delta_2[f], VL[f]));
									VR[f] = min(V_p, max(V_bar_j[f] - delta_2[f], VR[f]));
								} else {
									VL[f] = min(V_p, max(V_bar_i[f] - delta_2[f], VL[f]));
									VR[f] = max(V_m, min(V_bar_j[f] + delta_2[f], VR[f]));
								}
							}
						}
						const auto dx = pj.x - pi.x;
						const auto uij = (pi.u * (pj.x - xij).dot(dx) + pj.u * (xij - pi.x).dot(dx)) / (dx.dot(dx));
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
						VL = VL.boost_to(uij);
						VR = VR.boost_to(uij);
						if (!opts.first_order_time) {
							VL = VL + VL.dW_dt(grad[i]) * 0.5 * dt;
							VR = VR + VR.dW_dt(grad[j]) * 0.5 * dt;
						}
						auto VL1 = VL.rotate_to(norm);
						auto VR1 = VR.rotate_to(norm);
						auto VL2 = VL.rotate_to(-norm);
						auto VR2 = VR.rotate_to(-norm);
						auto F1 = riemann_solver(VL1, VR1);
						auto F2 = riemann_solver(VR2, VL2);
						F1 = F1.rotate_from(norm);
						F2 = F2.rotate_from(norm);
						flux_state F = (F1 + F2) * 0.5;
						F = F.boost_from(uij);
						if (i < nparts0) {
							dudt[i] = dudt[i] - F * abs(da) / pi.V;
						}
						if (j < nparts0) {
							dudt[j] = dudt[j] + F * abs(da) / pj.V;
						}
					}
				}
			}
		}
		parts.resize(nparts0);
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_time_derivatives_action>(children[ci], dt);
		}
		hpx::wait_all(futs);
	}
}

real tree::compute_timestep() const {
	real tmin = real::max();
	if (leaf) {
		PROFILE();
		for (int i = 0; i < nparts0; i++) {
			const auto &pi = parts[i];
			for (const auto &pj : parts) {
				const auto dx = pi.x - pj.x;
				const auto r = abs(dx);
				if (r > 0.0 && r < max(pi.h, pj.h)) {
					const auto Vi = pi.to_prim();
					const auto Vj = pj.to_prim();
					const auto ci = Vi.sound_speed() + abs(Vi.vel());
					const auto cj = Vj.sound_speed() + abs(Vj.vel());
					const real vsig = (ci + cj - min(0.0, (pi.u - pj.u).dot(dx) / r)) / 2.0;
					tmin = min(tmin, std::min(pi.h, pj.h) / vsig);
				}
			}
		}
	} else {
		std::array<hpx::future<real>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_timestep_action>(children[ci]);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			tmin = min(tmin, futs[ci].get());
		}
	}
	return tmin;
}

void tree::compute_interactions() {
	const auto toler = NNGB * 10.0 * real::eps();
	static auto opts = options::get();
	nparts0 = parts.size();
	if (leaf) {
		if (nparts0) {
			std::vector<vect> pos;
			pos.reserve(2 * parts.size());
			const auto h0 = pow(range_volume(box) / (CV * parts.size()), 1.0 / NDIM);
			for (auto &pi : parts) {
				pos.push_back(pi.x);
				pi.h = h0;
			}
			const auto hmax = box.max[0] - box.min[0];
			for (int pass = 0; pass < 2; pass++) {
				{
					PROFILE();
					for (auto &pi : parts) {
						if (!(pass == 1 && in_range(range_around(pi.x, pi.h), box))) {
							bool done = false;
							auto &h = pi.h;
//						real max_dh = real::max();
							int iters = 0;
							real dh = pi.h / 2.0;
							do {
								real N = 0.0;
								real Np = 0.0;
								real dNdh;
								const auto eps = 0.1 * abs(dh);
								for (const auto &pj : pos) {
									if (pj != pi.x) {
										const auto r = abs(pj - pi.x);
										if (r < h + eps) {
											Np += CV * pow(h + eps, NDIM) * W(r, h + eps);
											if (r < h) {
												N += CV * pow(h, NDIM) * W(r, h);
											}
										}
									}
								}
								if (abs(NNGB - N) < toler) {
									done = true;
								} else {
									dNdh = (Np - N) / eps;
									if (dNdh == 0.0) {
										h *= 1.2;
									} else {
										dh = -(N - NNGB) / dNdh;
										dh = min(h * 0.5, max(-0.5 * h, dh));
										//		max_dh = min(0.999 * max_dh, abs(dh));
										//		dh = copysign(min(max_dh, abs(dh)), dh);
										h += dh;
									}
								}
								iters++;
								if (pass == 0) {
									if (iters > 100 || h > hmax) {
										break;
									}
								} else {
									if (iters > 50) {
										printf("%e %e %e\n", h.get(), dh.get(), /*max_dh.get(), */N.get());
										if (iters == 100) {
											printf("Smoothing length failed to converge\n");
											abort();
										}
									}
								}
							} while (!done);
						}
					}
				}
				if (pass == 0) {
					range sbox = null_range();
					for (const auto &pi : parts) {
						for (int dim = 0; dim < NDIM; dim++) {
							sbox.min[dim] = min(sbox.min[dim], pi.x[dim] - pi.h);
							sbox.max[dim] = max(sbox.max[dim], pi.x[dim] + pi.h);
						}
					}
					std::vector<hpx::future<std::vector<vect>>> futs(siblings.size());
					for (int i = 0; i < siblings.size(); i++) {
						assert(siblings[i].id != hpx::invalid_id);
						futs[i] = hpx::async<get_particle_positions_action>(siblings[i].id, sbox, siblings[i].pshift);
					}
					for (int i = 0; i < siblings.size(); i++) {
						const auto tmp = futs[i].get();
						pos.insert(pos.end(), tmp.begin(), tmp.end());
					}
					if (opts.reflecting) {
						for (int i = 0; i < 2 * NDIM; i++) {
							std::vector<vect> rpos;
							const auto sz = pos.size();
							const auto dim = i / 2;
							range this_box;
							real axis = i % 2 ? root_box.max[dim] : root_box.min[dim];
							if (axis == (i % 2 ? box.max[dim] : box.min[dim])) {
								this_box = reflect_range(sbox, dim, axis);
								if (ranges_intersect(this_box, box)) {
									for (int j = 0; j < sz; j++) {
										const auto pix = pos[j];
										if (in_range(pix, this_box)) {
											auto this_x = pix;
											this_x[dim] = 2.0 * axis - this_x[dim];
											rpos.push_back(this_x);
										}
									}
								}
							}
							pos.insert(pos.end(), rpos.begin(), rpos.end());
						}
					}
				}
			}
			{
				PROFILE();
				Ncond.resize(parts.size());
				for (int i = 0; i < parts.size(); i++) {
					auto &pi = parts[i];
					pi.V = 0.0;
					for (const auto &pjx : pos) {
						const auto r = abs(pi.x - pjx);
						const auto h = pi.h;
						if (r < h) {
							pi.V += W(r, h);
						}
					}
					pi.V = 1.0 / pi.V;
					std::array<vect, NDIM> E, B;
					for (int n = 0; n < NDIM; n++) {
						E[n] = vect(0.0);
					}
					for (const auto &pjx : pos) {
						const auto r = abs(pi.x - pjx);
						const auto h = pi.h;
						if (r < h) {
							const auto psi_j = W(r, h) * pi.V;
							for (int n = 0; n < NDIM; n++) {
								for (int m = 0; m < NDIM; m++) {
									E[n][m] += (pjx[n] - pi.x[n]) * (pjx[m] - pi.x[m]) * psi_j;
								}
							}
						}
					}
					Ncond[i] = condition_number(E, pi.B);
					assert(Ncond[i] != 0.0);
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_interactions_action>(children[ci]);
		}
		hpx::wait_all(futs);
	}
}

void tree::compute_next_state(real dt) {
	static auto opts = options::get();
	if (leaf) {
		PROFILE();
		parts.resize(nparts0);
		if (opts.problem == "kepler") {
			for (int i = 0; i < nparts0; i++) {
				const auto &p = parts[i];
				const auto &x = p.x;
				const auto r = abs(x);
				for (int dim = 0; dim < NDIM; dim++) {
					const auto d = p.m / p.V * x[dim] / (r * r * r);
					dudt[i].mom()[dim] -= d;
					dudt[i].ene() -= d * p.u[dim];
				}
			}
		}
		for (int i = 0; i < nparts0; i++) {
			auto U = parts[i].to_con();
			U = U + dudt[i] * dt;
			parts[i] = parts[i].from_con(U);
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_next_state_action>(children[ci], dt);
		}
		hpx::wait_all(futs);
	}
}
