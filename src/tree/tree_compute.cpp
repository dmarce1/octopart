#include <octopart/math.hpp>
#include <octopart/tree.hpp>
#include <octopart/options.hpp>
#include <octopart/profiler.hpp>

#if(NDIM == 1 )
constexpr real CV = 2.0;
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

void tree::compute_drift(fixed_real dt) {
	static const auto opts = options::get();
	if (leaf) {
		std::vector<std::vector<particle>> send_parts(siblings.size());
		{
			PROFILE();
			int sz = parts.size();
			for (int i = 0; i < sz; i++) {
				auto &pi = parts[i];
				pi.x = pi.x + pi.vf * double(dt);
				if (!in_range(pi.x, box)) {
					bool found = false;
					for (int j = 0; j < siblings.size(); j++) {
						if (in_range(pi.x, shift_range(siblings[j].box, siblings[j].pshift))) {
							send_parts[j].push_back(pi);
							found = true;
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
		PROFILE();
		if (!opts.first_order_space) {
			for (int i = 0; i < nparts0; i++) {
				auto &pi = parts[i];
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
				const auto beta = max(1.0, min(2.0, 100.0 / Ncond[i]));
				//	const auto beta = 0.5;
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
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_gradients_action>(children[ci]);
		}
		hpx::wait_all(futs);
	}
}

void tree::compute_time_derivatives(fixed_real dt) {
	static auto opts = options::get();
	if (leaf) {
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
						const auto W_i = pi.W;
						const auto W_j = pj.W;
						auto WL = W_i;
						auto WR = W_j;
						const auto dx = pj.x - pi.x;
						const auto uij = (pi.vf * (pj.x - xij).dot(dx) + pj.vf * (xij - pi.x).dot(dx)) / (dx.dot(dx));
						if (!opts.first_order_time) {
							WL = WL.boost_to(uij);
							WR = WR.boost_to(uij);
							WL = WL + WL.dW_dt(pi.dW) * double(fixed_real(0.5) * dt);
							WR = WR + WR.dW_dt(pj.dW) * double(fixed_real(0.5) * dt);
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
						if (i < nparts0) {
							dudt[i] = dudt[i] - F * abs(da);
						}
						if (j < nparts0) {
							dudt[j] = dudt[j] + F * abs(da);
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

fixed_real tree::compute_timestep(fixed_real t) {
	const static auto opts = options::get();
	fixed_real tmin = fixed_real::max();
	if (leaf) {
		PROFILE();
		for (int i = 0; i < nparts0; i++) {
			const auto &pi = parts[i];
			for (const auto &pj : parts) {
				const auto dx = pi.x - pj.x;
				const auto r = abs(dx);
				if (r > 0.0 && r < max(pi.h, pj.h)) {
					const auto Wi = pi.W;
					const auto Wj = pj.W;
					const auto ci = Wi.sound_speed() + abs(Wi.v);
					const auto cj = Wj.sound_speed() + abs(Wj.v);
					const real vsig = (ci + cj - min(0.0, (pi.Q.p / pi.Q.m - pj.Q.p / pj.Q.m).dot(dx) / r)) / 2.0;
					if (vsig != 0.0) {
						tmin = min(tmin, fixed_real((min(pi.h, pj.h) / vsig).get()));
					}
					if (opts.gravity) {
						const auto a = abs(pi.g);
						if (a > 0.0) {
							tmin = min(tmin, fixed_real(sqrt(pi.h / a).get()));
						}
					}
					//			printf( "%li\n", tmin.get_int());
				}
			}
		}
		tmin *= opts.cfl;
		//	printf( "%li %li %li\n",t.next_bin().get_int(), tmin.get_int(), t.get_int());
		tmin = tmin.nearest_log2();
		tmin = min(tmin, t.next_bin() - t);
		parts.resize(nparts0);
//		sleep(10);
	} else {
		std::array<hpx::future<fixed_real>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_timestep_action>(children[ci], t);
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
					const bool rbc[3] = { opts.x_reflecting, opts.y_reflecting, opts.z_reflecting };
					for (int i = 0; i < 2 * NDIM; i++) {
						if (rbc[i / 2]) {
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

void tree::compute_next_state(fixed_real dt) {
	static auto opts = options::get();
	const auto use_grav = opts.gravity || opts.problem == "kepler";
	if (leaf) {
		PROFILE();
		parts.resize(nparts0);
		for (int i = 0; i < nparts0; i++) {
			auto &p = parts[i];
			auto &Q = p.Q;
			Q = Q + (dudt[i] * double(dt));
			if (Q.m <= 0.0) {
				printf("Negative density! %e\n", Q.m.get());
				abort();
			}
			p.con_to_prim();
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_next_state_action>(children[ci], dt);
		}
		hpx::wait_all(futs);
	}
}
void tree::set_drift_velocity() {
	if (leaf) {
		for (auto &p : parts) {
			p.vf = p.Q.p / p.Q.m;
			//		p.vf = vect(0.0);
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<set_drift_velocity_action>(children[ci]);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}
