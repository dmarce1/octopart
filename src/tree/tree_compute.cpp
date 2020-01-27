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

bool tree::adjust_timesteps(fixed_real t, int factor) {
	static const auto opts = options::get();
	bool rc = false;
	if (leaf) {
		for (int i = 0; i < nparts0; i++) {
			auto &pi = parts[i];
			pi.tmp = pi.dt;
			if (pi.t == t || opts.global_time) {
				for (const auto &pj : parts) {
					const auto r = abs(pi.x - pj.x);
					const auto &h = pi.h;
					if (r < h) {
						if (pi.tmp > pj.dt * fixed_real(factor)) {
							rc = true;
							pi.tmp = min(pi.tmp, fixed_real(factor) * pj.dt);
						}
					}
				}
			}
		}
		for (int i = 0; i < nparts0; i++) {
			auto &pi = parts[i];
			if (pi.t == t || opts.global_time) {
				pi.dt = pi.tmp;
			}
		}
		parts.resize(nparts0);
	} else {
		std::array<hpx::future<bool>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<adjust_timesteps_action>(children[ci], t, factor);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			rc = rc || futs[ci].get();
		}
	}
	return rc;
}

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

fixed_real tree::compute_timestep(fixed_real t) {
	const static auto opts = options::get();
	fixed_real tmin = fixed_real::max();
	if (leaf) {
		PROFILE();
		for (int i = 0; i < nparts0; i++) {
			auto &pi = parts[i];
			if (pi.t == t || opts.global_time) {
				pi.m0 = pi.Q.m;
				pi.g0 = pi.g;
				pi.dt = fixed_real::max();
				const auto Wi = pi.W;
				const auto ci = Wi.sound_speed();
				const auto ai = ci + abs(pi.Q.p / pi.Q.m - pi.vf);
				if (ai != 0.0) {
					pi.dt = min(pi.dt, fixed_real((pi.h / ai).get()));
				}
				//			printf( "%e %e %e\n", pi.h / ai, pi.h, ai);
				if (opts.gravity) {
					const auto a = abs(pi.g);
					if (a > 0.0) {
						pi.dt = double(min(pi.dt, fixed_real(sqrt(pi.h / a).get())));
					}
				}
				for (const auto &pj : parts) {
					const auto dx = pi.x - pj.x;
					const auto r = abs(dx);
					if (r > 0.0 && r < max(pi.h, pj.h)) {
						const auto Wj = pj.W;
						const auto ci = Wi.sound_speed() + abs(Wi.v);
						const auto cj = Wj.sound_speed() + abs(Wj.v);
						const real vsig = (ci + cj - min(0.0, (pi.vf - pj.vf).dot(dx) / r));
						if (vsig != 0.0) {
							pi.dt = min(pi.dt, fixed_real((r / vsig).get()));
						}
					}
				}
				pi.dt *= opts.cfl;
				if (pi.dt.nearest_log2() < fixed_real(0)) {
					printf("%e %e\n", (double) pi.dt, (double) pi.dt.nearest_log2());
				}
				pi.dt = pi.dt.nearest_log2();
				pi.dt = min(pi.dt, t.next_bin() - t);
				tmin = min(tmin, pi.dt);
			}
		}
		parts.resize(nparts0);
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

void tree::compute_interactions(fixed_real t, fixed_real dt) {
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
				if (pi.h == -1) {
					pi.h = h0;
				}
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
								const auto eps = abs(dh);
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
									if (iters >= 50) {
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
				for (int i = 0; i < parts.size(); i++) {
					auto &pi = parts[i];
					if (pi.t + pi.dt == t + dt || opts.global_time) {
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

						pi.Nc = condition_number(E, pi.B);
						assert(pi.Nc != 0.0);
					}
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_interactions_action>(children[ci], t, dt);
		}
		hpx::wait_all(futs);
	}
}


void tree::set_drift_velocity(fixed_real t) {
	static const auto opts = options::get();
	if (leaf) {
		for (int i = 0; i < nparts0; i++) {
			auto &pi = parts[i];
			if (pi.t == t || opts.global_time == 1) {
				auto w = pi.Q.p / pi.Q.m;
				if (opts.eulerian) {
					pi.vf = vect(0.0);
				} else {
					pi.vf = w;
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<set_drift_velocity_action>(children[ci], t);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}
