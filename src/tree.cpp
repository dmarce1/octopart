/*
 * tree.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#include <octopart/math.hpp>
#include <octopart/tree.hpp>

static constexpr int NPART_MAX = 1000;

#if(NDIM == 1 )
constexpr real CV = 1.0;
constexpr int NNGB = 2;
#else
#if( NDIM == 2 )
constexpr real CV = M_PI;
constexpr int NNGB = 16;
#else
constexpr real CV = 4.0 * M_PI / 3.0;
constexpr int NNGB = 32;
#endif
#endif

HPX_REGISTER_COMPONENT(hpx::components::component<tree>, tree);

tree::tree() {
	mtx = std::make_shared<hpx::lcos::local::mutex>();
	nparts0 = 0;
	dead = false;
	leaf = false;
}

tree::tree(std::vector<particle> &&these_parts, const range &box_) :
		box(box_), nparts0(0), dead(false) {
	const int sz = these_parts.size();

	mtx = std::make_shared<hpx::lcos::local::mutex>();

	/* Create initial box if root */
	if (box == null_range()) {
		for (const auto &part : these_parts) {
			box.min = min(box.min, part.x);
			box.max = max(box.max, part.x);
		}
		for (int dim = 0; dim < NDIM; dim++) {
			const auto dx = 0.01 * (box.max[dim] - box.min[dim]);
			box.min[dim] -= dx;
			box.max[dim] += dx;
		}
	}

	if (sz > NPART_MAX) {
		parts = std::move(these_parts);
		create_children();
	} else {
		leaf = true;
		parts = std::move(these_parts);
	}
}

void tree::compute_drift(real dt) {
	if (leaf) {
		std::vector<std::vector<particle>> send_parts(neighbors.size());
		int sz = parts.size();
		for (int i = 0; i < sz; i++) {
			auto &pi = parts[i];
			pi.x = pi.x + pi.u * dt;
			if (!in_range(pi.x, box)) {
				for (int j = 0; j < neighbors.size(); j++) {
					if (in_range(pi.x, neighbors[j].box)) {
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
		std::vector<hpx::future<void>> futs(neighbors.size());
		for (int j = 0; j < neighbors.size(); j++) {
			futs[j] = hpx::async<send_particles_action>(neighbors[j].id, std::move(send_parts[j]));
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_drift_action>(children[ci], dt);
		}
		hpx::wait_all(futs);
	}

}

void tree::compute_gradients() {
	if (leaf) {
		grad.resize(nparts0);
		grad_lim.resize(nparts0);
		range sbox;
		std::vector<hpx::future<std::vector<particle>>> futs(neighbors.size());
		for (const auto &pi : parts) {
			for (int dim = 0; dim < NDIM; dim++) {
				sbox.min[dim] = std::min(sbox.min[dim], pi.x[dim] - pi.h);
				sbox.max[dim] = std::max(sbox.max[dim], pi.x[dim] + pi.h);
			}
		}
		for (int i = 0; i < neighbors.size(); i++) {
			futs[i] = hpx::async<get_particles_action>(neighbors[i].id, sbox, box);
		}
		for (int i = 0; i < neighbors.size(); i++) {
			const auto these_parts = futs[i].get();
			for (const auto &pj : these_parts) {
				bool inrange = false;
				for (int i = 0; i < nparts0; i++) {
					const auto &pi = parts[i];
					if (abs(pi.x - pj.x) <= std::max(pi.h, pj.h)) {
						inrange = true;
						break;
					}
				}
				if (inrange) {
					parts.push_back(pj);
				}
			}
		}

		for (int i = 0; i < nparts0; i++) {
			const auto &pi = parts[i];
			primitive_state max_ngb;
			primitive_state min_ngb;
			const auto piV = pi.to_prim();
			for (int j = 0; j < STATE_SIZE; j++) {
				max_ngb[j] = min_ngb[i] = piV[i];
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
						grad[i][dim] = grad[i][dim] + (pjV - piV) * pj.psi_a[dim];
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
					max_dx = std::max(max_dx, std::sqrt(mid_dx.dot(mid_dx)));
					max_ngb = max(max_ngb, pjV);
					min_ngb = min(min_ngb, pjV);
				}
			}
			const auto beta = std::max(1.0, std::min(2.0, 100.0 / Ncond[i]));
			real alpha;
			for (int k = 0; k < STATE_SIZE; k++) {
				const auto dmax_ngb = max_ngb[k] - piV[k];
				const auto dmin_ngb = piV[k] - min_ngb[k];
				real grad_abs = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					grad_abs += grad_lim[i][dim][k] * grad_lim[i][dim][k];
				}
				grad_abs = std::sqrt(grad_abs);
				const real den = grad_abs * max_dx;
				if (den == 0.0) {
					alpha = 1.0;
				} else {
					alpha = std::min(1.0, beta * std::min(dmax_ngb / den, dmin_ngb / den));
				}
				for (int dim = 0; dim < NDIM; dim++) {
					grad_lim[i][dim][k] = grad_lim[i][dim][k] * alpha;
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
	if (leaf) {
		nparts0 = parts.size();
		range sbox;
		std::vector<hpx::future<std::vector<gradient>>> futs(neighbors.size());
		for (const auto &pi : parts) {
			for (int dim = 0; dim < NDIM; dim++) {
				sbox.min[dim] = std::min(sbox.min[dim], pi.x[dim] - pi.h);
				sbox.max[dim] = std::max(sbox.max[dim], pi.x[dim] + pi.h);
			}
		}
		for (int i = 0; i < neighbors.size(); i++) {
			futs[i] = hpx::async<get_gradients_action>(neighbors[i].id, sbox, box);
		}
		int cnt = 0;
		for (int i = 0; i < neighbors.size(); i++) {
			const auto these_parts = futs[i].get();
			for (const auto &pj_grad : these_parts) {
				const auto &pj = parts[nparts0 + cnt];
				bool inrange = false;
				for (int i = 0; i < nparts0; i++) {
					const auto &pi = parts[i];
					if (abs(pi.x - pj.x) <= std::max(pi.h, pj.h)) {
						inrange = true;
						break;
					}
				}
				if (inrange) {
					grad.push_back(pj_grad);
					grad_lim.push_back(pj_grad);
				}
				cnt++;
			}
		}

		constexpr auto psi1 = 0.5;
		constexpr auto psi2 = 0.25;
		dudt.resize(nparts0);
		for (auto &du : dudt) {
			for (int i = 0; i < STATE_SIZE; i++) {
				du[i] = 0.0;
			}
		}
		for (int i = 0; i < nparts0; i++) {
			const auto &pi = parts[i];
			for (int j = 0; j < parts.size(); j++) {
				if (i != j) {
					const auto &pj = parts[j];
					const auto r = abs(pj.x - pi.x);
					const auto h = pi.h;
					if (r < h) {
						const auto xij = pi.x + (pj.x - pi.x) * pi.h / (pi.h + pj.h);
						const auto V_i = pi.to_prim();
						const auto V_j = pj.to_prim();
						auto VL = V_i;
						auto VR = V_j;
						const auto dxi = xij - pi.x;
						const auto dxj = xij - pj.x;
						for (int dim = 0; dim < NDIM; dim++) {
							VL = VL + grad_lim[i][dim] * dxi[dim];
							VR = VR + grad_lim[i][dim] * dxj[dim];
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
							if ((V_min[f] - delta_1[f]) * V_min[f] > 0.0) {
								V_m = V_min[f] - delta_1[f];
							} else {
								V_m = V_min[f] * std::abs(V_min[f]) / (std::abs(V_min[f]) + delta_1[f]);
							}
							if ((V_max[f] + delta_1[f]) * V_max[f] > 0.0) {
								V_p = V_max[f] + delta_1[f];
							} else {
								V_p = V_max[f] * std::abs(V_max[f]) / (std::abs(V_max[f]) + delta_1[f]);
							}
							if (V_i[f] == V_j[f]) {
								VR[f] = VL[f] = V_i[f];
							} else if (V_i[f] < V_j[f]) {
								VL[f] = std::max(V_m, std::min(V_bar_i[f] + delta_2[f], VL[f]));
								VR[f] = std::max(V_p, std::min(V_bar_j[f] - delta_2[f], VR[f]));
							} else {
								VL[f] = std::max(V_p, std::min(V_bar_i[f] - delta_2[f], VL[f]));
								VR[f] = std::max(V_m, std::min(V_bar_j[f] + delta_2[f], VR[f]));
							}
						}
						const auto dx = pj.x - pi.x;
						const auto uij = pi.u + (pj.u - pi.u) * (xij - pi.x).dot(dx) / (dx.dot(dx));
						const auto da = pi.psi_a * pi.V - pj.psi_a * pj.V;
						const auto norm = da / abs(da);
						VL = VL.boost_to(uij);
						VR = VR.boost_to(uij);
						VL = VL + VL.dW_dt(grad[i]) * 0.5 * dt;
						VR = VR + VR.dW_dt(grad[j]) * 0.5 * dt;
						VL = VL.rotate_to(norm);
						VR = VR.rotate_to(norm);
						auto F = riemann_solver(VL, VR);
						F = F.rotate_from(norm).boost_from(uij);
						dudt[i] = dudt[i] - F * abs(da);
					}
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_time_derivatives_action>(children[ci], dt);
		}
		hpx::wait_all(futs);
	}
}

real tree::compute_timestep() const {
	real tmin = std::numeric_limits<real>::max();
	if (leaf) {
		for (int i = 0; i < nparts0; i++) {
			const auto &pi = parts[i];
			for (const auto &pj : parts) {
				const auto dx = pi.x - pj.x;
				const auto r = abs(dx);
				const auto &h = pi.h;
				if (r < h) {
					const auto ci = pi.to_prim().sound_speed();
					const auto cj = pj.to_prim().sound_speed();
					const real vsig = ci + cj - std::min(0.0, (pi.u - pj.u).dot(dx) / r);
					tmin = std::min(tmin, h / vsig);
				}
			}
		}
	} else {
		std::array<hpx::future<real>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_timestep_action>(children[ci]);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			tmin = std::min(tmin, futs[ci].get());
		}
	}
	return tmin;
}

void tree::compute_interactions() {
	nparts0 = parts.size();
	if (leaf) {
		std::vector<vect> pos;
		const auto h0 = std::pow(range_volume(box) / (CV * parts.size()), 1.0 / NDIM);
		for (int pass = 0; pass < 2; pass++) {
			for (auto &pi : parts) {
				pos.push_back(pi.x);
				pi.h = h0;
			}
			for (auto &pi : parts) {
				bool done = false;
				auto &h = pi.h;
				double max_dh = std::numeric_limits<real>::max();
				do {
					double N = 0.0;
					double dNdh = 0.0;
					for (const auto &pj : pos) {
						if (pj != pi.x) {
							const auto r = abs(pj - pi.x);
							if (r < h) {
								N += CV * std::pow(h, NDIM) * W(r, h);
								dNdh += CV * NDIM * std::pow(h, NDIM - 1) * W(r, h);
								dNdh += CV * std::pow(h, NDIM) * dW_dh(r, h);
							}
						}
					}
					if (dNdh == 0.0) {
						h *= 2.0;
					} else {
						auto dh = -(N - NNGB) / dNdh;
						dh = std::min(h, std::max(-h / 2.0, dh));
						max_dh = std::min(0.99 * max_dh, std::abs(dh));
						dh = std::copysign(std::min(max_dh, std::abs(dh)), dh);
						h += 0.99 * dh;
					}
					if (std::abs(NNGB - N) < 1e-6) {
						done = true;
					}
				} while (!done);
			}
			range sbox = null_range();
			if (pass == 0) {
				for (const auto &pi : parts) {
					for (int dim = 0; dim < NDIM; dim++) {
						sbox.min[dim] = std::min(sbox.min[dim], pi.x[dim] - pi.h);
						sbox.max[dim] = std::max(sbox.max[dim], pi.x[dim] + pi.h);
					}
				}
				std::vector<hpx::future<std::vector<vect>>> futs(neighbors.size());
				for (int i = 0; i < neighbors.size(); i++) {
					futs[i] = hpx::async<get_particle_positions_action>(neighbors[i].id, sbox);
				}
				for (int i = 0; i < neighbors.size(); i++) {
					const auto tmp = futs[i].get();
					pos.insert(pos.end(), tmp.begin(), tmp.end());
				}
			}
		}
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
				for (int m = 0; m < NDIM; m++) {
					E[n][m] = 0.0;
				}
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
			Ncond[i] = condition_number(E, B);
			for (const auto &pjx : pos) {
				const auto r = abs(pi.x - pjx);
				const auto h = pi.h;
				if (r < h) {
					const auto psi_j = W(r, h) * pi.V;
					for (int n = 0; n < NDIM; n++) {
						real psi_a_j = 0.0;
						for (int m = 0; m < NDIM; m++) {
							psi_a_j += B[n][m] * (pjx[m] - pi.x[m]) * psi_j;
						}
						pi.psi_a[n] = psi_a_j;
					}
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

void tree::compute_next_state(real dt, real beta) {
	if (leaf) {
		parts.resize(nparts0);
		for (int i = 0; i < nparts0; i++) {
			auto U = parts[i].to_con();
			U = U + dudt[i] * dt;
			parts[i] = parts[i].from_con(U);
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_next_state_action>(children[ci], dt, beta);
		}
		hpx::wait_all(futs);
	}
}

int tree::compute_workload() {
	int load;
	if (leaf) {
		std::fill(child_loads.begin(), child_loads.end(), 0);
		load = parts.size();
	} else {
		load = 0;
		std::array<hpx::future<int>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_workload_action>(children[ci]);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto this_load = futs[ci].get();
			child_loads[ci] = this_load;
			load += this_load;
		}
	}
	return load;
}

void tree::create_children() {
	leaf = false;
	nparts0 = 0;
	std::array<hpx::future<hpx::id_type>, NCHILD> futs;
	range child_box;
	for (int ci = 0; ci < NCHILD; ci++) {
		int m = ci;
		for (int dim = 0; dim < NDIM; dim++) {
			const auto &b = box.min[dim];
			const auto &e = box.max[dim];
			if (m & 1) {
				child_box.min[dim] = (e + b) * 0.5;
				child_box.max[dim] = e;
			} else {
				child_box.min[dim] = b;
				child_box.max[dim] = (e + b) * 0.5;
			}
			m >>= 1;
		}
		int this_sz = parts.size();
		std::vector<particle> child_parts;
		for (int i = 0; i < this_sz; i++) {
			auto &part = parts[i];
			if (in_range(part.x, child_box)) {
				child_parts.push_back(std::move(part));
				this_sz--;
				parts[i] = parts[this_sz];
				i--;
			}
		}
		parts.resize(this_sz);
		futs[ci] = hpx::async([child_box](std::vector<particle> child_parts) {
			return hpx::new_<tree>(hpx::find_here(), std::move(child_parts), child_box).get();
		}, std::move(child_parts));
	}
	for (int ci = 0; ci < NCHILD; ci++) {
		children[ci] = futs[ci].get();
	}
}

std::vector<particle> tree::destroy() {
	self = hpx::invalid_id;
	neighbors.clear();
	dead = true;
	return parts;
}

void tree::find_new_neighbors() {
	if (leaf) {
		std::vector<hpx::future<tree_attr>> futs(neighbors.size());
		std::vector<hpx::future<tree_attr>> afuts;
		std::vector<hpx::future<hpx::id_type>> pfuts;
		std::vector<hpx::future<std::array<hpx::id_type, NCHILD>>> cfuts;
		std::vector<neighbor_attr> new_neighbors;
		for (int i = 0; i < neighbors.size(); i++) {
			futs[i] = hpx::async<get_attributes_action>(neighbors[i].id);
		}
		int sz = neighbors.size();
		int new_sz = 0;
		for (int i = 0; i < sz; i++) {
			const auto attr = futs[i].get();
			if (attr.dead || !attr.leaf) {
				if (attr.dead) {
					pfuts.push_back(hpx::async<get_parent_action>(neighbors[i].id));
					new_sz++;
				} else {
					cfuts.push_back(hpx::async<get_children_action>(neighbors[i].id));
					new_sz += NCHILD;
				}
				sz--;
				neighbors[i] = std::move(neighbors[sz]);
				i--;
			}
		}
		neighbors.resize(sz);
		new_neighbors.reserve(new_sz);
		afuts.reserve(pfuts.size());
		for (auto &f : pfuts) {
			const auto id = f.get();
			afuts.push_back(hpx::async<get_attributes_action>(id));
			new_neighbors.push_back( { id, null_range() });
		}
		for (auto &f : cfuts) {
			const auto c = f.get();
			for (int ci = 0; ci < NCHILD; ci++) {
				const auto id = c[ci];
				afuts.push_back(hpx::async<get_attributes_action>(id));
				new_neighbors.push_back( { id, null_range() });
			}
		}
		for (int i = 0; i < afuts.size(); i++) {
			const auto attr = afuts[i].get();
			new_neighbors[i].box = attr.box;
		}
		sz = new_neighbors.size();
		for (int i = 0; i < sz; i++) {
			if (!ranges_intersect(box, new_neighbors[i].box)) {
				sz--;
				new_neighbors[i] = new_neighbors[sz];
				i--;
			}
		}
		new_neighbors.resize(sz);
		neighbors.insert(neighbors.end(), new_neighbors.begin(), new_neighbors.end());
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<find_new_neighbors_action>(children[ci]);
		}
		hpx::wait_all(futs);
	}
}

tree_attr tree::finish_drift() {
	if (leaf) {
		parts.insert(parts.end(), new_parts.begin(), new_parts.end());
		new_parts.clear();
		nparts0 = parts.size();
		if (nparts0 > NPART_MAX) {
			create_children();
		}
	} else {
		std::array<hpx::future<tree_attr>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<finish_drift_action>(children[ci]);
		}
		int cparts = 0;
		bool all_leaves = true;
		for (auto &f : futs) {
			const auto tmp = f.get();
			if (tmp.leaf) {
				cparts += tmp.nparts;
			} else {
				all_leaves = false;
				break;
			}
		}
		if (cparts < NPART_MAX) {
			std::array<hpx::future<std::vector<particle>>, NCHILD> dfuts;
			for (int ci = 0; ci < NCHILD; ci++) {
				dfuts[ci] = hpx::async<destroy_action>(children[ci]);
			}
			for (int ci = 0; ci < NCHILD; ci++) {
				const auto tmp = dfuts[ci].get();
				parts.insert(parts.end(), tmp.begin(), tmp.end());
			}
			leaf = true;
		}
	}
	return get_attributes();
}

void tree::finish_tree(std::vector<hpx::id_type> nids) {
	std::vector<hpx::future<tree_attr>> nfuts(nids.size());
	std::vector<hpx::future<std::array<hpx::id_type, NCHILD>>> cfuts;
	std::vector<tree_attr> attrs(nids.size());
	for (int i = 0; i < nids.size(); i++) {
		nfuts[i] = hpx::async<get_attributes_action>(nids[i]);
	}
	for (int i = 0; i < nids.size(); i++) {
		attrs[i] = nfuts[i].get();
	}
	for (int i = 0; i < nids.size(); i++) {
		if (ranges_intersect(attrs[i].box, box)) {
			if (attrs[i].leaf) {
				neighbors.push_back( { nids[i], attrs[i].box });
			} else {
				cfuts.push_back(hpx::async<get_children_action>(nids[i]));
			}
		}
	}
	nids.resize(0);
	nids.insert(nids.end(), children.begin(), children.end());
	for (auto &f : cfuts) {
		const auto list = f.get();
		nids.insert(nids.end(), list.begin(), list.end());
	}
	if (nids.size()) {
		finish_tree(std::move(nids));
	}

}

void tree::form_tree(const hpx::id_type &self_, const hpx::id_type &parent_, std::vector<hpx::id_type> nids) {
	self = self_;
	parent = parent_;
	std::vector<hpx::future<tree_attr>> nfuts(nids.size());
	std::vector<hpx::future<std::array<hpx::id_type, NCHILD>>> cfuts;
	std::vector<tree_attr> attrs(nids.size());
	for (int i = 0; i < nids.size(); i++) {
		nfuts[i] = hpx::async<get_attributes_action>(nids[i]);
	}
	for (int i = 0; i < nids.size(); i++) {
		attrs[i] = nfuts[i].get();
	}
	neighbors.clear();
	if (leaf) {
		for (int i = 0; i < nids.size(); i++) {
			if (ranges_intersect(attrs[i].box, box)) {
				if (attrs[i].leaf) {
					if (nids[i] != self) {
						neighbors.push_back( { nids[i], attrs[i].box });
					}
				} else {
					cfuts.push_back(hpx::async<get_children_action>(nids[i]));
				}
			}
		}
	} else {
		for (int i = 0; i < nids.size(); i++) {
			if (ranges_intersect(attrs[i].box, box)) {
				cfuts.push_back(hpx::async<get_children_action>(nids[i]));
			} else {
				nids[i] = nids[nids.size() - 1];
				nids.resize(nids.size() - 1);
				i--;
			}
		}
	}
	nids.resize(0);
	if (!leaf) {
		nids.insert(nids.end(), children.begin(), children.end());
	}
	for (auto &f : cfuts) {
		const auto list = f.get();
		nids.insert(nids.end(), list.begin(), list.end());
	}
	if (leaf) {
		if (nids.size()) {
			finish_tree(std::move(nids));
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<form_tree_action>(children[ci], children[ci], self, nids);
		}
		hpx::wait_all(futs);
	}
}

tree_attr tree::get_attributes() const {
	tree_attr attr;
	attr.dead = dead;
	attr.leaf = leaf;
	attr.box = box;
	attr.nparts = parts.size();
	return attr;
}

std::array<hpx::id_type, NCHILD> tree::get_children() const {
	return children;
}

std::vector<gradient> tree::get_gradients(const range &big, const range &small) const {
	std::vector<gradient> gj;
	for (int i = 0; i < nparts0; i++) {
		const auto &pi = parts[i];
		if (in_range(pi.x, big) || ranges_intersect(range_around(pi.x, pi.h), small)) {
			gj.push_back(grad[i]);
			gj.push_back(grad_lim[i]);
		}
	}
	return gj;
}

hpx::id_type tree::get_parent() const {
	return parent;
}

std::vector<vect> tree::get_particle_positions(const range &search) const {
	const int sz = parts.size();
	std::vector<vect> pos;
	for (int i = 0; i < sz; i++) {
		const auto &pi = parts[i];
		if (in_range(pi.x, search)) {
			pos.push_back(parts[i].x);
		}
	}
	return pos;
}

std::vector<particle> tree::get_particles(const range &big, const range &small) const {
	std::vector<particle> pj;
	for (const auto &pi : parts) {
		if (in_range(pi.x, big) || ranges_intersect(range_around(pi.x, pi.h), small)) {
			pj.push_back(pi);
		}
	}
	return pj;
}

void tree::redistribute_workload(int current, int total) {
	if (!leaf) {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			static const auto localities = hpx::find_all_localities();
			const int loc_id = current * localities.size() / total;
			if (localities[loc_id] != hpx::find_here()) {
				hpx::components::migrate<tree>(children[ci], localities[loc_id]);
			}
			futs[ci] = hpx::async<redistribute_workload_action>(children[ci], current, total);
			current += child_loads[ci];
		}
		hpx::wait_all(futs);
	}
}

void tree::send_particles(const std::vector<particle> &pj) {
	std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
	new_parts.insert(new_parts.end(), pj.begin(), pj.end());
}

tree_stats tree::tree_statistics() const {
	tree_stats stats;
	stats.max_level = 0;
	stats.nnodes = 1;
	if (leaf) {
		stats.nparts = parts.size();
		stats.nleaves = 1;
	} else {
		stats.nparts = 0;
		stats.nleaves = 0;
		std::array<hpx::future<tree_stats>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<tree_statistics_action>(children[ci]);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			tree_stats cstat = futs[ci].get();
			stats.max_level = std::max(stats.max_level, cstat.max_level + 1);
			stats.nleaves += cstat.nleaves;
			stats.nnodes += cstat.nnodes;
			stats.nparts += cstat.nparts;
		}
	}
	return stats;
}

void tree::write_checkpoint(const std::string &filename) {
	if (leaf) {
		FILE *fp = fopen(filename.c_str(), "wb");
		for (const auto &part : parts) {
			part.write(fp);
		}
		fclose(fp);
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			write_checkpoint_action()(children[ci], filename);
		}
		hpx::wait_all(futs);
	}
}
