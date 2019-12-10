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

HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);

tree::tree(std::vector<particle> &&these_parts, const range &box_) :
		box(box_) {
	const int sz = these_parts.size();

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

	/* Refine if too many particles for this box */
	if (sz > NPART_MAX) {
		leaf = false;
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
			int this_sz = these_parts.size();
			std::vector<particle> child_parts;
			for (int i = 0; i < this_sz; i++) {
				auto &part = these_parts[i];
				if (in_range(part.x, child_box)) {
					child_parts.push_back(std::move(part));
					this_sz--;
					these_parts[i] = these_parts[this_sz];
					i--;
				}
			}
			these_parts.resize(this_sz);
			futs[ci] = hpx::async([child_box](std::vector<particle> child_parts) {
				return hpx::new_<tree>(hpx::find_here(), std::move(child_parts), child_box).get();
			}, std::move(child_parts));
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci] = futs[ci].get();
		}

		/* No need to refine, this is a leaf */
	} else {
		leaf = true;
		parts = std::move(these_parts);
	}

}

void tree::compute_gradients() {
	if (leaf) {
		const int nparts0 = parts.size();
		auto tmp = get_neighbor_particles();
		parts.insert(parts.end(), tmp.begin(), tmp.end());

		for (int i = 0; i < nparts0; i++) {
			const auto &pi = parts[i];
			primitive_state max_ngb;
			primitive_state min_ngb;
			const auto piV = pi.to_prim();
			for (int j = 0; j < STATE_SIZE; j++) {
				max_ngb[j] = min_ngb[i] = piV[i];
				for (int dim = 0; dim < NDIM; dim++) {
					grads[i][dim][j] = 0.0;
				}
			}
			for (const auto &pj : parts) {
				const auto r = abs(pi.x - pj.x);
				const auto &h = pi.h;
				if (r < h) {
					const auto pjV = pj.to_prim();
					for (int dim = 0; dim < NDIM; dim++) {
						grads[i][dim] = grads[i][dim] + (pjV - piV) * pj.psi_a[dim];
					}
				}
			}
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
					grad_abs += grads[i][dim][k] * grads[i][dim][k];
				}
				grad_abs = std::sqrt(grad_abs);
				const real den = grad_abs * max_dx;
				if (den == 0.0) {
					alpha = 1.0;
				} else {
					alpha = std::min(1.0, beta * std::min(dmax_ngb / den, dmin_ngb / den));
				}
				for (int dim = 0; dim < NDIM; dim++) {
					grads[i][dim][k] = grads[i][dim][k] * alpha;
				}

			}

		}

		parts.resize(nparts0);
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_gradients_action>(children[ci]);
		}
		hpx::wait_all(futs);
	}
}

void tree::compute_interactions() {
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
					futs[i] = hpx::async<get_particle_positions_action>(neighbors[i], sbox);
				}
				for (int i = 0; i < neighbors.size(); i++) {
					const auto tmp = futs[i].get();
					pos.insert(pos.end(), tmp.begin(), tmp.end());
				}
			}
		}

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
				neighbors.push_back(nids[i]);
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
						neighbors.push_back(nids[i]);
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
	attr.leaf = leaf;
	attr.box = box;
	return attr;
}

std::vector<particle> tree::get_neighbor_particles() {
	const int nparts0 = parts.size();
	std::vector<particle> nparts;
	range sbox;
	std::vector<real> vol(nparts0);
	std::vector<hpx::future<std::vector<particle>>> futs(neighbors.size());
	for (const auto &pi : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			sbox.min[dim] = std::min(sbox.min[dim], pi.x[dim] - pi.h);
			sbox.max[dim] = std::max(sbox.max[dim], pi.x[dim] + pi.h);
		}
	}
	for (int i = 0; i < neighbors.size(); i++) {
		futs[i] = hpx::async<get_particles_action>(neighbors[i], sbox, box);
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
				nparts.push_back(pj);
			}
		}
	}
	return nparts;
}

std::vector<vect> tree::get_particle_positions(const range &search) const {
	const int sz = parts.size();
	std::vector<vect> pos(sz);
	for (int i = 0; i < sz; i++) {
		pos[i] = parts[i].x;
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

std::array<hpx::id_type, NCHILD> tree::get_children() const {
	return children;
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
