/*
 * tree.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#include <octopart/math.hpp>
#include <octopart/tree.hpp>

HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);

HPX_REGISTER_ACTION(tree::find_neighbors_action);
HPX_REGISTER_ACTION(tree::form_tree_action);
HPX_REGISTER_ACTION(tree::get_attributes_action);
HPX_REGISTER_ACTION(tree::get_particle_positions_action);
HPX_REGISTER_ACTION(tree::set_self_action);
HPX_REGISTER_ACTION(tree::smoothing_length_init_action);
HPX_REGISTER_ACTION(tree::smoothing_length_iter_action);
HPX_REGISTER_ACTION(tree::tree_statistics_action);
HPX_REGISTER_ACTION(tree::write_checkpoint_action);

static constexpr int NPART_MAX = 100;
#if(NDIM==1)
static constexpr int NNGB = 4;
static constexpr real CV = 1.0;
#elif(NDIM==2)
static constexpr int NNGB = 16;
static constexpr real CV = M_PI;
#else
static constexpr int NNGB = 32;
static constexpr real CV = 4.0 * M_PI / 3.0;
#endif

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
		std::vector<hpx::future<hpx::id_type>> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			int m = ci;
			range child_box;
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
			hpx::future<hpx::id_type> fut = hpx::new_<tree>(hpx::find_here(), std::move(child_parts), child_box);
			futs.push_back(std::move(fut));
		}
		hpx::wait_all(futs);

		/* No need to refine, this is a leaf */
	} else {
		leaf = true;
		parts = std::move(these_parts);
	}

}

void tree::find_neighbors(std::vector<hpx::id_type> list, bool init) {
	std::vector<hpx::future<tree_attr>> futs;
	std::vector<hpx::id_type> next_list;
	if (init) {
		neighbors.clear();
	}
	futs.reserve(list.size());
	for (const auto &id : list) {
		futs.push_back(hpx::async<get_attributes_action>(id));
	}
	for (int i = 0; i < list.size(); i++) {
		auto &fut = futs[i];
		const auto attr = fut.get();
		if (ranges_intersect(box, attr.reach) || ranges_intersect(attr.box, reach)) {
			if (leaf && attr.leaf) {
				neighbors.push_back(list[i]);
			} else {
				next_list.insert(next_list.end(), attr.children.begin(), attr.children.end());
			}
		}
	}
	if (leaf) {
		find_neighbors(std::move(next_list), false);
	} else {
		std::array<hpx::future<void>, NCHILD> cfuts;
		for (int ci = 0; ci < NCHILD; ci++) {
			cfuts[ci] = hpx::async<find_neighbors_action>(children[ci], next_list, true);
		}
		hpx: wait_all(cfuts);
	}

}

void tree::form_tree(const hpx::id_type &self_, const hpx::id_type &parent_) {
	self = self_;
	parent = parent_;
	if (!leaf) {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<form_tree_action>(children[ci], children[ci], self);
		}
		hpx::wait_all(futs);
	}
}

tree_attr tree::get_attributes() const {
	tree_attr attr;
	attr.box = box;
	attr.reach = reach;
	attr.leaf = leaf;
	if (!leaf) {
		attr.children.resize(NCHILD);
		std::copy(children.begin(), children.end(), attr.children.begin());
	}
	return attr;
}

std::vector<vect> tree::get_particle_positions(const range &other_reach) const {
	std::vector<vect> pos;
	for (const auto &pi : parts) {
		if (in_range(pi.x, other_reach)) {
			pos.push_back(pi.x);
		}
	}
	return pos;
}

void tree::set_self(const hpx::id_type &self_) {
	self = self_;
}

void tree::smoothing_length_init() {
	if (leaf) {
		reach = null_range();
		const auto h0 = std::pow(range_volume(box) / NNGB / CV, 1.0 / NDIM);
		for (auto &pi : parts) {
			pi.h = h0;
			for (int dim = 0; dim < NDIM; dim++) {
				reach.min[dim] = std::min(reach.min[dim], pi.x[dim] - pi.h);
				reach.max[dim] = std::max(reach.max[dim], pi.x[dim] + pi.h);
			}
		}
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<smoothing_length_init_action>(children[ci]);
		}
		hpx::wait_all(futs);
	}
}

bool tree::smoothing_length_iter() {
	if (leaf) {
		const int sz = parts.size();
		const int nsz = neighbors.size();
		bool done;
		int iters = 0;
		std::vector<real> max_dh(sz);
		std::vector<real> N(sz);
		std::vector<vect> pjxs;
		std::vector<hpx::future<std::vector<vect>>> futs(nsz);
		for (int ni = 0; ni < nsz; ni++) {
			futs[ni] = hpx::async<get_particle_positions_action>(neighbors[ni], reach);
		}
		for (const auto &pi : parts) {
			pjxs.push_back(pi.x);
		}
		for (int ni = 0; ni < nsz; ni++) {
			const auto pos = futs[ni].get();
			pjxs.insert(pjxs.end(), pos.begin(), pos.end());
		}
		do {
			std::fill(N.begin(), N.end(), 0.0);
			std::fill(max_dh.begin(), max_dh.end(), std::numeric_limits<real>::max());
			reach = null_range();
			done = true;
			for (int i = 0; i < sz; i++) {
				auto &pi = parts[i];
				if (std::abs(N[i] - NNGB) > 0.01) {
					done = false;
					N[i] = 0.0;
					double dNdh = 0.0;
					auto &h = pi.h;
					for (const auto &pjx : pjxs) {
						const auto r = abs(pi.x - pjx);
						if (r < h) {
							const auto w = W(r, h);
							const auto hpow = std::pow(h, NDIM);
							N[i] += CV * hpow * w;
							dNdh += CV * NDIM * hpow * w / h;
							dNdh += CV * hpow * dW_dh(r, h);
						}
					}
					if (dNdh == 0.0) {
						h *= 2.0;
					} else {
						auto dh = -(N[i] - NNGB) / dNdh;
						dh = std::min(h, std::max(-h / 2.0, dh));
						max_dh[i] = std::min(0.99 * max_dh[i], std::abs(dh));
						dh = std::copysign(std::min(max_dh[i], std::abs(dh)), dh);
						h += 0.99 * dh;
					}
					for (int dim = 0; dim < NDIM; dim++) {
						reach.min[dim] = std::min(reach.min[dim], pi.x[dim] - h);
						reach.max[dim] = std::max(reach.max[dim], pi.x[dim] + h);
					}
				}
			}
			iters++;
		} while (!done);
		return iters < 2;
	} else {
		bool done = true;
		std::array<hpx::future<bool>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<smoothing_length_iter_action>(children[ci]);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			done = done && futs[ci].get();
		}
		return done;
	}
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
