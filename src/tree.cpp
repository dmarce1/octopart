/*
 * tree.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#include <octopart/math.hpp>
#include <octopart/tree.hpp>

using write_checkpoint_action_type = tree::write_checkpoint_action;
using first_sweep_action_type = tree::first_sweep_action;
using particle_search_action_type = tree::particle_search_action;
using particle_gather_action_type = tree::particle_gather_action;
using form_tree_action_type = tree::form_tree_action;
using set_self_action_type = tree::set_self_action;

HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);

HPX_REGISTER_ACTION(write_checkpoint_action_type);
HPX_REGISTER_ACTION(first_sweep_action_type);
HPX_REGISTER_ACTION(particle_search_action_type);
HPX_REGISTER_ACTION(particle_gather_action_type);
HPX_REGISTER_ACTION(form_tree_action_type);
HPX_REGISTER_ACTION(set_self_action_type);

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

tree::tree(std::vector<particle> &&these_parts, const range &box_) :
		box(box_) {
	const int sz = these_parts.size();

	/* Create initial box if root */
	if (box == null_range()) {
		for (const auto &part : these_parts) {
			box.first = min(box.first, part.x);
			box.second = max(box.second, part.x);
		}
		for (int dim = 0; dim < NDIM; dim++) {
			const auto dx = 0.01 * (box.second[dim] - box.first[dim]);
			box.first[dim] -= dx;
			box.second[dim] += dx;
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
				const auto &b = box.first[dim];
				const auto &e = box.second[dim];
				if (m & 1) {
					child_box.first[dim] = (e + b) * 0.5;
					child_box.second[dim] = e;
				} else {
					child_box.first[dim] = b;
					child_box.second[dim] = (e + b) * 0.5;
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
					parts[i] = parts[this_sz];
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

real tree::first_sweep() {
	real dt = std::numeric_limits<real>::max();
	if (!leaf) {
		std::array<hpx::future<real>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<first_sweep_action>(children[ci]);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			dt = std::min(dt, futs[ci].get());
		}
	}

	const auto h0 = std::pow(range_volume(box) / NNGB / CV, 1.0 / NDIM);
	for (auto &part : parts) {
		part.h = h0;
	}

	bool done = false;
	const int sz = parts.size();
	std::vector<real> max_dh(sz, std::numeric_limits<real>::max());
	std::vector<real> N(sz, 0.0);
	range search_box(null_range());
	for (auto &part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			auto &min = search_box.first[dim];
			auto &max = search_box.second[dim];
			min = std::min(min, part.x[dim] - part.h);
			max = std::max(max, part.x[dim] + part.h);
		}
	}
	do {
		std::fill(N.begin(), N.end(), 0.0);
		std::fill(max_dh.begin(), max_dh.end(), std::numeric_limits<real>::max());
		nparts = particle_search(search_box, self, hpx::find_here());
		search_box = null_range();
		done = true;
		for (int i = 0; i < sz; i++) {
			auto &pi = parts[i];
			if (std::abs(N[i] - NNGB) > 0.01) {
				done = false;
				N[i] = 0.0;
				double dNdh = 0.0;
				auto &h = pi.h;
				for (auto &pj_ptr : nparts) {
					const auto &pj = *pj_ptr;
					const auto r = abs(pi.x - pj.x);
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
					auto &min = search_box.first[dim];
					auto &max = search_box.second[dim];
					min = std::min(min, pi.x[dim] - h);
					max = std::max(max, pi.x[dim] + h);
				}
			}
		}
	} while (!done);
	return dt;

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

void tree::set_self(const hpx::id_type &self_) {
	self = self_;
}

std::vector<std::shared_ptr<particle>> tree::particle_search(const range &search_box, const hpx::id_type &caller, const hpx::id_type &for_locality) {
	std::vector<std::shared_ptr<particle>> ret_parts;
	if (!leaf) {
		std::array<hpx::future<std::vector<std::shared_ptr<particle>>>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			if (children[ci] != caller) {
				futs[ci] = hpx::async<particle_gather_action>(children[ci], search_box, for_locality);
			} else {
				futs[ci] = hpx::make_ready_future<std::vector<std::shared_ptr<particle>>>();
			}
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto these_parts = futs[ci].get();
			ret_parts.insert(ret_parts.begin(), these_parts.begin(), these_parts.end());
		}
	} else {
		ret_parts = particle_gather(search_box, for_locality);
	}
	if (!in_range(search_box, box) && parent != hpx::invalid_id) {
		const auto these_parts = hpx::async<particle_search_action>(parent, search_box, self, for_locality).get();
		ret_parts.insert(ret_parts.begin(), these_parts.begin(), these_parts.end());
	}
	return ret_parts;
}

std::vector<std::shared_ptr<particle>> tree::particle_gather(const range &search_box, const hpx::id_type &for_locality) {
	std::vector<std::shared_ptr<particle>> ret_parts;
	if (ranges_intersect(box, search_box)) {
		if (!leaf) {
			std::array<hpx::future<std::vector<std::shared_ptr<particle>>>, NCHILD> futs;
			for (int ci = 0; ci < NCHILD; ci++) {
				futs[ci] = hpx::async<particle_gather_action>(children[ci], search_box, for_locality);
			}
			for (int ci = 0; ci < NCHILD; ci++) {
				const auto these_parts = futs[ci].get();
				ret_parts.insert(ret_parts.end(), these_parts.begin(), these_parts.end());
			}
		} else {
			for (auto &part : parts) {
				if (in_range(part.x, search_box)) {
					std::shared_ptr<particle> ptr;
					if (hpx::find_here() == for_locality) {
						ptr = std::shared_ptr<particle>(&part, [](particle*) {
						});
					} else {
						ptr = std::make_shared<particle>(part);
					}
					ret_parts.push_back(ptr);
				}
			}
		}
	}
	return ret_parts;
}
