/*
 * tree.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#include <octopart/math.hpp>
#include <octopart/tree.hpp>

static constexpr int NPART_MAX = 1000;

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

std::array<hpx::id_type, NCHILD> tree::get_children() const {
	return children;
}

void tree::set_self(const hpx::id_type &self_) {
	self = self_;
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
