/*
 * tree.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#include <octopart/tree.hpp>

using write_checkpoint_action_type = tree::write_checkpoint_action;
using form_tree_action_type = tree::form_tree_action;
using set_self_action_type = tree::set_self_action;

HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);

HPX_REGISTER_ACTION(write_checkpoint_action_type);
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
