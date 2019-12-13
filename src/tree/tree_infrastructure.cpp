#include <octopart/initialize.hpp>
#include <octopart/tree.hpp>

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
			const auto dx = 1.0e-10 * (box.max[dim] - box.min[dim]);
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
			const auto mid = (e + b) * 0.5;
			if (m & 1) {
				child_box.min[dim] = mid;
				child_box.max[dim] = e;
			} else {
				child_box.min[dim] = b;
				child_box.max[dim] = mid;
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
		futs[ci] = hpx::async([child_box, this](std::vector<particle> child_parts) {
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
		neighbors.clear();
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
			}
		}
		if (cparts <= NPART_MAX && all_leaves) {
			std::array<hpx::future<std::vector<particle>>, NCHILD> dfuts;
			for (int ci = 0; ci < NCHILD; ci++) {
				dfuts[ci] = hpx::async<destroy_action>(children[ci]);
			}
			for (int ci = 0; ci < NCHILD; ci++) {
				const auto tmp = dfuts[ci].get();
				children[ci] = hpx::invalid_id;
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
	for (auto &f : cfuts) {
		const auto list = f.get();
		nids.insert(nids.end(), list.begin(), list.end());
	}
	if (nids.size()) {
		finish_tree(std::move(nids));
	}

}

void tree::form_tree(std::vector<hpx::id_type> nids) {
	std::vector<hpx::future<tree_attr>> nfuts(nids.size());
	std::vector<hpx::future<std::array<hpx::id_type, NCHILD>>> cfuts;
	std::vector<tree_attr> attrs(nids.size());
	std::vector<hpx::id_type> next_nids;
	for (int i = 0; i < nids.size(); i++) {
		assert(nids[i] != hpx::invalid_id);
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
		for (auto &f : cfuts) {
			const auto list = f.get();
			next_nids.insert(next_nids.end(), list.begin(), list.end());
		}
		if (nids.size()) {
			finish_tree(std::move(next_nids));
		}
	} else {
		for (int i = 0; i < nids.size(); i++) {
			if (ranges_intersect(attrs[i].box, box)) {
				if (nids[i] != self) {
					if (attrs[i].leaf) {
						next_nids.push_back(nids[i]);
					} else {
						cfuts.push_back(hpx::async<get_children_action>(nids[i]));
					}
				}
			}
		}
		next_nids.insert(next_nids.end(), children.begin(), children.end());
		for (auto &f : cfuts) {
			const auto list = f.get();
			next_nids.insert(next_nids.end(), list.begin(), list.end());
		}
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<form_tree_action>(children[ci], next_nids);
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

void tree::initialize(const std::string &init_name) {
	if (leaf) {
		const auto f = get_initialization_function(init_name);
		for (auto &pi : parts) {
			f(pi);
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<initialize_action>(children[ci], init_name);
		}
		hpx::wait_all(futs);
	}
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


void tree::set_self_and_parent(const hpx::id_type s, const hpx::id_type p) {
	assert(!dead);
	self = std::move(s);
	parent = std::move(p);
	if (!leaf) {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			assert(children[ci] != hpx::invalid_id);
			auto this_p = children[ci];
			futs[ci] = hpx::async<set_self_and_parent_action>(children[ci], std::move(this_p), self);
		}
		hpx::wait_all(futs);
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

void tree::write_checkpoint(const std::string &filename) const {
	FILE *fp;
	if (parent == hpx::invalid_id) {
		fp = fopen(filename.c_str(), "wb");
		fclose(fp);
	}
	if (leaf) {
		fp = fopen(filename.c_str(), "ab");
		for (const auto &part : parts) {
			part.write(fp);
		}
		fclose(fp);
	} else {
		for (int ci = 0; ci < NCHILD; ci++) {
			write_checkpoint_action()(children[ci], filename);
		}
	}
}

void tree::write_silo(int num) const {
	std::string base_name = "Y" + std::to_string(num);
	std::string command = "./check2silo " + base_name;
	write_checkpoint(base_name);
	system(command.c_str());
}
