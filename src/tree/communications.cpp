#include <octopart/tree.hpp>

using particle_search_action_type = tree::particle_search_action;
using particle_gather_action_type = tree::particle_gather_action;

HPX_REGISTER_ACTION(particle_search_action_type);
HPX_REGISTER_ACTION(particle_gather_action_type);

std::vector<std::shared_ptr<particle>> tree::particle_search(const range &big_box, const range &small_box, const hpx::id_type &caller,
		const hpx::id_type &for_locality) {
	std::vector<std::shared_ptr<particle>> ret_parts;
	if (!leaf) {
		std::array<hpx::future<std::vector<std::shared_ptr<particle>>>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			if (children[ci] != caller) {
				futs[ci] = hpx::async<particle_gather_action>(children[ci], big_box, small_box, for_locality);
			} else {
				futs[ci] = hpx::make_ready_future<std::vector<std::shared_ptr<particle>>>();
			}
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto these_parts = futs[ci].get();
			ret_parts.insert(ret_parts.begin(), these_parts.begin(), these_parts.end());
		}
	} else {
		ret_parts = particle_gather(big_box, small_box, for_locality);
	}
	if (!(in_range(big_box, box) && in_range(small_box, reach)) && parent != hpx::invalid_id) {
		const auto these_parts = hpx::async<particle_search_action>(parent, big_box, small_box, self, for_locality).get();
		ret_parts.insert(ret_parts.begin(), these_parts.begin(), these_parts.end());
	}
	return ret_parts;
}

std::vector<std::shared_ptr<particle>> tree::particle_gather(const range &big_box, const range &small_box, const hpx::id_type &for_locality) {
	std::vector<std::shared_ptr<particle>> ret_parts;
	if (ranges_intersect(box, big_box) || ranges_intersect(reach, small_box)) {
		if (!leaf) {
			std::array<hpx::future<std::vector<std::shared_ptr<particle>>>, NCHILD> futs;
			for (int ci = 0; ci < NCHILD; ci++) {
				futs[ci] = hpx::async<particle_gather_action>(children[ci], big_box, small_box, for_locality);
			}
			for (int ci = 0; ci < NCHILD; ci++) {
				const auto these_parts = futs[ci].get();
				ret_parts.insert(ret_parts.end(), these_parts.begin(), these_parts.end());
			}
		} else {
			for (auto &part : parts) {
				const auto prange = range_around(part.x, part.h);
				if (in_range(part.x, big_box) || ranges_intersect(prange, small_box)) {
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

