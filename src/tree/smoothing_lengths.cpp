
#include <octopart/math.hpp>
#include <octopart/tree.hpp>

using find_smoothing_lengths_action_type = tree::find_smoothing_lengths_action;

HPX_REGISTER_ACTION(find_smoothing_lengths_action_type);

range tree::find_smoothing_lengths() {
	if (!leaf) {
		reach = null_range();
		std::array<hpx::future<range>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<find_smoothing_lengths_action>(children[ci]);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto this_reach = futs[ci].get();
			for (int dim = 0; dim < NDIM; dim++) {
				reach.min[dim] = std::min(reach.min[dim], this_reach.min[dim]);
				reach.max[dim] = std::max(reach.max[dim], this_reach.max[dim]);
			}
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

	reach = null_range();
	for (auto &part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			auto &min = reach.min[dim];
			auto &max = reach.max[dim];
			min = std::min(min, part.x[dim] - part.h);
			max = std::max(max, part.x[dim] + part.h);
		}
	}
	do {
		std::fill(N.begin(), N.end(), 0.0);
		std::fill(max_dh.begin(), max_dh.end(), std::numeric_limits<real>::max());
		nparts = particle_search(reach, box, self, hpx::find_here());
		reach = null_range();
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
					auto &min = reach.min[dim];
					auto &max = reach.max[dim];
					min = std::min(min, pi.x[dim] - h);
					max = std::max(max, pi.x[dim] + h);
				}
			}
		}
	} while (!done);
	return reach;

}
