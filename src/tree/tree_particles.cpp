/*
 * tree_particles.cpp
 *
 *  Created on: Jan 26, 2020
 *      Author: dmarce1
 */

#include <octopart/initialize.hpp>
#include <octopart/options.hpp>
#include <octopart/profiler.hpp>
#include <octopart/tree.hpp>



void tree::get_neighbor_particles(tree::bnd_ex_type type) {
	static auto opts = options::get();
	if (leaf) {
		parts.resize(nparts0);
		range sbox = null_range();
		for (const auto &pi : parts) {
			for (int dim = 0; dim < NDIM; dim++) {
				sbox.min[dim] = min(sbox.min[dim], pi.x[dim] - pi.h);
				sbox.max[dim] = max(sbox.max[dim], pi.x[dim] + pi.h);
			}
		}
		int k = nparts0;

		switch (type) {
		case NESTING: {
			std::vector<hpx::future<std::vector<nesting_particle>>> futs(siblings.size());
			for (int i = 0; i < siblings.size(); i++) {
				futs[i] = hpx::async<get_nesting_particles_action>(siblings[i].id, sbox, box, siblings[i].pshift);
			}
			for (int i = 0; i < siblings.size(); i++) {
				const auto these_parts = futs[i].get();
				std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
				parts.insert(parts.end(), these_parts.begin(), these_parts.end());
			}
			break;
		}
		case TIMESTEP: {
			std::vector<hpx::future<std::vector<timestep_particle>>> futs(siblings.size());
			for (int i = 0; i < siblings.size(); i++) {
				futs[i] = hpx::async<get_timestep_particles_action>(siblings[i].id, sbox, box, siblings[i].pshift);
			}
			for (int i = 0; i < siblings.size(); i++) {
				const auto these_parts = futs[i].get();
				std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
				parts.insert(parts.end(), these_parts.begin(), these_parts.end());
			}
			break;
		}
		case HYDRO: {
			std::vector<hpx::future<std::vector<hydro_particle>>> futs(siblings.size());
			for (int i = 0; i < siblings.size(); i++) {
				futs[i] = hpx::async<get_hydro_particles_action>(siblings[i].id, sbox, box, siblings[i].pshift);
			}
			for (int i = 0; i < siblings.size(); i++) {
				const auto these_parts = futs[i].get();
				std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
				parts.insert(parts.end(), these_parts.begin(), these_parts.end());
			}
			break;
		}
		case PRIMITIVE: {
			std::vector<hpx::future<std::vector<primitive_particle>>> futs(siblings.size());
			for (int i = 0; i < siblings.size(); i++) {
				futs[i] = hpx::async<get_primitive_particles_action>(siblings[i].id, sbox, box, siblings[i].pshift);
			}
			for (int i = 0; i < siblings.size(); i++) {
				const auto these_parts = futs[i].get();
				std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
				parts.insert(parts.end(), these_parts.begin(), these_parts.end());
			}
			break;
		}
		case ALL: {
			std::vector<hpx::future<std::vector<particle>>> futs(siblings.size());
			for (int i = 0; i < siblings.size(); i++) {
				futs[i] = hpx::async<get_particles_action>(siblings[i].id, sbox, box, siblings[i].pshift);
			}
			for (int i = 0; i < siblings.size(); i++) {
				const auto these_parts = futs[i].get();
				std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
				parts.insert(parts.end(), these_parts.begin(), these_parts.end());
			}
			break;
		}
		default:
			assert(false);
		}

		const bool rbc[3] = { opts.x_reflecting, opts.y_reflecting, opts.z_reflecting };
		for (int i = 0; i < 2 * NDIM; i++) {
			if (rbc[i / 2]) {
				std::vector<particle> rparts;
				const auto sz = parts.size();
				const auto dim = i / 2;
				real axis = i % 2 ? root_box.max[dim] : root_box.min[dim];
				if (axis == (i % 2 ? box.max[dim] : box.min[dim])) {
					const auto rsbox = reflect_range(sbox, dim, axis);
					const auto rbox = reflect_range(box, dim, axis);
					for (int j = 0; j < sz; j++) {
						auto pj = parts[j];
						if (in_range(pj.x, rsbox) || ranges_intersect(range_around(pj.x, pj.h), rbox)) {
							pj.x[dim] = 2.0 * axis - pj.x[dim];
							pj.Q.p[dim] = -pj.Q.p[dim];
							for (int n = 0; n < NDIM; n++) {
								pj.B[n][dim] = -pj.B[n][dim];
								pj.B[dim][n] = -pj.B[dim][n];
							}
							rparts.push_back(pj);
						}
					}
				}
				std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
				parts.insert(parts.end(), rparts.begin(), rparts.end());
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<get_neighbor_particles_action>(children[ci], type);
		}
		hpx::wait_all(futs);
	}
}


std::vector<particle> tree::get_particles(range big, range small, const vect &shift) const {
	std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
	PROFILE();
	std::vector<particle> pj;
	big = shift_range(big, -shift);
	small = shift_range(small, -shift);
	for (int i = 0; i < nparts0; i++) {
		auto pi = parts[i];
		if (in_range(pi.x, big) || ranges_intersect(range_around(pi.x, pi.h), small)) {
			pi.x = pi.x + shift;
			pj.push_back(pi);
		}
	}
	return pj;
}

std::vector<nesting_particle> tree::get_nesting_particles(range big, range small, const vect &shift) const {
	std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
	PROFILE();
	std::vector<nesting_particle> pj;
	big = shift_range(big, -shift);
	small = shift_range(small, -shift);
	for (int i = 0; i < nparts0; i++) {
		auto pi = parts[i];
		if (in_range(pi.x, big) || ranges_intersect(range_around(pi.x, pi.h), small)) {
			pi.x = pi.x + shift;
			pj.push_back(pi);
		}
	}
	return pj;
}


std::vector<timestep_particle> tree::get_timestep_particles(range big, range small, const vect &shift) const {
	std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
	PROFILE();
	std::vector<timestep_particle> pj;
	big = shift_range(big, -shift);
	small = shift_range(small, -shift);
	for (int i = 0; i < nparts0; i++) {
		auto pi = parts[i];
		if (in_range(pi.x, big) || ranges_intersect(range_around(pi.x, pi.h), small)) {
			pi.x = pi.x + shift;
			pj.push_back(pi);
		}
	}
	return pj;
}

std::vector<hydro_particle> tree::get_hydro_particles(range big, range small, const vect &shift) const {
	std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
	PROFILE();
	std::vector<hydro_particle> pj;
	big = shift_range(big, -shift);
	small = shift_range(small, -shift);
	for (int i = 0; i < nparts0; i++) {
		auto pi = parts[i];
		if (in_range(pi.x, big) || ranges_intersect(range_around(pi.x, pi.h), small)) {
			pi.x = pi.x + shift;
			pj.push_back(pi);
		}
	}
	return pj;
}


std::vector<primitive_particle> tree::get_primitive_particles(range big, range small, const vect &shift) const {
	std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
	PROFILE();
	std::vector<primitive_particle> pj;
	big = shift_range(big, -shift);
	small = shift_range(small, -shift);
	for (int i = 0; i < nparts0; i++) {
		auto pi = parts[i];
		if (in_range(pi.x, big) || ranges_intersect(range_around(pi.x, pi.h), small)) {
			pi.x = pi.x + shift;
			pj.push_back(pi);
		}
	}
	return pj;
}


std::vector<vect> tree::get_particle_positions(range search, const vect &shift) const {
	PROFILE();
	const int sz = parts.size();
	std::vector<vect> pos;
	search = shift_range(search, -shift);
	for (int i = 0; i < sz; i++) {
		const auto &pi = parts[i];
		if (in_range(pi.x, search)) {
			pos.push_back(parts[i].x + shift);
		}
	}
	return pos;
}


void tree::send_particles(const std::vector<particle> &pj, const vect &shift) {
	std::lock_guard<hpx::lcos::local::mutex> lock(*mtx);
	PROFILE();
	for (auto p : pj) {
		p.x = p.x - shift;
		new_parts.push_back(std::move(p));
	}
}


