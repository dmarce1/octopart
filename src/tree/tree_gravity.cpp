#include <octopart/math.hpp>
#include <octopart/tree.hpp>

mass_attr tree::compute_mass_attributes() {
	auto& Xcom = mass.com;
	auto& mtot = mass.mtot;
	auto& rmaxs = mass.rmaxs;
	auto& rmaxb = mass.rmaxb;
	Xcom = vect(0);
	mtot = 0.0;
	rmaxs = 0.0;
	rmaxb = 0.0;
	mass.leaf = leaf;
	if (leaf) {
		for (const auto &p : parts) {
			Xcom = Xcom + p.x * p.m;
			mtot += p.m;
		}
		Xcom = Xcom / mtot;
		for (const auto &p : parts) {
			rmaxb = max(rmaxb, abs(p.x - Xcom));
		}
	} else {
		std::array<hpx::future<mass_attr>, NCHILD> futs;
		std::array<mass_attr, NCHILD> child_attr;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_mass_attributes_action>(children[ci]);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			child_attr[ci] = futs[ci].get();
			const auto &c = child_attr[ci];
			Xcom = Xcom + c.com * c.mtot;
			mtot += c.mtot;
		}
		Xcom = Xcom / mtot;
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto &c = child_attr[ci];
			rmaxb = max(rmaxb, c.rmaxb + abs(Xcom - c.com));
		}
	}
	vect Xv;
	for (int i = 0; i < NCHILD; i++) {
		int m = i;
		for (int dim = 0; dim < NDIM; dim++) {
			if (m & 1) {
				Xv[dim] = box.min[dim];
			} else {
				Xv[dim] = box.max[dim];
			}
			m >>= 1;
		}
		rmaxs = max(rmaxs, abs(Xv - Xcom));
	}
	return mass;
}

std::vector<gravity_part> tree::get_gravity_particles() const {
	std::vector<gravity_part> gparts(parts.size());
	for (int i = 0; i < parts.size(); i++) {
		gparts[i].m = parts[i].m;
		gparts[i].x = parts[i].x;
		gparts[i].h = parts[i].h;
	}
	return gparts;
}

mass_attr tree::get_mass_attributes() const {
	return mass;
}

void tree::compute_gravity(std::vector<hpx::id_type> nids, std::vector<mass_attr> masses) {
	constexpr real theta = 0.25;
	constexpr real G = 1.0;
	std::vector<hpx::future<mass_attr>> futs;
	std::vector<hpx::future<std::array<hpx::id_type, NCHILD>>> ncfuts;
	for (const auto &n : nids) {
		futs.push_back(hpx::async<get_mass_attributes_action>(n));
	}
	const auto rmaxA = min(mass.rmaxb, mass.rmaxs);
	const auto ZA = mass.com;
	if (leaf) {
		std::vector<hpx::id_type> near;
		ncfuts.clear();
		while (nids.size()) {
			for (int i = 0; i < nids.size(); i++) {
				const auto tmp = futs[i].get();
				const auto rmaxB = min(tmp.rmaxb, tmp.rmaxs);
				const auto ZB = tmp.com;
				if (abs(ZA - ZB) > (rmaxA + rmaxB) / theta) {
					masses.push_back(tmp);
				} else if (tmp.leaf) {
					near.push_back(nids[i]);
				} else {
					ncfuts.push_back(hpx::async<get_children_action>(nids[i]));
				}
				nids.clear();
				for (auto &f : ncfuts) {
					const auto tmp = f.get();
					nids.insert(nids.end(), tmp.begin(), tmp.end());
				}
			}
		}
		std::vector<hpx::future<std::vector<gravity_part>>> gfuts(near.size());
		for (int i = 0; i < near.size(); i++) {
			gfuts[i] = hpx::async<get_gravity_particles_action>(near[i]);
		}
		gforce.resize(parts.size(), vect(0.0));
		for (int i = 0; i < parts.size(); i++) {
			const auto &pi = parts[i];
			for (int j = 0; j < masses.size(); j++) {
				const auto r = masses[j].com - pi.x;
				const auto r3inv = pow(abs(r), -3);
				gforce[i] = gforce[i] - r * (G * masses[i].mtot * r3inv);
			}
		}
		for (auto &n : gfuts) {
			const auto list = n.get();
			for (int i = 0; i < parts.size(); i++) {
				const auto &pi = parts[i];
				for (int j = 0; j < list.size(); j++) {
					const auto r = list[j].x - pi.x;
					const auto h = 0.5 * (pi.h + list[j].h);
					const auto f = grav_force(abs(r), h);
					gforce[i] = gforce[i] + r * f * G * list[i].m;
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> cfuts;
		std::vector<hpx::id_type> leaf_nids;
		for (int i = 0; i < nids.size(); i++) {
			const auto tmp = futs[i].get();
			const auto rmaxB = min(tmp.rmaxb, tmp.rmaxs);
			const auto ZB = tmp.com;
			if (abs(ZA - ZB) > (rmaxA + rmaxB) / theta) {
				masses.push_back(tmp);
			} else if (tmp.leaf) {
				leaf_nids.push_back(nids[i]);
			} else {
				ncfuts.push_back(hpx::async<get_children_action>(nids[i]));
			}
		}
		nids.clear();
		nids.insert(nids.end(), leaf_nids.begin(), leaf_nids.end());
		for (auto &c : ncfuts) {
			const auto tmp = c.get();
			nids.insert(nids.end(), tmp.begin(), tmp.end());
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			cfuts[ci] = hpx::async<compute_gravity_action>(children[ci], nids, masses);
		}
		hpx::wait_all(cfuts);
	}
}
