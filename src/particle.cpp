/*
 * particle.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */
#include <octopart/options.hpp>
#include <octopart/particle.hpp>
#include <octopart/rand.hpp>

std::vector<particle> disc_particle_set(int N) {
	std::vector<particle> rparts;
	const auto cparts = cartesian_particle_set(N);
	for (auto p : cparts) {
		const auto r = abs(p.x);
		if (r < 0.5 && r > 0.0) {
			rparts.push_back(p);
		}
	}
	return rparts;
}

std::vector<particle> cartesian_particle_set(int N) {
	std::vector<particle> parts;
	parts.reserve(std::pow(N, NDIM));
	particle part;
#if(NDIM==3)
	for (int l = 0; l < N; l++) {
		part.x[2] = (l + 0.5) / N - 0.5;
#endif
#if(NDIM>=2)
	for (int j = 0; j < N; j++) {
		part.x[1] = (real(j) + 0.5) / real(N) - real(0.5);
#endif
	for (int i = 0; i < N; i++) {
		part.x[0] = (i + 0.5) / N - 0.5;
		parts.push_back(part);
	}
#if(NDIM>=2)
	}
#endif
#if(NDIM==3)
	}
#endif
	return std::move(parts);
}

std::vector<particle> random_particle_set(int N) {
	std::vector<particle> parts(N);
	for (auto &part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			part.x[dim] = rand_unit_box();
		}
	}
	return std::move(parts);
}

void particle::write(FILE *fp) const {
	real r;
	fwrite(&Q, sizeof(qcon_state), NDIM, fp);
	fwrite(&x, sizeof(real), NDIM, fp);
	fwrite(&g, sizeof(real), NDIM, fp);
	fwrite(&V, sizeof(real), 1, fp);
	fwrite(&c, sizeof(real), 1, fp);
	fwrite(&h, sizeof(real), 1, fp);
	fwrite(&t, sizeof(real), 1, fp);
	fwrite(&dt, sizeof(real), 1, fp);
	fwrite(&B, sizeof(real), NDIM * NDIM, fp);
}

int particle::read(FILE *fp) {
	int cnt = 0;
	cnt += fread(&Q, sizeof(qcon_state), NDIM, fp);
	cnt += fread(&x, sizeof(real), NDIM, fp);
	cnt += fread(&g, sizeof(real), NDIM, fp);
	cnt += fread(&V, sizeof(real), 1, fp);
	cnt += fread(&c, sizeof(real), 1, fp);
	cnt += fread(&h, sizeof(real), 1, fp);
	cnt += fread(&t, sizeof(real), 1, fp);
	cnt += fread(&dt, sizeof(real), 1, fp);
	cnt += fread(&B, sizeof(real), NDIM * NDIM, fp);
	return cnt;
}

void particle::con_to_prim() {
	static const auto opts = options::get();
	W.rho = Q.m / V;
	W.p = (opts.fgamma - 1.0) * (Q.E - Q.p.dot(Q.p) / (2.0 * Q.m)) / V;
	W.v = Q.p / Q.m;
}

//primitive_state particle::to_prim() const {
//	return to_con().to_prim();
//}
//
//conserved_state particle::to_con() const {
//	static const auto fgamma = options::get().fgamma;
//	conserved_state U0;
//	U0.den() = Q.m / V;
//	U0.mom() = Q.p / V;
//	U0.ene() = Q.E / V;
//	return U0;
//}
//
//particle particle::from_con(const conserved_state &U) const {
//	static const auto fgamma = options::get().fgamma;
//	particle p;
//	p.V = V;
//	p.Q.E = U.ene() * V;
//	p.h = h;
//	p.Q.m = U.den() * V;
//	p.Q.p = U.mom() * V;
//	p.B = B;
//	p.x = x;
//	p.g = g;
//	p.vf = vf;
//	return p;
//}

