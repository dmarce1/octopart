/*
 * particle.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */
#include <octopart/particle.hpp>
#include <octopart/rand.hpp>

std::vector<particle> disc_particle_set(int N) {
	std::vector<particle> rparts;
	const auto cparts = cartesian_particle_set(N);
	for( auto p : cparts) {
		const auto r = abs(p.x);
		if( r < 0.45 && r > 0.05 ) {
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
		part.x[1] = (j + 0.5) / N - 0.5;
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
	fwrite(&x, sizeof(real), NDIM, fp);
	fwrite(&u, sizeof(real), NDIM, fp);
	fwrite(&g, sizeof(real), NDIM, fp);
	fwrite(&m, sizeof(real), 1, fp);
	fwrite(&e, sizeof(real), 1, fp);
	fwrite(&V, sizeof(real), 1, fp);
	fwrite(&h, sizeof(real), 1, fp);
	fwrite(&B, sizeof(real), NDIM * NDIM, fp);
}

int particle::read(FILE *fp) {
	int cnt = 0;
	cnt += fread(&x, sizeof(real), NDIM, fp);
	cnt += fread(&u, sizeof(real), NDIM, fp);
	cnt += fread(&g, sizeof(real), NDIM, fp);
	cnt += fread(&m, sizeof(real), 1, fp);
	cnt += fread(&e, sizeof(real), 1, fp);
	cnt += fread(&V, sizeof(real), 1, fp);
	cnt += fread(&h, sizeof(real), 1, fp);
	cnt += fread(&B, sizeof(real), NDIM * NDIM, fp);
	return cnt;
}

primitive_state particle::to_prim() const {
	return to_con().to_prim();
}

conserved_state particle::to_con() const {
	conserved_state U;
	U.den() = m / V;
	U.mom() = u * (m / V);
	U.ene() = e / V;
	return U;
}

particle particle::from_con(const conserved_state &U) const {
	particle p;
	p.V = V;
	p.e = U.ene() * V;
	p.h = h;
	p.m = U.den() * V;
	p.u = U.mom() / U.den();
	p.B = B;
	p.x = x;
	return p;
}

