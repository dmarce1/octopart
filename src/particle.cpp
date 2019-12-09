/*
 * particle.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */
#include <octopart/particle.hpp>
#include <octopart/rand.hpp>


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
	for (int i = 0; i < NDIM; i++) {
		r = x[i];
		fwrite(&r, sizeof(real), 1, fp);
	}
	fwrite(&V, sizeof(real), 1, fp);
	fwrite(&h, sizeof(real), 1, fp);
	for (int i = 0; i < STATE_SIZE; i++) {
		r = U[i];
		fwrite(&r, sizeof(real), 1, fp);
	}
}

int particle::read(FILE *fp) {
	real r;
	int rc = 0;
	for (int i = 0; i < NDIM; i++) {
		if (fread(&r, sizeof(real), 1, fp) == 0) {
			rc = -1;
		}
		x[i] = r;
	}
	if (fread(&V, sizeof(real), 1, fp) == 0) {
		rc = -1;
	}
	if (fread(&h, sizeof(real), 1, fp) == 0) {
		rc = -1;
	}
	for (int i = 0; i < STATE_SIZE; i++) {
		if (fread(&r, sizeof(real), 1, fp) == 0) {
			rc = -1;
		}
		U[i] = r;
	}
	return rc;
}
