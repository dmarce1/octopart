/*
 * particle.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */

#include <octopart/particle.hpp>

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
