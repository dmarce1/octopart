/*
 * math.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */

#include <octopart/range.hpp>

bool in_range(const vect &x, const range &r) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (x[dim] < r.min[dim] || x[dim] >= r.max[dim]) {
			return false;
		}
	}
	return true;
}

range shift_range(const range& r_, int d, real s) {
	range r = r_;
	r.min[d] += s;
	r.max[d] += s;
	return r;
}


bool in_range(const range &a, const range &b) {
	return in_range(a.min, b) && in_range(a.max, b);
}

real range_volume(const range &r) {
	real v = 1.0;
	for (int dim = 0; dim < NDIM; dim++) {
		v *= r.max[dim] - r.min[dim];
	}
	return v;
}

bool ranges_intersect(const range &a, const range &b) {
	for (int dim = 0; dim < NDIM; dim++) {
		const auto front = max(a.min[dim], b.min[dim]);
		const auto back = min(a.max[dim], b.max[dim]);
		if (front > back) {
			return false;
		}
	}
	return true;
}

range null_range() {
	range null;
	for (int dim = 0; dim < NDIM; dim++) {
		null.min[dim] = real::max();
		null.max[dim] = -real::max();
	}
	return null;
}

range range_around(const vect &p, real h) {
	range r;
	for (int dim = 0; dim < NDIM; dim++) {
		r.max[dim] = p[dim] + h;
		r.min[dim] = p[dim] - h;
	}
	return r;
}

range expand_range(const range &rs, real h) {
	range rb;
	for (int dim = 0; dim < NDIM; dim++) {
		rb.min[dim] = rs.min[dim] + h;
		rb.max[dim] = rs.max[dim] - h;
	}
	return rb;
}

bool operator==(const range &a, const range &b) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (a.min[dim] != b.min[dim]) {
			return false;
		}
		if (a.max[dim] != b.max[dim]) {
			return false;
		}
	}
	return true;
}

bool operator!=(const range &a, const range &b) {
	return !(a == b);
}
