/*
 * math.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */

#include <octopart/range.hpp>

bool in_range(const vect &x, const range &r) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (x[dim] < r.first[dim] || x[dim] >= r.second[dim]) {
			return false;
		}
	}
	return true;
}

bool in_range(const range &a, const range &b) {
	return in_range(a.first, b) && in_range(a.second, b);
}

real range_volume(const range& r) {
	real v = 1.0;
	for( int dim = 0; dim < NDIM; dim++) {
		v *= r.second[dim] - r.first[dim];
	}
	return v;
}


bool ranges_intersect(const range &a, const range &b) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (a.first[dim] < b.first[dim]) {
			if (a.second[dim] < b.first[dim]) {
				return false;
			}
		} else {
			if (b.second[dim] < a.first[dim]) {
				return false;
			}
		}
	}
	return true;
}

range null_range() {
	range null;
	for (int dim = 0; dim < NDIM; dim++) {
		null.first[dim] = std::numeric_limits<real>::max();
		null.second[dim] = -std::numeric_limits<real>::max();
	}
	return null;
}
