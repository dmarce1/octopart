#include <octopart/nindex.hpp>
#include <array>


nindex::operator int() {
	return i;
}

nindex nindex::begin() {
	nindex I;
	I.i = 0;
	return I;
}

nindex nindex::end() {
	nindex I;
	I.i = NNEIGHBOR;
	return I;
}

int nindex::get_dim(int d) const {
	int j = i;
	while (d > 0) {
		j /= 3;
		d--;
	}
	return j % 3 - 1;
}

void nindex::set_dim(int d, int v) {
	std::array<int, NDIM> dims;
	for (int dim = 0; dim < NDIM; dim++) {
		dims[dim] = get_dim(dim);
	}
	dims[d] = v;
	int j = 0;
	for (int dim = 0; dim < NDIM; dim++) {
		j *= 3;
		j += dims[dim] + 1;
	}
}
