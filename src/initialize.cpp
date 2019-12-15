#include <octopart/initialize.hpp>
#include <octopart/rand.hpp>

void drift(particle &p) {
	p.m = 1.0;
	p.e = 0.0;
	const auto r = abs(p.x);
	p.u = vect(0);
	p.u[0] = 1.0;

}

void kh(particle &p) {
	for (int dim = 0; dim < NDIM; dim++) {
		p.u[dim] = rand_unit_box() * 0.01;
	}
	if (abs(p.x[1]) > 0.25) {
		p.u[0] += -0.5;
		p.m = p.V;
		p.e = p.V + p.u.dot(p.u) * 0.5 * p.m;
	} else {
		p.u[0] += +0.5;
		p.m = 2.0 * p.V;
		p.e = p.V + p.u.dot(p.u) * 0.5 * p.m; ;
	}
}

void sod(particle &p) {
	if (p.x[0] > 0.0) {
		p.m = 1.0 * p.V;
		p.e = 2.5 * p.V;
	} else {
		p.m = 0.125 * p.V;
		p.e = 0.25 * p.V;
	}
	for (int dim = 0; dim < NDIM; dim++) {
		p.u[dim] = 0.0;
	}
}

void blast(particle &p) {
	const auto r = sqrt(p.x.dot(p.x));
	p.m = p.V;
	p.e = exp(-500.0 * r) * p.V;
	for (int dim = 0; dim < NDIM; dim++) {
		p.u[dim] = 0.0;
	}
}

init_func_type get_initialization_function(const std::string &name) {
	if (name == "drift") {
		return drift;
	} else if (name == "kh") {
		return kh;
	} else if (name == "blast") {
		return blast;
	} else if (name == "sod") {
		return sod;
	} else {
		printf("Error: Initialization function %s is not known\n", name.c_str());
		abort();
	}
}
