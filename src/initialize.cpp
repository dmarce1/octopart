#include <octopart/initialize.hpp>

void drift(particle &p) {
	p.m = 1.0;
	p.e = 0.0;
	const auto r = abs(p.x);
	p.u = -p.x;
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

init_func_type get_initialization_function(const std::string &name) {
	if (name == "drift") {
		return drift;
	} else if (name == "sod") {
		return sod;
	} else {
		printf("Error: Initialization function %s is not known\n", name.c_str());
		abort();
	}
}
