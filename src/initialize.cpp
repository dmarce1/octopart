#include <octopart/initialize.hpp>
#include <octopart/options.hpp>
#include <octopart/rand.hpp>

void drift(particle &p) {
	p.m = 1.0;
	p.E = 0.0;
	const auto r = abs(p.x);
	p.v = vect(0);
	p.v[0] = rand_unit_box();
	p.v[1] = rand_unit_box();

}

void kh(particle &p) {
	const real x = p.x[0];
	const real y = p.x[1];
	constexpr real rho1 = 1.0;
	constexpr real rho2 = 2.0;
	const real drho = 0.5 * (rho2 - rho1);
	constexpr real dy = 0.025;
	constexpr real E = 3.75;
	real rho, vx, vy;
	if (y < -0.25) {
		rho = rho1 + drho * exp((y + 0.25) / dy);
		vx = -0.5 + 0.5 * exp((y + 0.25) / dy);
	} else if (y < 0.) {
		rho = rho2 - drho * exp((-0.25 - y) / dy);
		vx = +0.5 - 0.5 * exp((-0.25 - y) / dy);
	} else if (y < 0.25) {
		rho = rho2 - drho * exp((y - 0.25) / dy);
		vx = +0.5 - 0.5 * exp((y - 0.25) / dy);
	} else {
		rho = rho1 + drho * exp((0.25 - y) / dy);
		vx = -0.5 + 0.5 * exp((0.25 - y) / dy);
	}
	vy = 0.01 * sin(4.0 * M_PI * (x + 0.5));
	p.v[0] = -vx;
	p.v[1] = vy;
#if(NDIM==3)
	p.v[2] = 0.0;
#endif
	p.m = rho * p.V;
	p.E = E * p.V + p.v.dot(p.v) * 0.5 * p.m;
	;
//	for (int dim = 0; dim < NDIM; dim++) {
//		p.v[dim] = rand_unit_box() * 0.01;
//	}
//	if (abs(p.x[1]) > 0.25) {
//		p.v[0] += +0.5;
//		p.m = p.V;
//		p.E = 6.25 * p.V + p.v.dot(p.v) * 0.5 * p.m;
//	} else {
//		p.v[0] += -0.5;
//		p.m = 2.0 * p.V;
	//		p.E = 6.25 * p.V + p.v.dot(p.v) * 0.5 * p.m; ;
//	}
}

void kepler(particle &p) {
	static const auto opts = options::get();
	static const auto eps = opts.kep_eps;
	const auto r = abs(p.x);
	const auto y = p.x[1];
	const auto x = p.x[0];
	p.v = vect(0);
	const auto tmp = pow(eps * eps + r * r, -0.75);
	p.v[0] = -y * tmp;
	p.v[1] = +x * tmp;
	if (r < 1.0 / 12.0) {
		p.m = 0.01 + 12 * r * r * r;
	} else if (r < 1.0 / 3.0) {
		p.m = 0.01 + 1.0;
	} else {
		p.m = 0.01 + 1.0 / pow(1.0 + (r - 1.0 / 3.0) / .1, 3.0);
	}
	p.m *= p.V;
	p.E = 1.0e-6 * p.V + 0.5 * p.v.dot(p.v) * p.m;
}

void sod(particle &p) {
	for (int dim = 0; dim < NDIM; dim++) {
		p.v[dim] = 0.0;
	}
	if (p.x[0] > 0.0) {
		p.m = 1.0 * p.V;
		p.v[0] = 1;
		p.E = 2.5 + p.v.dot(p.v) / 2.0 * p.m * p.V;
	} else {
		p.m = 1.0 * p.V;
		p.v[0] = -1;
		p.E = 0.25 + p.v.dot(p.v) / 2.0 * p.m * p.V;
	}
}

void blast(particle &p) {
	const auto r = sqrt(p.x.dot(p.x));
	p.m = p.V;
	for (int dim = 0; dim < NDIM; dim++) {
		p.v[dim] = 0.0;
	}
}

void collapse(particle &p) {
	if (abs(p.x) < 0.4) {
		p.m = 1.0e+6 * p.V;
	} else {
		p.m = 1.0 * p.V;
	}
	p.v = vect(0);
}

init_func_type get_initialization_function(const std::string &name) {
	if (name == "drift") {
		return drift;
#if(NDIM!=1)
	} else if (name == "kh") {
		return kh;
#endif
	} else if (name == "collapse") {
		return collapse;
	} else if (name == "blast") {
		return blast;
	} else if (name == "sod") {
		return sod;
	} else if (name == "kepler") {
		return kepler;
	} else {
		printf("Error: Initialization function %s is not known\n", name.c_str());
		abort();
	}
}
