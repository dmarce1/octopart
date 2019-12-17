#include <octopart/initialize.hpp>
#include <octopart/rand.hpp>

void drift(particle &p) {
	p.m = 1.0;
	p.e = 0.0;
	const auto r = abs(p.x);
	p.u = vect(0);
	p.u[0] = rand_unit_box();
	p.u[1] = rand_unit_box();

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
	p.u[0] = -vx;
	p.u[1] = vy;
#if(NDIM==3)
	p.u[2] = 0.0;
#endif
	p.m = rho * p.V;
	p.e = E * p.V + p.u.dot(p.u) * 0.5 * p.m;
	;
//	for (int dim = 0; dim < NDIM; dim++) {
//		p.u[dim] = rand_unit_box() * 0.01;
//	}
//	if (abs(p.x[1]) > 0.25) {
//		p.u[0] += +0.5;
//		p.m = p.V;
//		p.e = 6.25 * p.V + p.u.dot(p.u) * 0.5 * p.m;
//	} else {
//		p.u[0] += -0.5;
//		p.m = 2.0 * p.V;
	//		p.e = 6.25 * p.V + p.u.dot(p.u) * 0.5 * p.m; ;
//	}
}

void kepler(particle &p) {
	const auto r = abs(p.x);
	const auto y = p.x[1];
	const auto x = p.x[0];
	p.u = vect(0);
	p.u[0] = -y / r / sqrt(r);
	p.u[1] = +x / r / sqrt(r);
	if (r < 1.0 / 12.0) {
		p.m = 0.01 + 12 * r * r * r;
	} else if (r < 1.0 / 3.0) {
		p.m = 0.01 + 1.0;
	} else {
		p.m = 0.01 + 1.0 / pow(1.0 + (r - 1.0 / 3.0) / .1, 3.0);
	}
	p.m *= p.V;
	p.e = 1.0e-6 * p.V + 0.5 * p.u.dot(p.u) * p.m;
}

void sod(particle &p) {
	if (p.x[0] < 0.0) {
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
#if(NDIM!=1)
	} else if (name == "kh") {
		return kh;
#endif
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
