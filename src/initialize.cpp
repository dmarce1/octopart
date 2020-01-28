#include <octopart/initialize.hpp>
#include <octopart/options.hpp>
#include <octopart/physcon.hpp>
#include <octopart/polytrope.hpp>
#include <octopart/rand.hpp>

void drift(particle &p) {
	printf("Must re-write %s\n", __FUNCTION__);
//	p.Q.m = 1.0;
//	p.E = 0.0;
//	const auto r = abs(p.x);
//	p.v = vect(0);
//	p.v[0] = rand_unit_box();
//	p.v[1] = rand_unit_box();

}

void rt(particle &p) {
	printf("Must re-write %s\n", __FUNCTION__);
//	const auto opts = options::get();
//	const auto fgamma = opts.fgamma;
//	const auto rho1 = 1.0;
//	const auto rho2 = 2.0;
//	auto y = p.x[1];
//	const auto delta = 0.025;
//	const auto rho = rho1 + (rho2 - rho1) / (1 + exp(-y / delta));
//	y = 0.5 - y;
//	const auto P = 0.5 * (rho1 * y + (rho2 - rho1) * delta * log(1 + exp(y / delta)));
//
//	p.v = vect(0);
//	if (y > -0.2 && y < 0.2) {
//		p.v[1] = 0.025 * (1.0 + cos(8 * M_PI * (p.x[0] + 0.25))) * (1.0 + cos(5 * M_PI * y));
//	}
//
//	p.Q.m = rho * p.V;
//	p.E = P / (fgamma - 1) * p.V + p.v.dot(p.v) * p.Q.m / 2.0;
//
}

void kh(particle &p) {
//	printf("Must re-write %s\n", __FUNCTION__);
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
	p.Q.m = rho * p.V;
	p.Q.p[0] = -vx * p.Q.m;
	p.Q.p[1] = vy * p.Q.m;
#if(NDIM==3)
	p.Q.p[2] = 0.0;
#endif
	p.Q.E = E * p.V + p.Q.p.dot(p.Q.p) / (2.0 * p.Q.m);
}

void kepler(particle &p) {
#if(NDIM==1)
	printf("Cannot do Kepler problem in 1D\n");
	abort();
#else
	static const auto opts = options::get();
	static const auto eps = opts.kep_eps;
	const auto r = abs(p.x);
	const auto y = p.x[1];
	const auto x = p.x[0];
	p.Q.p = vect(0);
	const auto tmp = pow(eps * eps + r * r, -0.75);
	p.Q.p[0] = -y * tmp;
	p.Q.p[1] = +x * tmp;
	if (r < 1.0 / 12.0) {
		p.Q.m = 0.01 + 12 * r * r * r;
	} else if (r < 1.0 / 3.0) {
		p.Q.m = 0.01 + 1.0;
	} else {
		p.Q.m = 0.01 + 1.0 / pow(1.0 + (r - 1.0 / 3.0) / .1, 3.0);
	}
	p.Q.m *= p.V;
	p.Q.p = p.Q.p * p.Q.m;
	p.Q.E = 1.0e-6 * p.V + 0.5 * p.Q.p.dot(p.Q.p) / p.Q.m;
#endif
}

void sod(particle &p) {
	for (int dim = 0; dim < NDIM; dim++) {
		p.Q.p[dim] = 0.0;
	}
	if (p.x[0] < 0) {
		p.Q.m = 1.0 * p.V;
		p.Q.E = 2.5 * p.V;
	} else {
		p.Q.m = 0.125 * p.V;
		p.Q.E = 0.25 * p.V;
	}
	p.Q.E = p.Q.E + p.Q.p.dot(p.Q.p) / 2.0 / p.Q.m;
}

void blast(particle &p) {
	const auto r = sqrt(p.x.dot(p.x));
	p.Q.m = p.V;
	for (int dim = 0; dim < NDIM; dim++) {
		p.Q.p[dim] = 0.0;
	}
	p.Q.E = max( 1.0e-10, p.V * exp(-r * r * 10000));
}

void single_polytrope(particle &p) {
	static const auto opts = options::get();
	const auto r = abs(p.x);
	const auto rho = max(1.0e-6, polytrope(r));
	const auto c0 = 4.0 * M_PI * physcon::G;
	p.Q.m = rho * p.V;
	p.Q.E = c0 * pow(rho, opts.fgamma) * p.V;
	p.Q.p = vect(0);
}

void collapse(particle &p) {
	if (abs(p.x) < 0.4) {
		p.Q.m = 1.0e+6 * p.V;
	} else {
		p.Q.m = 1.0 * p.V;
	}
	p.Q.p = vect(0);
	p.Q.E = p.Q.m * 1.0e-6;
}

init_func_type get_initialization_function(const std::string &name) {
	init_func_type init;
	if (name == "drift") {
		init = drift;
#if(NDIM!=1)
	} else if (name == "kh") {
		init =kh;
#endif
	} else if (name == "collapse") {
		init = collapse;
	} else if (name == "blast") {
		init = blast;
	} else if (name == "sod") {
		init = sod;
	} else if (name == "kepler") {
		init = kepler;
	} else if (name == "rt") {
		init = rt;
	} else if (name == "polytrope") {
		init = single_polytrope;
	} else {
		printf("Error: Initialization function %s is not known\n", name.c_str());
		abort();
	}
	return [init](particle &p) {
		init(p);
		p.t = p.dt = fixed_real(0);
		p.con_to_prim();
	};
}
