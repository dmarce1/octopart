#include <octopart/state.hpp>
#include <octopart/math.hpp>
#include <octopart/options.hpp>

flux_state primitive_state::to_flux() const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	flux_state f;
	const auto vx = v[0];
	f[0] = rho * vx;
	for (int dim = 0; dim < NDIM; dim++) {
		f[2 + dim] = rho * vx * v[dim];
	}
	const auto P = max(real(0.0), p);
	f[2] += P;
	f[1] = (p / (fgamma - 1.0) + P + rho * v.dot(v) * 0.5) * vx;
	return f;
}

conserved_state primitive_state::to_con() const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	conserved_state U;
	U.den() = rho;
	U.ene() = p / (fgamma - 1.0) + v.dot(v) * rho * 0.5;
	U.mom() = v * rho;
	return U;
}

real primitive_state::sound_speed() const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	if (rho < 0.0) {
		printf("Negative density on recon %e\n", rho.get());
	}
	return sqrt(fgamma * max(p, real(0.0)) / rho);
}

primitive_state primitive_state::boost_to(const vect &vf) const {
	primitive_state W = *this;
	W.v = W.v - vf;
	return W;
}

primitive_state primitive_state::dW_dt(const gradient &dW_dx) const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	primitive_state dW;
	dW.zero();
	for (int dim = 0; dim < NDIM; dim++) {
		const auto u = v[dim];
		dW.rho -= u * dW_dx[dim].rho;
		dW.rho -= rho * dW_dx[dim].v[dim];
		dW.p -= u * dW_dx[dim].p;
		dW.p -= fgamma * p * dW_dx[dim].v[dim];
		for (int n = 0; n < NDIM; n++) {
			dW.v[n] -= u * dW_dx[dim].v[n];
		}
		dW.v[dim] -= dW_dx[dim].p / rho;
	}
	return dW;
}

primitive_state primitive_state::operator*(const real &other) const {
	primitive_state rc;
	rc.rho = rho * other;
	rc.p = p * other;
	rc.v = v * other;
	return rc;
}



real primitive_state::operator[](int i) const {
	switch (i) {
	case 0:
		return rho;
	case 1:
		return p;
	default:
		return v[i - 2];
	}
}

real& primitive_state::operator[](int i) {
	switch (i) {
	case 0:
		return rho;
	case 1:
		return p;
	default:
		return v[i - 2];
	}
}

primitive_state primitive_state::operator+(const primitive_state &other) const {
	primitive_state rc;
	rc.rho = rho + other.rho;
	rc.p = p + other.p;
	rc.v = v + other.v;
	return rc;
}

primitive_state abs(const primitive_state &a, const primitive_state &b) {
	primitive_state rc;
	rc.rho = abs(a.rho - b.rho);
	rc.p = abs(a.p - b.p);
	rc.v = abs(a.v - b.v);
	return rc;
}

primitive_state primitive_state::operator-() const {
	primitive_state rc;
	rc.rho = -rho;
	rc.p = -p;
	rc.v = -v;
	return rc;
}

primitive_state max(const primitive_state &a, const primitive_state &b) {
	primitive_state rc;
	rc.rho = max(a.rho, b.rho);
	rc.p = max(a.p, b.p);
	rc.v = max(a.v, b.v);
	return rc;
}

primitive_state min(const primitive_state &a, const primitive_state &b) {
	primitive_state rc;
	rc.rho = min(a.rho, b.rho);
	rc.p = min(a.p, b.p);
	rc.v = min(a.v, b.v);
	return rc;
}

primitive_state primitive_state::rotate_to(const vect &norm) const {
	primitive_state W = *this;
	W.v = ::rotate_to(W.v, norm);
	return W;
}

void primitive_state::zero() {
	v = vect(0.0);
	p = 0.0;
	rho = 0.0;
}

