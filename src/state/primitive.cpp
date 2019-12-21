
#include <octopart/state.hpp>
#include <octopart/math.hpp>
#include <octopart/options.hpp>



flux_state primitive_state::to_flux() const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	flux_state f;
	const auto v = (*this)[v_i];
	f[d_i] = den() * v;
	f[t_i] = den() * v;
	for (int dim = 0; dim < NDIM; dim++) {
		f[v_i + dim] = den() * v * (*this)[v_i + dim];
	}
	const auto P = max(real(0.0), pre());
	f[v_i] += P;
	f[p_i] = (pre() / (fgamma - 1.0) + P + den() * vel().dot(vel()) * 0.5) * v;
	return f;
}

conserved_state primitive_state::to_con() const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	conserved_state U;
	U.den() = den();
	U.ene() = pre() / (fgamma - 1.0) + vel().dot(vel()) * den() * 0.5;
	U.mom() = vel() * den();
	U.tau() = tau();
	return U;
}

real primitive_state::sound_speed() const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	return sqrt(fgamma * max(pre(), real(0.0)) / den());
}

primitive_state primitive_state::boost_to(const vect &vf) const {
	primitive_state V = *this;
	auto &v = V.vel();
	v = v - vf;
	return V;
}

primitive_state primitive_state::dW_dt(const gradient &dW_dx) const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	primitive_state V;
	for (int i = 0; i < STATE_SIZE; i++) {
		V[i] = 0;
	}
	for (int dim = 0; dim < NDIM; dim++) {
		const auto u = (*this)[v_i + dim];
		V[d_i] -= u * dW_dx[dim][d_i];
		V[d_i] -= (*this)[d_i] * dW_dx[dim][v_i + dim];
		V[t_i] -= u * dW_dx[dim][t_i];
		V[t_i] -= (*this)[t_i] * dW_dx[dim][v_i + dim];
		V[p_i] -= u * dW_dx[dim][p_i];
		V[p_i] -= fgamma * (*this)[p_i] * dW_dx[dim][v_i + dim];
		for (int n = 0; n < NDIM; n++) {
			V[v_i + n] -= u * dW_dx[dim][v_i + n];
			V[v_i + n] -= dW_dx[dim][p_i] / (*this)[d_i];
		}
	}
	return V;
}

primitive_state primitive_state::rotate_to(const vect &norm) const {
	primitive_state V = *this;
	V.vel() = ::rotate_to(V.vel(), norm);
	return V;
}

