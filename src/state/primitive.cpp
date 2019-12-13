
#include <octopart/state.hpp>
#include <octopart/math.hpp>

conserved_state primitive_state::to_con() const {
	conserved_state U;
	U.den() = den();
	U.ene() = pre() / (FGAMMA - 1.0) + vel().dot(vel()) * den() * 0.5;
	U.mom() = vel() * den();
	return U;
}

real primitive_state::sound_speed() const {
	return sqrt(FGAMMA * max(pre(), real(0.0)) / den());
}

primitive_state primitive_state::boost_to(const vect &vf) const {
	primitive_state V = *this;
	auto &v = V.vel();
	v = v - vf;
	return V;
}

primitive_state primitive_state::dW_dt(const gradient &dW_dx) const {
	primitive_state V;
	for (int i = 0; i < STATE_SIZE; i++) {
		V[i] = 0;
	}
	for (int dim = 0; dim < NDIM; dim++) {
		const auto u = (*this)[v_i + dim];
		V[d_i] -= u * dW_dx[dim][d_i];
		V[d_i] -= (*this)[d_i] * dW_dx[dim][v_i + dim];
		V[p_i] -= u * dW_dx[dim][p_i];
		V[p_i] -= FGAMMA * (*this)[p_i] * dW_dx[dim][v_i + dim];
		for (int n = 0; n < NDIM; n++) {
			V[v_i + n] -= u * (*this)[v_i + n];
			V[v_i + n] -= u * dW_dx[dim][p_i] * (*this)[d_i + n];
		}
	}
	return V;
}

primitive_state primitive_state::rotate_to(const vect &norm) const {
	primitive_state V = *this;
	V.vel() = ::rotate_to(V.vel(), norm);
	return V;
}

