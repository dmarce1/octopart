
#include <octopart/math.hpp>
#include <octopart/state.hpp>

flux_state primitive_state::to_flux() const {
	flux_state f;
	const auto v = (*this)[v_i];
	f[d_i] = den() * v;
	for (int dim = 0; dim < NDIM; dim++) {
		f[v_i] = den() * v * (*this)[v_i + dim];
	}
	const auto P = max(real(0.0), pre());
	f[v_i] += P;
	f[p_i] = (pre() / (FGAMMA - 1.0) + P + den() * vel().dot(vel()) * 0.5) * v;
	return f;
}

flux_state flux_state::boost_from(const vect &vf) const {
	flux_state F;
	F.mass() = mass();
	F.energy() = energy() + vf.dot(vf) * 0.5 * mass() + vf.dot(momentum());
	F.momentum() = momentum() + vf * mass();
	return F;
}

flux_state flux_state::rotate_from(const vect &norm) const {
	flux_state F;
	F.mass() = mass();
	F.energy() = energy();
	F.momentum() = ::rotate_from(momentum(), norm);
	return F;
}

