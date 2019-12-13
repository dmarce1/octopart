
#include <octopart/math.hpp>
#include <octopart/state.hpp>

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

