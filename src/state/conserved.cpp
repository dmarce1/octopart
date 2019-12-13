
#include <octopart/state.hpp>


primitive_state conserved_state::to_prim() const {
	primitive_state V;
	V.den() = den();
	V.pre() = (FGAMMA - 1.0) * (ene() - mom().dot(mom()) / den() * 0.5);
	V.vel() = mom() / den();
	return V;
}
