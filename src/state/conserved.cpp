#include <octopart/options.hpp>
#include <octopart/state.hpp>

primitive_state conserved_state::to_prim() const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	primitive_state V;
	V.den() = den();
	auto ei = (ene() - mom().dot(mom()) / den() * 0.5);
	V.pre() = (fgamma - 1.0) * ei;
	V.vel() = mom() / den();
	return V;
}
