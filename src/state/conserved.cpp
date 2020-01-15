#include <octopart/options.hpp>
#include <octopart/state.hpp>

primitive_state conserved_state::to_prim() const {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	primitive_state W;
	W.rho = den();
	auto ei = (ene() - mom().dot(mom()) / den() * 0.5);
	W.p = (fgamma - 1.0) * ei;
	W.v = mom() / den();
	return W;
}
