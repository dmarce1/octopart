/*
 * state.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */

#include <octopart/state.hpp>

primitive_state conserved_state::to_prim() const {
	primitive_state V;
	V.den() = den();
	V.pre() = (FGAMMA - 1.0) * (ene() - mom().dot(mom()) / den() * 0.5);
	V.vel() = mom() / den();
	return V;
}

conserved_state primitive_state::to_con() const {
	conserved_state U;
	U.den() = den();
	U.ene() = pre() / (FGAMMA - 1.0) + vel().dot(vel()) * den() * 0.5;
	U.mom() = vel() * den();
	return U;
}

