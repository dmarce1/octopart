/*
 * state.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */

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
	F.momentum() = ::rotate_from(F.momentum(), norm);

}

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

real primitive_state::sound_speed() const {
	return std::sqrt(FGAMMA * std::max(pre(), 0.0) / den());
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

flux_state riemann_solver(const primitive_state &VL, const primitive_state &VR) {
	flux_state F, FR, FL;
	const auto UL = VL.to_con();
	const auto UR = VR.to_con();
	const auto &rhoL = UL.den();
	const auto &rhoR = UR.den();
	const auto &EL = UL.ene();
	const auto &ER = UR.ene();
	const auto &PL = std::max(VL.pre(), 0.0);
	const auto &PR = std::max(VR.pre(), 0.0);
	const auto hL = (PL + EL) / rhoL;
	const auto hR = (PR + ER) / rhoR;
	const auto uL = VL.vel()[0];
	const auto uR = VR.vel()[0];
	const auto wR = std::sqrt(rhoR);
	const auto wL = std::sqrt(rhoL);
	const auto uroe = (VR.vel() * wR + VL.vel() * wL) / (wR + wL);
	const auto hroe = (hR * wR + hL * wL) / (wR + wL);
	const auto aroe = std::sqrt(std::max(0.0, (FGAMMA - 1.0) * (hroe - uroe.dot(uroe) / 2.0)));
	const auto sR = uroe[0] + aroe;
	const auto sL = uroe[0] - aroe;

	if (sR > 0.0 && sL < 0.0) {
		const auto sstar_den = rhoL * (sL - uL) - rhoR * (sR - uR);
		const auto sstar_num = PR - PL + rhoL * uL * (sL - uL) - rhoR * uR * (sR - uR);
		const auto sstar = std::max(sL, std::min(sR, sstar_num / sstar_den));
		if (sstar > 0.0) {

		} else {

		}
	} else {
		if (sR < 0.0 && sL < 0.0) {

		} else {

		}
	}
	return F;
}
