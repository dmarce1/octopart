#include <octopart/math.hpp>
#include <octopart/state.hpp>


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
	if (sL > 0.0) {
		F = VL.to_flux();
	} else if (sR < 0.0) {
		F = VR.to_flux();
	} else {
		const auto s0_den = rhoL * (sL - uL) - rhoR * (sR - uR);
		const auto s0_num = PR - PL + rhoL * uL * (sL - uL) - rhoR * uR * (sR - uR);
		const auto s0 = std::max(sL, std::min(sR, s0_num / s0_den));
		conserved_state U0;
		if (s0 > 0.0) {
			const auto rho0 = rhoL * (sL - uL) / (sL - s0);
			U0.den() = rho0;
			U0.mom()[0] = rho0 * s0;
			for (int dim = 1; dim < NDIM; dim++) {
				U0.mom()[dim] = rho0 * VL.vel()[dim];
			}
			U0.ene() = EL / rhoL + (s0 - uL) * (s0 + PL / (rhoL * (sL - uL)));
			F = VL.to_flux() + (U0 - UL) * s0;
		} else {
			const auto rho0 = rhoR * (sR - uR) / (sR - s0);
			U0.den() = rho0;
			U0.mom()[0] = rho0 * s0;
			for (int dim = 1; dim < NDIM; dim++) {
				U0.mom()[dim] = rho0 * VR.vel()[dim];
			}
			U0.ene() = ER / rhoR + (s0 - uR) * (s0 + PR / (rhoR * (sR - uR)));
			F = VR.to_flux() + (U0 - UR) * s0;
		}
	}
	return F;
}
