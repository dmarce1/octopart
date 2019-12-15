#include <octopart/math.hpp>
#include <octopart/state.hpp>

flux_state KT(const primitive_state &VL, const primitive_state &VR) {
	flux_state F;
	const auto UR = VR.to_con();
	const auto UL = VL.to_con();
	const auto vR = VR.vel()[0];
	const auto vL = VL.vel()[0];
	const auto aR = VR.sound_speed() + abs(vR);
	const auto aL = VL.sound_speed() + abs(vL);
	const auto a = max(aR, aL);
	F = (VR.to_flux() + VL.to_flux() - (UR - UL) * a) * 0.5;
	return F;
}

static flux_state HLLC(const primitive_state &VL, const primitive_state &VR) {
	flux_state F;
	flux_state FR;
	flux_state FL;
	const auto UL = VL.to_con();
	const auto UR = VR.to_con();
	const auto aL = VL.sound_speed();
	const auto aR = VR.sound_speed();
	const auto rhoL = UL.den();
	const auto rhoR = UR.den();
	const auto EL = UL.ene();
	const auto ER = UR.ene();
	const auto PL = VL.pre();
	const auto PR = VR.pre();
	const auto rho_bar = 0.5 * (rhoR + rhoL);
	const auto a_bar = (aL + aR) * 0.5;
	const auto uL = VL.vel()[0];
	const auto uR = VR.vel()[0];
	const auto P0 = std::max(real(0), 0.5 * (PL + PR) - 0.5 * (uR - uL) * rho_bar * a_bar);
	if (PL <= 0.0 || PR <= 0.0 || P0 <= 0.0) {
		F = KT(VL, VR);
	} else {
		real qR, qL;
		if (P0 < PR) {
			qR = 1.0;
		} else {
			qR = sqrt(1.0 + (FGAMMA + 1.0) / (2.0 * FGAMMA) * (P0 / PR - 1.0));
		}
		if (P0 < PL) {
			qL = 1.0;
		} else {
			qL = sqrt(1.0 + (FGAMMA + 1.0) / (2.0 * FGAMMA) * (P0 / PL - 1.0));
		}
		const auto sR = uR + aR * qR;
		const auto sL = uL - aL * qL;
		if (sL > 0.0) {
			F = VL.to_flux();
		} else if (sR < 0.0) {
			F = VR.to_flux();
		} else {
			const auto s0_den = rhoL * (sL - uL) - rhoR * (sR - uR);
			const auto s0_num = PR - PL + rhoL * uL * (sL - uL) - rhoR * uR * (sR - uR);
			const auto s0 = s0_num / s0_den;
			conserved_state U0;
			if (s0 > 0.0) {
				const auto rho0 = rhoL * (sL - uL) / (sL - s0);
				U0.den() = rho0;
				U0.mom()[0] = rho0 * s0;
				for (int dim = 1; dim < NDIM; dim++) {
					U0.mom()[dim] = rho0 * VL.vel()[dim];
				}
				U0.ene() = rho0 * (EL / rhoL + (s0 - uL) * (s0 + PL / (rhoL * (sL - uL))));
				const auto FL = VL.to_flux();
				F = FL + (U0 - UL) * sL;
			} else {
				const auto rho0 = rhoR * (sR - uR) / (sR - s0);
				U0.den() = rho0;
				U0.mom()[0] = rho0 * s0;
				for (int dim = 1; dim < NDIM; dim++) {
					U0.mom()[dim] = rho0 * VR.vel()[dim];
				}
				U0.ene() = rho0 * (ER / rhoR + (s0 - uR) * (s0 + PR / (rhoR * (sR - uR))));
				const auto FR = VR.to_flux();
				F = FR + (U0 - UR) * sR;
			}
		}
		return F;
	}
}


flux_state riemann_solver(const primitive_state &VL, const primitive_state &VR) {
	return HLLC(VL, VR);
}
