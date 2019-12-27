#include <octopart/math.hpp>
#include <octopart/options.hpp>
#include <octopart/state.hpp>
#include <cassert>

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

bool exact_Riemann(flux_state &F, const primitive_state &VL, const primitive_state &VR) {

	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	const auto rhoL = VL.den();
	const auto rhoR = VR.den();
	real PL, PR;
	const auto UR = VR.to_con();
	const auto UL = VL.to_con();
	PR = max(VR.pre(), 1.0e-10);
	PL = max(VL.pre(), 1.0e-10);
	flux_state FR, FL;

	const auto gam1 = fgamma - 1.0;
	const auto gam2 = fgamma + 1.0;
	const auto gam3 = gam1 / (2.0 * fgamma);
	const auto gam4 = gam2 / (2.0 * fgamma);

	const auto AL = 2.0 / (gam2 * rhoL);
	const auto AR = 2.0 / (gam2 * rhoR);
	const auto BL = gam1 / gam2 * PL;
	const auto BR = gam1 / gam2 * PR;
	const auto uL = VL.vel()[0];
	const auto uR = VR.vel()[0];
	const auto aR = sqrt(fgamma * PR / rhoR);
	const auto aL = sqrt(fgamma * PL / rhoL);

	const auto sOL = uL + 2.0 * aL / gam1;
	const auto sOR = uR - 2.0 * aR / gam1;

	real dP, P0, err, s0;

	const auto fL = [&](real P) {
		if (P > PL) {
			return (P - PL) * sqrt(AL / (P + BL));
		} else {
			if (PL == 0.0) {
				return real(0.0);
			} else {
				return 2.0 * aL / gam1 * (pow(P / PL, gam3) - 1);
			}
		}
	};

	const auto fR = [&](real P) {
		if (P > PR) {
			return (P - PR) * sqrt(AR / (P + BR));
		} else {
			return 2.0 * aR / gam1 * (pow(P / PR, gam3) - 1);
		}
	};

	const auto dfLdP = [&](real P) {
		if (P > PL) {
			return sqrt(AL / (BL + P)) * (1.0 - (P - PL) / (2.0 * (BL + P)));
		} else {
			return 1.0 / (rhoL * aL) * pow(PL / P, gam4);
		}
	};

	const auto dfRdP = [&](real P) {
		if (P > PR) {
			return sqrt(AR / (BR + P)) * (1.0 - (P - PR) / (2.0 * (BR + P)));
		} else {
			return 1.0 / (rhoR * aR) * pow(PR / P, gam4);
		}
	};

	const auto f = [&](real P) {
		return fL(P) + fR(P) + uR - uL;
	};

	const auto dfdP = [&](real P) {
		return dfLdP(P) + dfRdP(P);
	};

	const auto VrL = [=]() {
		primitive_state V;
		const auto tmp = 2.0 / gam2 + (gam1 / gam2 / aL * uL);
		V.den() = rhoL * pow(tmp, 2.0 / gam1);
		V.vel()[0] = 2.0 / gam2 * (aL + gam1 / 2.0 * uL);
		V.pre() = PL * pow(2.0 / gam2 + (gam1 / gam2 / aL * uL), 2.0 * fgamma / gam1);
		for (int dim = 1; dim < NDIM; dim++) {
			V.vel()[dim] = VL.vel()[dim];
		}
		return V;
	};

	const auto VrR = [=]() {
		primitive_state V;
		const auto tmp = 2.0 / gam2 - (gam1 / gam2 / aR * uR);
		V.den() = rhoR * pow(tmp, 2.0 / gam1);
		V.vel()[0] = 2.0 / gam2 * (-aR + gam1 / 2.0 * uR);
		V.pre() = PR * pow(2.0 / gam2 - (gam1 / gam2 / aR * uR), 2.0 * fgamma / gam1);
		for (int dim = 1; dim < NDIM; dim++) {
			V.vel()[dim] = VR.vel()[dim];
		}
		return V;
	};

	const auto VO = [=]() {
		primitive_state V;
		for (int i = 0; i < STATE_SIZE; i++) {
			V[i] = 0.0;
		}
		return V;
	};

	const auto V0L = [&]() {
		primitive_state V;
		real rho0L;
		if (P0 > PL) {
			const auto num = P0 / PL + gam1 / gam2;
			const auto den = (gam1 / gam2) * P0 / PL + 1.0;
			rho0L = rhoL * num / den;
		} else {
			rho0L = rhoL * pow(P0 / PL, 1.0 / fgamma);
		}
		V.pre() = P0;
		V.den() = rho0L;
		V.vel()[0] = s0;
		for (int dim = 1; dim < NDIM; dim++) {
			V.vel()[dim] = VL.vel()[dim];
		}
		return V;
	};

	const auto V0R = [&]() {
		primitive_state V;
		real rho0R;
		if (P0 > PR) {
			const auto num = P0 / PR + gam1 / gam2;
			const auto den = (gam1 / gam2) * P0 / PR + 1.0;
			rho0R = rhoR * num / den;
		} else {
			rho0R = rhoR * pow(P0 / PR, 1.0 / fgamma);
		}
		V.pre() = P0;
		V.den() = rho0R;
		V.vel()[0] = s0;
		for (int dim = 1; dim < NDIM; dim++) {
			V.vel()[dim] = VR.vel()[dim];
		}
		return V;
	};

	primitive_state Vi;
	if (sOL < sOR) {
		if (sOL > 0.0) {
			Vi = VrL();
		} else if (sOR < 0.0) {
			Vi = VrR();
		} else {
			Vi = VO();
		}
		F = Vi.to_flux();
	} else {

		P0 = (PL + PR) / 2.0;
		if (P0 == 0.0) {
			P0 = 1.0e-6;
		}
		dP = (PR - PL) / 2.0;
		int iter = 0;
		real Plast;
		do {
			const auto g = f(P0);
			if (g == 0.0) {
				break;
			}
			const auto dgdp = dfdP(P0);
			assert(dgdp != 0.0);
			dP = -g / dgdp;
			real c0 = pow(0.1, double(iter) / double(1000));
			P0 = min(max(P0 + c0 * dP, P0 / 2.0), P0 * 2.0);
			if (iter > 1) {
				err = abs(P0 - Plast) / (P0 + Plast) / 2.0;
			} else {
				err = real::max();
			}
			iter++;
			//	printf("%i %e %e %e %e %e %e %e %e\n", iter, c0, PL, PR, uL, uR, P0, dP, err);
			if (iter >= 1000) {
				printf("Riemann solver failed to converge\n");
				abort();
				;
			}
			Plast = P0;
		} while (err > 1.0e-6 && P0 > real::min());

		s0 = 0.5 * (uL + uR) + 0.5 * (fR(P0) - fL(P0));
		const auto QL = sqrt((P0 + BL) / AL);
		const auto QR = sqrt((P0 + BR) / AR);
		real sL, sR, sHL, sTL, sHR, sTR;

		if (s0 >= 0.0) {
			if (P0 > PL) {
				sL = uL - QL / rhoL;
				if (sL > 0.0) {
					Vi = VL;
				} else {
					Vi = V0L();
				}
			} else {
				const auto a0L = aL * pow(P0 / PL, (fgamma - 1) / (2 * fgamma));
				sHL = uL - aL;
				sTL = s0 - a0L;
				assert(sHL <= sTL);
				if (sHL > 0.0) {
					Vi = VL;
				} else if (sTL < 0.0) {
					Vi = V0L();
				} else {
					Vi = VrL();
				}
			}
			FL = Vi.to_flux();
		}
		if (s0 <= 0.0) {
			if (P0 > PR) {
				sR = uR + QR / rhoR;
				if (sR < 0.0) {
					Vi = VR;
				} else {
					Vi = V0R();
				}
			} else {
				const auto a0R = aR * pow(P0 / PR, (fgamma - 1) / (2 * fgamma));
				sHR = uR + aR;
				sTR = s0 + a0R;
				assert(sHR >= sTR);
				if (sHR < 0.0) {
					Vi = VR;
				} else if (sTR > 0.0) {
					Vi = V0R();
				} else {
					Vi = VrR();
				}
			}
			FR = Vi.to_flux();
		}
		if (s0 > 0.0) {
			F = FL;
		} else if (s0 < 0.0) {
			F = FR;
		} else {
			F = (FR + FL) / 2.0;
		}
	}
	return true;
}

static bool HLLC(flux_state &F, const primitive_state &VL, const primitive_state &VR) {
	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	flux_state FR, F0R;
	flux_state FL, F0L;
	const auto UL = VL.to_con();
	const auto UR = VR.to_con();
	const auto aL = VL.sound_speed();
	const auto aR = VR.sound_speed();
	const auto rhoL = UL.den();
	const auto rhoR = UR.den();
	const auto EL = UL.ene();
	const auto ER = UR.ene();
	const auto PL = max(VL.pre(), 1.0e-10);
	const auto PR = max(VR.pre(), 1.0e-10);
	const auto uL = VL.vel()[0];
	const auto uR = VR.vel()[0];
	int iter = 0;
	real sR;
	real sL;
	const auto wR = sqrt(rhoR);
	const auto wL = sqrt(rhoL);
	const auto hR = ER / rhoR;
	const auto hL = EL / rhoL;
	const auto ubar = (VL.vel() * wL + VR.vel() * wR) / (wL + wR);
	const auto hbar = (hR * wR + hL * wL) / (wL + wR);
	const auto abar = sqrt(max((fgamma - 1.0) * (hbar - ubar.dot(ubar) / 2.0), real(0.0)));

	sR = ubar[0] + abar;
	sL = ubar[0] - abar;

	if (sL == 0.0 && sR == 0.0) {
		F = (VL.to_flux() + VR.to_flux()) / 2.0;
	} else if (sL >= 0.0) {
		F = VL.to_flux();
	} else if (sR <= 0.0) {
		F = VR.to_flux();
	} else {
		const auto s0_den = rhoL * (sL - uL) - rhoR * (sR - uR);
		const auto s0_num = PR - PL + rhoL * uL * (sL - uL) - rhoR * uR * (sR - uR);
		const auto s0 = s0_num / s0_den;
		auto P0 = 0.5 * (PR + rhoR * (sR - uR) * (s0 - uR) + PL + rhoL * (sL - uL) * (s0 - uL));
		if (P0 <= 0.0) {
			sR = max(uR, uL) + max(aL, aR);
			sL = min(uR, uL) - max(aL, aR);
			if (sL == 0.0 && sR == 0.0) {
				F = (VL.to_flux() + VR.to_flux()) / 2.0;
				goto RETURN;
			} else if (sL >= 0.0) {
				F = VL.to_flux();
				goto RETURN;
			} else if (sR <= 0.0) {
				F = VR.to_flux();
				goto RETURN;
			}
			P0 = 0.5 * (PR + rhoR * (sR - uR) * (s0 - uR) + PL + rhoL * (sL - uL) * (s0 - uL));
			if (P0 <= 0.0) {
				sR = max(uR + aR, uL + aL);
				sL = min(uR - aR, uL - aL);
				if (sL == 0.0 && sR == 0.0) {
					F = (VL.to_flux() + VR.to_flux()) / 2.0;
					goto RETURN;
				} else if (sL >= 0.0) {
					F = VL.to_flux();
					goto RETURN;
				} else if (sR <= 0.0) {
					F = VR.to_flux();
					goto RETURN;
				}
				P0 = 0.5 * (PR + rhoR * (sR - uR) * (s0 - uR) + PL + rhoL * (sL - uL) * (s0 - uL));
				if (P0 <= 0.0) {
					return exact_Riemann(F, VL, VR);
				}
			}
		}
		conserved_state U0R, U0L;
		if (s0 >= 0.0) {
			const auto rho0L = rhoL * (sL - uL) / (sL - s0);
			U0L.den() = rho0L;
			U0L.mom()[0] = rho0L * s0;
			for (int dim = 1; dim < NDIM; dim++) {
				U0L.mom()[dim] = rho0L * VL.vel()[dim];
			}
			if (rho0L != 0.0) {
				U0L.ene() = rho0L * (EL / rhoL + (s0 - uL) * (s0 + PL / (rhoL * (sL - uL))));
			} else {
				U0L.ene() = 0.0;
			}
			const auto FL = VL.to_flux();
			F0L = FL + (U0L - UL) * sL;
		}
		if (s0 <= 0.0) {
			const auto rho0R = rhoR * (sR - uR) / (sR - s0);
			U0R.den() = rho0R;
			U0R.mom()[0] = rho0R * s0;
			for (int dim = 1; dim < NDIM; dim++) {
				U0R.mom()[dim] = rho0R * VR.vel()[dim];
			}
			if (rho0R != 0.0) {
				U0R.ene() = rho0R * (ER / rhoR + (s0 - uR) * (s0 + PR / (rhoR * (sR - uR))));
			} else {
				U0R.ene() = 0.0;
			}
			const auto FR = VR.to_flux();
			F0R = FR + (U0R - UR) * sR;
		}
		if (s0 == 0.0) {
			F = (F0R + F0L) * 0.5;
		} else if (s0 < 0.0) {
			F = F0R;
		} else {
			F = F0L;
		}

	}
	RETURN: return true;
}

bool riemann_solver(flux_state &f, const primitive_state &VL, const primitive_state &VR) {
	return HLLC(f, VL, VR);
}
