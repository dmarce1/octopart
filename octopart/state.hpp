/*
 * state.hpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */

#ifndef SRC_STATE_HPP_
#define SRC_STATE_HPP_

#include "vect.hpp"

static constexpr int STATE_SIZE = NDIM + 2;
static constexpr real FGAMMA = 5.0 / 3.0;

using state = general_vect<real, STATE_SIZE>;

class primitive_state;


class flux_state : public state{
	static constexpr int mas_i = 0;
	static constexpr int ene_i = 1;
	static constexpr int mom_i = 2;
public:
	inline real& mass() {
		return (*this)[mas_i];
	}
	inline real& energy() {
		return (*this)[ene_i];
	}
	inline real mass() const {
		return (*this)[mas_i];
	}
	inline real energy() const {
		return (*this)[ene_i];
	}
	inline vect& momentum() {
		auto *ptr = reinterpret_cast<vect*>(&((*this)[mom_i]));
		return *ptr;
	}
	inline vect momentum() const {
		vect v;
		for (int dim = 0; dim < NDIM; dim++) {
			v[dim] = (*this)[mom_i + dim];
		}
		return v;
	}
	flux_state boost_from(const vect& v) const;
	flux_state rotate_from(const vect& norm) const;
};


class conserved_state: public state {
	static constexpr int d_i = 0;
	static constexpr int e_i = 1;
	static constexpr int s_i = 2;
public:
	inline real& den() {
		return (*this)[d_i];
	}
	inline real& ene() {
		return (*this)[e_i];
	}
	inline real den() const {
		return (*this)[d_i];
	}
	inline real ene() const {
		return (*this)[e_i];
	}
	inline vect& mom() {
		auto *ptr = reinterpret_cast<vect*>(&((*this)[s_i]));
		return *ptr;
	}
	inline vect mom() const {
		vect v;
		for (int dim = 0; dim < NDIM; dim++) {
			v[dim] = (*this)[s_i + dim];
		}
		return v;
	}
	primitive_state to_prim() const;
	flux_state to_flux() const;
	conserved_state boost_to(const vect& v) const;
	conserved_state rotate_to(const vect& norm) const;
};

class primitive_state: public state {
	static constexpr int d_i = 0;
	static constexpr int p_i = 1;
	static constexpr int v_i = 2;
public:
	inline real& den() {
		return (*this)[d_i];
	}
	inline real& pre() {
		return (*this)[p_i];
	}
	inline real den() const {
		return (*this)[d_i];
	}
	inline real pre() const {
		return (*this)[p_i];
	}
	inline vect& vel() {
		auto *ptr = reinterpret_cast<vect*>(&((*this)[v_i]));
		return *ptr;
	}
	inline vect vel() const {
		vect v;
		for (int dim = 0; dim < NDIM; dim++) {
			v[dim] = (*this)[v_i + dim];
		}
		return v;
	}
	conserved_state to_con() const;
};

flux_state riemann_solver(const conserved_state& UL, const conserved_state& UR);


#endif /* SRC_STATE_HPP_ */
