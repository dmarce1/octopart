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

#endif /* SRC_STATE_HPP_ */
