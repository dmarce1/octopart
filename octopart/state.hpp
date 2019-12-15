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


using state = general_vect<real, STATE_SIZE>;

class primitive_state;
using gradient = general_vect<primitive_state,NDIM>;


class flux_state: public state {
	static constexpr int mas_i = 0;
	static constexpr int ene_i = 1;
	static constexpr int mom_i = 2;
public:
	flux_state& operator=(const general_vect<real, STATE_SIZE> &other) {
		for( int i = 0; i < STATE_SIZE; i++) {
			(*this)[i] = other[i];
		}
		return *this;
	}
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
	flux_state boost_from(const vect &v) const;
	flux_state rotate_from(const vect &norm) const;
};

class conserved_state: public state {
	static constexpr int d_i = 0;
	static constexpr int e_i = 1;
	static constexpr int s_i = 2;
public:
	conserved_state& operator=(const general_vect<real, STATE_SIZE> &other) {
		for( int i = 0; i < STATE_SIZE; i++) {
			(*this)[i] = other[i];
		}
		return *this;
	}
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
	primitive_state& operator=(const general_vect<real, STATE_SIZE> &other) {
		for( int i = 0; i < STATE_SIZE; i++) {
			(*this)[i] = other[i];
		}
		return *this;
	}
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
	real sound_speed() const;
	conserved_state to_con() const;
	primitive_state boost_to(const vect &v) const;
	primitive_state rotate_to(const vect &norm) const;
	primitive_state dW_dt(const gradient&) const;
	flux_state to_flux() const;
};

flux_state riemann_solver(const primitive_state &UL, const primitive_state &UR);

#endif /* SRC_STATE_HPP_ */
