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

struct qcon_state {
	real m;
	real E;
	general_vect<real, NDIM> p;
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & m;
		arc & E;
		arc & p;
	}
};

class primitive_state;
using gradient = general_vect<primitive_state,NDIM>;

class flux_state: public state {
	static constexpr int mas_i = 0;
	static constexpr int ene_i = 1;
	static constexpr int mom_i = 2;
public:
	flux_state() = default;
	flux_state(const general_vect<real, STATE_SIZE> &other) {
		for (int i = 0; i < STATE_SIZE; i++) {
			(*this)[i] = other[i];
		}
	}
	flux_state& operator=(const general_vect<real, STATE_SIZE> &other) {
		for (int i = 0; i < STATE_SIZE; i++) {
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
		for (int i = 0; i < STATE_SIZE; i++) {
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

struct primitive_state {
	real rho;
	real p;
	vect v;
	real sound_speed() const;
	conserved_state to_con() const;
	primitive_state boost_to(const vect &v) const;
	primitive_state rotate_to(const vect &norm) const;
	primitive_state dW_dt(const gradient&) const;
	flux_state to_flux() const;
	void zero();
	primitive_state operator*(const real &other) const;
	primitive_state operator+(const primitive_state &other) const;
	primitive_state operator-() const;
	inline primitive_state operator/(const real &other) const {
		return (*this) * (1.0 / other);
	}
	inline primitive_state operator-(const primitive_state &other) const {
		return (*this) + -other;
	}
	real operator[](int i) const;
	real& operator[](int i);
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & rho;
		a & p;
		a & v;
	}
};

primitive_state max(const primitive_state &a, const primitive_state &b);
primitive_state min(const primitive_state &a, const primitive_state &b);
primitive_state abs(const primitive_state &a, const primitive_state &b);


flux_state KT(const primitive_state &UL, const primitive_state &UR);
bool riemann_solver(flux_state&, const primitive_state &UL, const primitive_state &UR);
bool exact_Riemann(flux_state &F, const primitive_state &VL, const primitive_state &VR);

#endif /* SRC_STATE_HPP_ */
