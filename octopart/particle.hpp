/*
 * particle.hpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#ifndef SRC_PARTICLE_HPP_
#define SRC_PARTICLE_HPP_

#include <vector>

#include "state.hpp"



struct particle {
	vect x;
	real m;
	vect u;
	real e;
	real V;
	real h;
	vect psi_a;
	void write(FILE*) const;
	int read(FILE*);
	template<class Arc>
	void serialize(Arc&& a, unsigned) {
		a & x;
		a & m;
		a & u;
		a & e;
		a & V;
		a & h;
		a & psi_a;
	}
	primitive_state to_prim() const;
	conserved_state to_con() const;
	particle from_con(const conserved_state&) const;
};



std::vector<particle> cartesian_particle_set(int);
std::vector<particle> random_particle_set(int);

#endif /* SRC_PARTICLE_HPP_ */
