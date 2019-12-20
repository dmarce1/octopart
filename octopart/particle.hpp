/*
 * particle.hpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#ifndef SRC_PARTICLE_HPP_
#define SRC_PARTICLE_HPP_

#include <octopart/vect.hpp>
#include <vector>

#include "state.hpp"

struct particle {
	vect x;
	real m;
	vect u;
	vect vf;
	real e;
	real V;
	real h;
	std::array<vect, NDIM> B;
	void write(FILE*) const;
	int read(FILE*);
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & x;
		a & m;
		a & u;
		a & vf;
		a & e;
		a & V;
		a & h;
		a & B;
	}
	primitive_state to_prim() const;
	conserved_state to_con() const;
	particle from_con(const conserved_state&) const;
};

std::vector<particle> cartesian_particle_set(int);
std::vector<particle> random_particle_set(int);
std::vector<particle> disc_particle_set( int N);

#endif /* SRC_PARTICLE_HPP_ */
