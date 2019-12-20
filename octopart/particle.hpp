/*
 * particle.hpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#ifndef SRC_PARTICLE_HPP_
#define SRC_PARTICLE_HPP_

#include <octopart/containers.hpp>
#include <octopart/vect.hpp>
#include <vector>

#include "state.hpp"

struct particle {
	vect x;
	real m;
	vect u;
	vect vf;
	vect g;
	real E;
	real U;
	real V;
	real h;
	array<vect, NDIM> B;
	void write(FILE*) const;
	int read(FILE*);
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & vf;
		a & g;
		a & x;
		a & m;
		a & u;
		a & E;
		a & U;
		a & V;
		a & h;
		a & B;
	}
	primitive_state to_prim() const;
	conserved_state to_con() const;
	particle from_con(const conserved_state&) const;
};

vector<particle> cartesian_particle_set(int);
vector<particle> random_particle_set(int);
vector<particle> disc_particle_set( int N);

#endif /* SRC_PARTICLE_HPP_ */
